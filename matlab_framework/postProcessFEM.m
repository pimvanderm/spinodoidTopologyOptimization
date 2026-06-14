function [results_summary] = postProcessFEM(tr, meshData, febioPaths, settings)
    % POSTPROCESSFEM Extracts and compiles specific multiscale validation metrics
    %
    % INPUTS:
    %   tr         - Original triangulation object (for exact volume)
    %   meshData   - Struct from runExplicitFEM containing V, E, elemVol, bcPrescribeList
    %   febioPaths - File paths for logs
    %   settings   - Must contain: .V_predicted, .c_TopOpt, .U_topopt_scaled
    %
    % OUTPUTS:
    %   results_summary - Structure containing DIC profiles, stress ratios, and displacements

    fprintf('--- POST PROCESSING: Harvesting Verification Metrics ---\n');
    
    % 1. EXACT GEOMETRIC VOLUME & RATIO
    v = tr.Points; f = tr.ConnectivityList;
    v1 = v(f(:,1), :); v2 = v(f(:,2), :); v3 = v(f(:,3), :);
    signed_volumes = sum(v1 .* cross(v2, v3, 2), 2) / 6;
    V_rendered = abs(sum(signed_volumes));
    
    results_summary.V_rendered = V_rendered;
    results_summary.V_ratio = V_rendered / settings.V_predicted;

    % 2. READ SOLVER LOGS
    dataStruct_disp = importFEBio_logfile(febioPaths.disp, 0, 1);
    N_disp_mat = dataStruct_disp.data;
    
    dataStruct_stress = importFEBio_logfile(febioPaths.stress, 0, 1);
    E_stress_mat = dataStruct_stress.data;
    E_stress_mat(isnan(E_stress_mat)) = 0;

    % 3. LOADING LINE VERTICAL (Z) DISPLACEMENTS
    % --- FEBio Mean Z-Displacement at Load Row ---
    febio_U_z_load = N_disp_mat(meshData.bcPrescribeList, 3, end);
    results_summary.mean_U_z_load_FEBio = mean(febio_U_z_load);
    
    % --- TopOpt Mean Z-Displacement at Load Row ---
    % Find the corresponding load line nodes on the TopOpt grid
    % Matching the loading line at x = 1 (1st node), z = nelz+1 (top node row)
    % =========================================================================
    % 3. LOADING LINE VERTICAL (Z) DISPLACEMENTS (Corrected Coordinate Traps)
    % =========================================================================
    nelx_to = settings.nelx;
    nely_to = settings.nely;
    nelz_to = settings.nelz;
    nnodz = nelz_to + 1;
    nnody = nely_to + 1;
    nnodx = nelx_to + 1;
    
    % Reconstruct the identical index grid topology used by your optimizer
    % Elements vary Z-fastest, then Y, then X: [Z, Y, X]
    [GridZ, GridY, GridX] = ndgrid(1:nnodz, 1:nnody, 1:nnodx);
    
    % Track if we are analyzing a half-symmetry model using robust strcmp
    if strcmp(settings.problemType, 'three_point_bending_half')
        % Target index specified by user: [z = 1, y = all, x = 2*nelz + 1]
        target_X_coord = 2 * nelz_to + 1;
        topopt_load_nodes = find(GridZ == 1 & GridX == target_X_coord);
        
        % The optimization load was applied UPWARDS (+Z). 
        % We take the absolute value or invert the sign to match FEBio's downward push.
        topopt_U_z_load = settings.U_topopt_scaled(3 * topopt_load_nodes);
        results_summary.mean_U_z_load_TopOpt = abs(mean(topopt_U_z_load));
        
    else
        % Full Three-Point Bending: [z = nelz + 1, y = all, x = nelx/2 + 1]
        target_X_coord = nelx_to / 2 + 1;
        topopt_load_nodes = find(GridZ == nnodz & GridX == target_X_coord);
        
        topopt_U_z_load = settings.U_topopt_scaled(3 * topopt_load_nodes);
        results_summary.mean_U_z_load_TopOpt = abs(mean(topopt_U_z_load));
    end
    
    % Harvest FEBio vertical displacement (always absolute to calculate true structural work)
    febio_U_z_load = N_disp_mat(meshData.bcPrescribeList, 3, end);
    results_summary.mean_U_z_load_FEBio = abs(mean(febio_U_z_load));

    % 4. GLOBAL MECHANICS: Displacement Max & Energy Mismatch
    results_summary.max_U = max(sqrt(sum(N_disp_mat(:,:,end).^2, 2)));
    
    loaded_U_z = N_disp_mat(meshData.bcPrescribeList, 3, end);
    nodalForce = -1 / length(meshData.bcPrescribeList); % Aligned with your 1N convergence runs
    c_Explicit = sum(loaded_U_z .* nodalForce);
    
    results_summary.c_Explicit = c_Explicit;
    results_summary.energy_mismatch = abs(c_Explicit - settings.c_TopOpt) / settings.c_TopOpt;

    % 5. CONSTITUTIVE STRESS METRICS (Mean, Peak, and Ratio)
    S_vm_ND = sqrt(0.5 * ((E_stress_mat(:,1,end) - E_stress_mat(:,2,end)).^2 + ...
                          (E_stress_mat(:,2,end) - E_stress_mat(:,3,end)).^2 + ...
                          (E_stress_mat(:,1,end) - E_stress_mat(:,3,end)).^2));
                      
    % --- Volume-Weighted Mean Stress ---
    totalVol = sum(meshData.elemVol);
    S_vm_mean = sum(S_vm_ND .* meshData.elemVol) / totalVol;
    results_summary.mean_stress = S_vm_mean;
    
    % --- 99th Percentile Peak Stress ---
    sorted_stress = sort(S_vm_ND);
    idx_99 = round(0.99 * length(sorted_stress));
    results_summary.stress_99th_pct = sorted_stress(idx_99);
    
    % --- Stress Concentration Ratio ---
    results_summary.stress_peak_to_mean_ratio = results_summary.stress_99th_pct / S_vm_mean;

    % 6. PLANAR DIC SURFACE EXTRACTION (y = 0 plane)
    y_min_bound = min(meshData.V(:, 2));
    dic_face_indices = find(meshData.V(:, 2) <= y_min_bound + 1e-3);
    
    results_summary.dic_coords = meshData.V(dic_face_indices, :);
    results_summary.dic_disp_x = N_disp_mat(dic_face_indices, 1, end);
    results_summary.dic_disp_y = N_disp_mat(dic_face_indices, 2, end);
    results_summary.dic_disp_z = N_disp_mat(dic_face_indices, 3, end);

    fprintf('    -> FEBio Mean Z-Load Disp: %.6f mm | TopOpt Mean Z-Load Disp: %.6f mm\n', ...
            results_summary.mean_U_z_load_FEBio, results_summary.mean_U_z_load_TopOpt);
    fprintf('    -> Mean Stress: %.4f MPa | Peak Stress: %.2f MPa | Ratio: %.3f\n', ...
            results_summary.mean_stress, results_summary.stress_99th_pct, results_summary.stress_peak_to_mean_ratio);
end