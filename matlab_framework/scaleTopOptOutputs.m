function [scaledResults] = scaleTopOptOutputs(U_final, c_history, finalIter, params, settings)
    fprintf('--- DIMENSIONAL SCALING: Mapping TopOpt to Physical Domain ---\n');

    % 1. DEFINE PHYSICAL CONTEXT BENCHMARKS
    height_physical = 20;         % Target physical height of the structure (mm)
    
    % 2. COMPUTE LENGTH SCALE
    % Length of a single voxel element in real space (mm per element)
    L_0 = height_physical / params.nelz; 
    
    % 3. COMPUTE PURE GEOMETRIC SCALING MULTIPLIERS
    % Because C_scale (2800 MPa) and totalLoad (25 N) are already 
    % embedded in the TopOpt K and F matrices, the only missing scaling is geometry.
    S_u = 1 / L_0; 
    
    % Compliance (Work) scales with 1/L_0, BUT we also multiply by 2 
    % to map the 25N half-model energy to the 50N full-model energy!
    S_c = 2 / L_0;

    % 4. TRANSLATE TO REAL PHYSICAL QUANTITIES
    scaledResults.U_physical = U_final * S_u;
    scaledResults.c_history_physical = c_history * S_c;
    scaledResults.c_final_physical = c_history(finalIter) * S_c;
    scaledResults.L_0 = L_0;
    
    % Calculate the maximum deflection component present on the continuum model in mm
    scaledResults.max_U_topopt_mm = max(sqrt(scaledResults.U_physical(1:3:end).^2 + ...
                                              scaledResults.U_physical(2:3:end).^2 + ...
                                              scaledResults.U_physical(3:3:end).^2));

    fprintf('    -> Voxel Length Scale (L_0): %.4f mm/voxel\n', L_0);
    fprintf('    -> Continuum Max Displacement: %.6f mm\n', scaledResults.max_U_topopt_mm);
end