% =========================================================================
% MASTER MULTISCALE SPINODOID PIPELINE
% Sequentially executes TopOpt, Rendering, FEA Validation, and Data Harvest
% =========================================================================


%% 1. INITIALIZE GLOBAL PARAMETERS
% TopOpt parameters
params.nelx = 44;
params.nely = 10;
params.nelz = 10;
params.problemType = 'three_point_bending';
params.rho_min = 0.3;
params.rho_max = 0.6;
params.theta_min = 15;
params.lambda1 = 50;
params.lambda2 = 50;
params.rho_nn_min = 0;
params.rho_nn_max = 1;
params.theta_nn_min = 0;
params.theta_nn_max = 90;
params.C_scale = 2800;
params.dlnetF_path = 'dlnetF.mat';
params.eps = 1e-3;
params.maxiter = 1000;
params.tol = 5e-4;
params.stagnation_tol = 5e-3;
params.volfrac = 0.3;
params.move = 0.15;
params.beta = 5;
params.eta = 0.5;
params.r_filter = 1.5;
params.alpha_max = 180;
params.nele = params.nelx * params.nely * params.nelz;
params.rho_init = params.volfrac;
params.theta_init = [0 15 15];
params.alpha_init = 0;
%%
% Rendering and Validation Settings
settings.nelx = params.nelx;
if strcmp(params.problemType, 'three_point_bending_half')
    settings.nelx = 2 * params.nelx;
end
settings.nely = params.nely;
settings.nelz = params.nelz;
settings.domainSize = [settings.nelx, settings.nely, settings.nelz];
settings.volfrac = params.volfrac;
settings.problemType = params.problemType;
settings.numberOfWaves = 10; % Baseline test wavenumber
settings.sampleSize = 20; % Assuming the Z-height maps to 20mm
settings.res = 8 * settings.numberOfWaves / params.nelz; 
settings.k_wave = settings.numberOfWaves * 1 / params.nelz * 2 * pi;
settings.writeSTL = 0;
%%
% Setup save directory for FEBio
defaultFolder = fileparts(mfilename('fullpath'));
settings.savePath = fullfile(defaultFolder, 'data', 'temp');
if ~exist(settings.savePath, 'dir')
    mkdir(settings.savePath);
end

%% =========================================================================
% STAGE 1: MACRO-SCALE TOPOLOGY OPTIMIZATION
% =========================================================================
fprintf('\n>>> STAGE 1: Running Topology Optimization...\n');
[chi_raw_history, chi_history, c_history, U_final] = optimizeTopology(params);

% Extract the final valid iteration
finalIter = find(all(chi_raw_history(:,:,1) ~= 0, 2), 1, 'last');
chi_final_raw = squeeze(chi_raw_history(finalIter,:,:));

% Snap to bounds
chi_final = snapToValidRanges(chi_final_raw, params);
if strcmp(params.problemType, 'three_point_bending_half')
    [chi_final, params] = mirrorHalfStructure(chi_final,params);
end


% --- NEW INSERTION: EXECUTE PHYSICAL RESCALING ---
scaledTopOpt = scaleTopOptOutputs(U_final, c_history, finalIter, params, settings);

% Record the TRUE PHYSICALLY SCALED compliance for downstream energy verification
settings.c_TopOpt = scaledTopOpt.c_final_physical; 
settings.U_topopt_scaled = scaledTopOpt.U_physical;
settings.c_TopOpt = scaledTopOpt.c_final_physical;

% Calculate predicted physical volume of the continuum model
% Total domain volume * target volume fraction
domain_vol_mm3 = (params.nelx * params.nely * params.nelz) * (settings.sampleSize / params.nelz)^3;
settings.V_predicted = domain_vol_mm3 * mean(chi_final(:,1));


%% =========================================================================
% STAGE 2: EXPLICIT SPINODOID RENDERING
% =========================================================================
fprintf('\n>>> STAGE 2: Rendering Explicit Spinodoid GRF (N = %d)...\n', settings.numberOfWaves);
tr = renderSpinodoidStructure(chi_final, settings);

%% =========================================================================
% STAGE 3: HIGH-FIDELITY FEM VALIDATION
% =========================================================================
fprintf('\n>>> STAGE 3: Generating TetMesh and Running FEBio...\n');
[runFlag, meshData, febioPaths] = runExplicitFEM(tr, settings);

%% =========================================================================
% STAGE 4: POST-PROCESSING & METRICS HARVESTING
% =========================================================================
if runFlag == 1
    fprintf('\n>>> STAGE 4: Harvesting Data...\n');
    results_summary = postProcessFEM(tr, meshData, febioPaths, settings);
    
    % Display Final Summary
    fprintf('\n====================================================\n');
    fprintf('PIPELINE COMPLETE - FINAL METRICS\n');
    fprintf('====================================================\n');
    fprintf('Predicted TopOpt Volume:      %.2f mm^3\n', settings.V_predicted);
    fprintf('Actual Rendered Volume:       %.2f mm^3\n', results_summary.V_rendered);
    fprintf('Volume Ratio:                 %.3f\n', results_summary.V_ratio);
    fprintf('----------------------------------------------------\n');
    fprintf('TopOpt Mean Z-Load Disp:      %.6f mm\n', results_summary.mean_U_z_load_TopOpt);
    fprintf('FEBio Mean Z-Load Disp:       %.6f mm\n', abs(results_summary.mean_U_z_load_FEBio));
    fprintf('----------------------------------------------------\n');
    fprintf('Volume-Weighted Mean Stress:  %.4f MPa\n', results_summary.mean_stress);
    fprintf('99th Pct Stress (Peak):       %.2f MPa\n', results_summary.stress_99th_pct);
    fprintf('Stress Concentration Ratio:   %.3f\n', results_summary.stress_peak_to_mean_ratio);
    fprintf('====================================================\n');
else
    fprintf('\n[!] Pipeline halted: FEBio failed to converge.\n');
end

% Note: If snapToValidRanges is currently bundled inside optimizeTopology.m as a local function, 
% you will need to extract it into its own file named `snapToValidRanges.m` so this master script can call it!