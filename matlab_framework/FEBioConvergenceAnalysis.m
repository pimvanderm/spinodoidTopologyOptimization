% --- FEBio Mesh Convergence Analysis Script ---

clear all
% Load the raw explicit STL structure
tr = stlread('TopOpt_three_point_bending_half_88x20x20_v03_k10_2026-06-05_103315.stl');
F_raw = tr.ConnectivityList;
V_raw = tr.Points;

% Define 5 levels of structural mesh density (Coarse to Fine)
mesh_levels = [15, 20, 25, 30, 35, 40, 50, 60]; 
num_steps = length(mesh_levels);

% Preallocate convergence logs
convergence_results = struct('num_elements', [], 'max_displacement', [], ...
                             'mean_stress', [], 'vol_fraction', []);
%%

%% =========================================================================
% INITIALIZE GLOBAL SIMULATION CONTROLS (Missed variables)
% =========================================================================
E_youngs = 2800; % Scaled Young's Modulus (GPa)
nu = 0.33;             % Poisson's ratio

numTimeSteps = 1;      % Number of time steps
max_refs = 15;         % Max reforms
max_ups = 0;           % Set to zero for full-Newton iterations
opt_iter = 6;          % Optimum number of iterations
max_retries = 5;       % Maximum number of retries
dtmin = (1/numTimeSteps)/100; % Minimum time step size
dtmax = 1/numTimeSteps;       % Maximum time step size

runMode = 'external';  % FEBio run mode
sampleSize = 20; 

% File naming definitions
savePath = fullfile(fileparts(mfilename('fullpath')), 'data', 'temp');
if ~exist(savePath, 'dir'), mkdir(savePath); end

febioFebFileNamePart = '3PB_Convergence';
febioFebFileName = fullfile(savePath, [febioFebFileNamePart, '.feb']);
febioLogFileName = [febioFebFileNamePart, '.txt'];
febioLogFileName_disp = [febioFebFileNamePart, '_disp_out.txt'];
febioLogFileName_stress_prin = [febioFebFileNamePart, '_stress_prin_out.txt'];
febioLogFileName_stress_full = [febioFebFileNamePart, '_stress_full_out.txt'];

%% =========================================================================
% THE PARAMETRIC CONVERGENCE LOOP
% =========================================================================
for m = 1:num_steps
    fprintf('\n==================================================\n');
    fprintf('RUNNING MESH LEVEL %d/%d (Resolution Factor: %d)\n', m, num_steps, mesh_levels(m));
    fprintf('==================================================\n');
    
    % Clear old loop structures to avoid memory bleed
    clear inputStruct meshStruct optionStruct febio_spec;

    % 1. Set Surface Point Spacing dynamically (Fixed variable references)
    res_factor = mesh_levels(m);
    pointSpacing = sampleSize / res_factor;
    fprintf('   [Stage 1/5] Executing Geogram surface remeshing (pointSpacing = %.5f mm)...\n', pointSpacing);
    optionStruct.pointSpacing = pointSpacing;
    [F_remesh, V_remesh] = ggremesh(F_raw, V_raw, optionStruct); % Fixed mapping link
    C_marker = zeros(size(F_remesh, 1), 1);
    fprintf('   [Stage 2/5] Initializing TetGen volume element generation...\n');
    % 2. Set Volume Constraint dynamically
    inputStruct.stringOpt = '-pq2AY'; % Standardized quality ratio constraint
    inputStruct.Faces = F_remesh;
    inputStruct.Nodes = V_remesh;
    inputStruct.holePoints = [];
    inputStruct.faceBoundaryMarker = C_marker;
    inputStruct.regionPoints = getInnerPoint(F_remesh, V_remesh);
    inputStruct.regionA = 10 * tetVolMeanEst(F_remesh, V_remesh); % Standardized multiplier
    inputStruct.minRegionMarker = 2;

    % Execute TetGen
    [meshOutput] = runTetGen(inputStruct);

    % Access generated mesh coordinates
    Fb = meshOutput.facesBoundary;
    V = meshOutput.nodes;
    E = meshOutput.elements;
    num_elements = size(E, 1);
    elemVol = tetVol(E, V);
    fprintf('   [Stage 3/5] Mapping spatial coordinates to track support and loading node boundaries...\n');
    % 3. DYNAMIC BOUNDARY CONDITION SELECTION (CRITICAL MISSING STEP)
    % Coordinates vary with mesh resolution, recalculate for current grid:
    tolDir = 1;
    C_vertex = zeros(size(V, 1), 1);
    dimensions = max(V_raw);
    xmin = min(V_raw(:, 1));

    constraint1 = [(dimensions(1)+xmin)/2-2*dimensions(3), 0];
    constraint2 = [(dimensions(1)+xmin)/2+2*dimensions(3), 0];
    load1 = [(dimensions(1)+xmin)/2, dimensions(3)];

    logic1 = vecnorm(V(:, [1 3]) - constraint1, 2, 2) < tolDir; % Left Support
    logic2 = vecnorm(V(:, [1 3]) - constraint2, 2, 2) < tolDir; % Right Support
    logic3 = vecnorm(V(:, [1 3]) - load1, 2, 2) < tolDir;       % Loading Line

    C_vertex(logic1) = 1;
    C_vertex(logic2) = 1;
    C_vertex(logic3) = 2;
    
    % Define boundary condition lists based on C_vertex
    bcSupportList   = find(C_vertex == 1); % Both roller support lines combined
    bcPrescribeList = find(C_vertex == 2); % The entire top loading line

    % Find the center node along the loading axis to eliminate axial movement
    y_center = (max(V(:, 2)) + min(V(:, 2))) / 2;
    [~, center_idx_within_load] = min(abs(V(bcPrescribeList, 2) - y_center));
    bcCenterNode = bcPrescribeList(center_idx_within_load);

    fprintf('   [Stage 4/5] Compiling febio_spec XML tree and dispatching to external FEBio solver engine...\n');
    %% 4. Setup FEBio Specification XML
    febio_spec = febioStructTemplate; 
febio_spec.ATTR.version = '4.0';
febio_spec.Module.ATTR.type = 'solid';
febio_spec.Globals.Constants.T = 0;
febio_spec.Globals.Constants.R = 0;
%Control section
febio_spec.Control.analysis='STATIC';
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=1/numTimeSteps;
febio_spec.Control.solver.max_refs=max_refs;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax; 
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;

%Material section
materialName1='Material1';
febio_spec.Material.material{1}.ATTR.name=materialName1;
febio_spec.Material.material{1}.ATTR.type='isotropic elastic';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.E=E_youngs;
febio_spec.Material.material{1}.v=nu;

% Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='Object1'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='tet4'; %Element type
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E; %The element matrix

% -> NodeSets (Configuring the 3 distinct groups)
nodeSetName1 = 'bcSupportRollers';
nodeSetName2 = 'bcLoadingLine';
nodeSetName3 = 'bcCenterStabilizerNode';

febio_spec.Mesh.NodeSet{1}.ATTR.name = nodeSetName1;
febio_spec.Mesh.NodeSet{1}.VAL       = mrow(bcSupportList);

febio_spec.Mesh.NodeSet{2}.ATTR.name = nodeSetName2;
febio_spec.Mesh.NodeSet{2}.VAL       = mrow(bcPrescribeList);

febio_spec.Mesh.NodeSet{3}.ATTR.name = nodeSetName3;
febio_spec.Mesh.NodeSet{3}.VAL       = mrow(bcCenterNode);

%MeshDomains section
febio_spec.MeshDomains.SolidDomain.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain.ATTR.mat=materialName1;

% ---- Boundary conditions -------------------------------------------------

% 1. Roller Supports: Only constrain vertical movement (z-direction)
febio_spec.Boundary.bc{1}.ATTR.name     = 'RollerSupports_Z';
febio_spec.Boundary.bc{1}.ATTR.type     = 'zero displacement';
febio_spec.Boundary.bc{1}.ATTR.node_set = nodeSetName1;
febio_spec.Boundary.bc{1}.x_dof         = 0;
febio_spec.Boundary.bc{1}.y_dof         = 0;
febio_spec.Boundary.bc{1}.z_dof         = 1;

% 2. Top Loading Line: Lock the x-direction to prevent the model from sliding left/right
febio_spec.Boundary.bc{2}.ATTR.name     = 'LoadingLine_X';
febio_spec.Boundary.bc{2}.ATTR.type     = 'zero displacement';
febio_spec.Boundary.bc{2}.ATTR.node_set = nodeSetName2;
febio_spec.Boundary.bc{2}.x_dof         = 1;
febio_spec.Boundary.bc{2}.y_dof         = 0;
febio_spec.Boundary.bc{2}.z_dof         = 0;

% 3. Middle Node: Lock the y-direction to stop rigid body axial telescoping
febio_spec.Boundary.bc{3}.ATTR.name     = 'CenterNode_Y';
febio_spec.Boundary.bc{3}.ATTR.type     = 'zero displacement';
febio_spec.Boundary.bc{3}.ATTR.node_set = nodeSetName3;
febio_spec.Boundary.bc{3}.x_dof         = 0;
febio_spec.Boundary.bc{3}.y_dof         = 1;
febio_spec.Boundary.bc{3}.z_dof         = 0;

% ---- External Loads ------------------------------------------------------

totalLoad = 1;
nodalForce = totalLoad/length(bcPrescribeList);
febio_spec.Loads.nodal_load{1}.ATTR.type     = 'nodal_force';
febio_spec.Loads.nodal_load{1}.ATTR.node_set = nodeSetName2;
febio_spec.Loads.nodal_load{1}.value.ATTR.lc = 1;
febio_spec.Loads.nodal_load{1}.value.VAL = [0,0,-nodalForce];

% ---- Load controller -----------------------------------------------------
febio_spec.LoadData.load_controller{1}.ATTR.name = 'LC1';
febio_spec.LoadData.load_controller{1}.ATTR.id   = 1;
febio_spec.LoadData.load_controller{1}.ATTR.type = 'loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate = 'LINEAR';
febio_spec.LoadData.load_controller{1}.extend      = 'CONSTANT';
febio_spec.LoadData.load_controller{1}.points.pt.VAL = [0 0; 1 1];

%Output section
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress_prin;
febio_spec.Output.logfile.element_data{1}.ATTR.data='s1;s2;s3';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';

febio_spec.Output.logfile.element_data{2}.ATTR.file=febioLogFileName_stress_full;
febio_spec.Output.logfile.element_data{2}.ATTR.data='sx;sy;sz;sxy;syz;sxz';
febio_spec.Output.logfile.element_data{2}.ATTR.delim=',';

febio_spec.Output.plotfile.compression = 0;

febioStruct2xml(febio_spec, febioFebFileName);
        
febioAnalysis.run_filename=febioFebFileName; %The input file name
febioAnalysis.run_logname=febioLogFileName; %The name for the log file

% CRITICAL: Point directly to the .exe here to bypass the GUI
febio_exe_path = 'C:\Program Files\FEBioStudio\bin\febio4.exe'; % Ensure this is correct
febioAnalysis.FEBioPath = febio_exe_path;
febioAnalysis.disp_on=1; %Display information on the command window
febioAnalysis.runMode=runMode;
febioAnalysis.maxLogCheckTime = 300; % Maximum wait time for log file creation (seconds)
febioAnalysis.mesh_scale = 1.0;

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!

    %% 5. Harvest Convergence Results
    if runFlag == 1
        fprintf('   [Stage 5/5] FEBio run successful. Harvesting log metrics and calculating fields...\n');
        % Extract Displacement profiles
        dataStruct = importFEBio_logfile(fullfile(savePath, febioLogFileName_disp), 0, 1);
        N_disp_mat = dataStruct.data;
        
        % Extract Principal Stress profiles
        dataStruct_stress = importFEBio_logfile(fullfile(savePath, febioLogFileName_stress_prin), 0, 1);
        E_stress_mat = dataStruct_stress.data;
        E_stress_mat(isnan(E_stress_mat)) = 0;
        
        % Calculate local Von Mises stresses
        S_vm_ND = sqrt(0.5 * ((E_stress_mat(:,1,:) - E_stress_mat(:,2,:)).^2 + ...
                             (E_stress_mat(:,2,:) - E_stress_mat(:,3,:)).^2 + ...
                             (E_stress_mat(:,1,:) - E_stress_mat(:,3,:)).^2));
        
        totalVol = sum(elemVol);
        S_vm_mean = squeeze(sum(S_vm_ND .* repmat(elemVol, [1 1 size(S_vm_ND, 3)]), 1) ./ totalVol);
        
        % Commit structural data metrics to the master summary log
        convergence_results(m).num_elements = num_elements;
        convergence_results(m).vol_fraction = totalVol / (max(V_raw(:,1))*max(V_raw(:,2))*max(V_raw(:,3)));
        convergence_results(m).max_displacement = max(sqrt(sum(N_disp_mat(:, :, end).^2, 2)));
        convergence_results(m).mean_stress = S_vm_mean(end);
    else
        error('Mesh level %d failed to converge in FEBio.', mesh_levels(m));
    end
end

%% 6. Generate Verification Plots
figure('Name', 'Mesh Convergence Study', 'Color', 'w');
subplot(2,1,1);
plot(mesh_levels, [convergence_results.max_displacement], 'b-o', 'LineWidth', 2);
grid on; xlabel('Number of Elements'); ylabel('Max Displacement [mm]');
title('Displacement Convergence (Plateau verification)');

subplot(2,1,2);
plot(mesh_levels, [convergence_results.mean_stress], 'r-s', 'LineWidth', 2);
grid on; xlabel('Number of Elements'); ylabel('Mean Von Mises Stress [MPa]');
title('Stress Convergence Profile');