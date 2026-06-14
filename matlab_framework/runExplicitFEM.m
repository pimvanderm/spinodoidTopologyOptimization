function [runFlag, meshData, febioPaths] = runExplicitFEM(tr, settings)
    % RUNEXPLICITFEM Meshes a spinodoid STL and runs a 3PB FEBio analysis
    %
    % INPUTS:
    %   tr       - MATLAB triangulation object of the explicit spinodoid
    %   settings - Struct containing: .numberOfWaves, .sampleSize, .savePath
    %
    % OUTPUTS:
    %   runFlag    - 1 if FEBio converged, 0 otherwise
    %   meshData   - Struct containing V, E, elemVol, bcPrescribeList
    %   febioPaths - Struct containing paths to generated log files
    %%
    fprintf('\n--- EXPLICIT FEM PIPELINE: Meshing & Solving ---\n');
    F_raw = tr.ConnectivityList;
    V_raw = tr.Points;

    % 1. DYNAMIC POINT SPACING (5 elements per wavelength rule)
    pointSpacing = settings.sampleSize / (5 * settings.numberOfWaves);
    fprintf('1. Geogram Remeshing (Wavenumber: %d -> Spacing: %.4f mm)...\n', ...
            settings.numberOfWaves, pointSpacing);
    
    optionStruct.pointSpacing = pointSpacing;
    [F_remesh, V_remesh] = ggremesh(F_raw, V_raw, optionStruct);
    C_marker = zeros(size(F_remesh, 1), 1);

    % 2. VOLUME TETRAHEDRALIZATION
    fprintf('2. TetGen Volume Generation...\n');
    inputStruct.stringOpt = '-pq2AY'; 
    inputStruct.Faces = F_remesh;
    inputStruct.Nodes = V_remesh;
    inputStruct.holePoints = [];
    inputStruct.faceBoundaryMarker = C_marker;
    % --- Flat Top Isocap Inward Step Method ---

    % 1. Find the maximum Z-coordinate of your mesh to identify the top isocap
    zMax = max(V_remesh(:, 3));
    
    % 2. Find all vertices that sit on this top cap (using a small tolerance for floating point)
    tol = 1e-5; 
    isocapIndices = find(abs(V_remesh(:, 3) - zMax) < tol);
    
    if isempty(isocapIndices)
        error('Could not find vertices on the specified isocap plane.');
    end
    
    % 3. Pick a vertex near the middle of this isocap population to avoid edge effects
    medianIdx = round(length(isocapIndices) / 2);
    targetVertexIdx = isocapIndices(medianIdx);
    isocapPoint = V_remesh(targetVertexIdx, :);
    
    % 4. Step safely inside the volume
    % Since it is the TOP cap, the material is strictly BELOW it (-Z direction)
    stepDistance = 0.1; % 0.1 mm step
    innerPoint = isocapPoint;
    innerPoint(3) = innerPoint(3) - stepDistance; % Subtract from Z
    
    % Assign directly to your input structure
    inputStruct.regionPoints = innerPoint;
    inputStruct.regionA = 10 * tetVolMeanEst(F_remesh, V_remesh); 
    inputStruct.minRegionMarker = 2;

    [meshOutput] = runTetGen(inputStruct);
    V = meshOutput.nodes;
    E = meshOutput.elements;
    elemVol = tetVol(E, V);
    
    % 3. DYNAMIC Z-ANCHORED BOUNDARY CONDITIONS
    z_min_active = min(V(:, 3)); 
    z_max_active = max(V(:, 3)); 
    dimensions = max(V);
    xmin = min(V(:, 1));

    constraint1 = [(dimensions(1)+xmin)/2 - 2*dimensions(3), z_min_active];
    constraint2 = [(dimensions(1)+xmin)/2 + 2*dimensions(3), z_min_active];
    load1       = [(dimensions(1)+xmin)/2,                   z_max_active];

    tolDir = max(1.2 * pointSpacing, 0.5); 

    logic1 = vecnorm(V(:, [1 3]) - constraint1, 2, 2) < tolDir; 
    logic2 = vecnorm(V(:, [1 3]) - constraint2, 2, 2) < tolDir; 
    logic3 = vecnorm(V(:, [1 3]) - load1,       2, 2) < tolDir; 

    bcSupportList   = find(logic1 | logic2);
    bcPrescribeList = find(logic3);

    y_center = (max(V(:, 2)) + min(V(:, 2))) / 2;
    x_center = (dimensions(1) + xmin) / 2;
    [~, bcCenterNode] = min(vecnorm(V - [x_center, y_center, z_max_active], 2, 2));

    bcSupportList   = setdiff(bcSupportList, bcCenterNode);
    bcPrescribeList = setdiff(bcPrescribeList, bcCenterNode);

    if isempty(bcSupportList) || isempty(bcPrescribeList)
        error('Adaptive BC mapping failed: Zero nodes found on boundary lines.');
    end

    % 4. BUILD FEBIO XML & EXECUTE

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
    % File naming definitions
    savePath = fullfile(fileparts(mfilename('fullpath')), 'data', 'temp');
    if ~exist(savePath, 'dir'), mkdir(savePath); end
    
    febioFebFileNamePart = '3PB_Convergence';
    febioFebFileName = fullfile(savePath, [febioFebFileNamePart, '.feb']);
    febioLogFileName = [febioFebFileNamePart, '.txt'];
    febioLogFileName_disp = [febioFebFileNamePart, '_disp_out.txt'];
    febioLogFileName_stress_prin = [febioFebFileNamePart, '_stress_prin_out.txt'];
    febioLogFileName_stress_full = [febioFebFileNamePart, '_stress_full_out.txt'];

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

    totalLoad = 50;
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
    
    febioPaths.xml    = febioFebFileName;
    febioPaths.log    = fullfile(savePath, febioLogFileName);
    febioPaths.disp   = fullfile(savePath, febioLogFileName_disp);
    febioPaths.stress = fullfile(savePath, febioLogFileName_stress_prin);

    % Package mesh data for post-processing
    meshData.V = V;
    meshData.E = E;
    meshData.elemVol = elemVol;
    meshData.bcPrescribeList = bcPrescribeList;
end