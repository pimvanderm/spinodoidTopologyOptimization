%% DEMO_febio_0010_trabeculae_compression
% Below is a demonstration for:
%
% * Building geometry for trabecular structure with tetrahedral elements
% * Defining the boundary conditions
% * Coding the febio structure
% * Running the model
% * Importing and visualizing the displacement and stress results

%% Keywords:
% * febio_spec version 4.0
% * febio, FEBio
% * compression, tension, compressive, tensile
% * displacement control, displacement boundary condition
% * trabecular
% * tetgen, meshing
% * tetrahedral elements, tet4
% * static, solid
% * hyperelastic, Ogden
% * displacement logfile
% * Stress logfile

clear all; close all;
profile on % Start the clock
%% INITIALIZE PROGRESS BAR
totalSimulations = 3 * 11 * 6; % 198 total load cases
simCounter = 0;
hWait = waitbar(0, sprintf('Initializing FEBio Pipeline... (0/%d)', totalSimulations), ...
                'Name', 'Homogenization Sweep Progress');
%% 0. INITIALIZE FAIL-SAFE CSV LOGGER
defaultFolder = fileparts(mfilename('fullpath'));
settings.savePath = fullfile(defaultFolder, 'data', 'temp');
if ~exist(settings.savePath, 'dir')
    mkdir(settings.savePath);
end

csvFileName = fullfile(settings.savePath, 'Homogenization_Stiffness_Ledger.csv');

% If the file doesn't exist yet, create it and write the header row
if ~isfile(csvFileName)
    fid = fopen(csvFileName, 'w');
    fprintf(fid, 'Seed,LoopIndex,Weight1,Rho1,T1_X,T1_Y,T1_Z,Rho2,T2_X,T2_Y,T2_Z,C11,C12,C13,C22,C23,C33,C44,C55,C66\n');
    fclose(fid);
end

%% 1. NESTED SEED LOOP
seedList = [1, 2, 3];

for s_idx = 1:length(seedList)
    currentSeed = seedList(s_idx);
    rng(currentSeed); % Enforce exact stochasticity for this run
    
    fprintf('\n=========================================================\n');
    fprintf('>>> STARTING MACRO ITERATION FOR SEED %d/3 (RNG: %d) <<<\n', s_idx, currentSeed);
    fprintf('=========================================================\n');
    
    tic;
    results = zeros(6,6,11); % Reinitialize the stiffness container for this seed

for loop = 0:10
    fprintf('\n Starting loop %d \n',(loop+1))
    clear inputs
    clear inputStruct
    clear meshStruct
    
    dimX = 1;
    dimY = 1;
    dimZ = 1;
    domainSize = [dimX,dimY,dimZ];
    k = 10;
    res = 8*k;
    resolution = res * domainSize;
    nEle = 2;
    inputs(nEle) = struct();
    num_waves = 500;
    padding = 1;
    neighborhoodRadius = padding / 2; % From your original function
    kappa = 3;
    k_matrix = k*2*pi*ones(nEle);
    %% Example design variable matrices (these would come from the TopOpt scheme)
    
    theta1 = [15 0 0];
    theta2 = [0 15 15];
    theta_matrix = [theta1; theta2];
    rho1 = 0.6;
    rho2 = 0.6;
    rho_matrix = [rho1 rho2];
    orient_matrix = [0 0 0;
                     0 0 0];
                     
    weight1 = 0.1*loop;
    weight2 = 1 - weight1;
    weights = [weight1 weight2];
    
    %%
    
    % Generate all element centers in one line
    [Z, Y, X] = ndgrid(1:dimZ, 1:dimY, 1:dimX);
    
    % Element centers (subtract 0.5 for center)
    cx = X(:) - 0.5; cy = Y(:) - 0.5; cz = Z(:) - 0.5;
    centers = [cx, cy, cz];
    centers = repmat(centers,2,1);
    
    for idx = 1:nEle
        inputs(idx).relativeDensity = rho_matrix(idx); % Use valid density for GRF generation
        inputs(idx).thetas = theta_matrix(idx,:);
        inputs(idx).waveNumber = k_matrix(idx);
        inputs(idx).numWaves = num_waves;
        inputs(idx).isocap = true;
        padding = 1;
        inputs(idx).domainSize = [1, 1, 1] * padding; % Do not render for the full volume; instead, render only the element, with a small padding around it for the interpolation
        inputs(idx).resolution = round(res * padding + 1) * [1, 1, 1];
        inputs(idx).ignoreChecks = false;
    end
    
    
    
    % --- Pre-calculate Weights and Global Grid ---
    
    % Generate the TRUE global grid for the whole beam
    xv = linspace(0.5/res, dimX - 0.5/res, resolution(1));
    yv = linspace(0.5/res, dimY - 0.5/res, resolution(2));
    zv = linspace(0.5/res, dimZ - 0.5/res, resolution(3));
    [X, Y, Z] = ndgrid(xv, yv, zv);
    
    
    % --- Setup Global Grid Containers ---
    % We initialize the global container as empty zeros
    graded_GRF = zeros(dimX*res, dimY*res, dimZ*res); 
    
    % --- The Optimized Accumulation Loop ---
    fprintf('Starting optimized generation...\n');
    start_time = toc;
    
    for idx = 1:nEle
        
        [X_q, Y_q, Z_q, weights_i,limits] = calcWeights(centers(idx,:),res,resolution,neighborhoodRadius,kappa);
        [GRF_src, X_src, Y_src, Z_src, levelset_i] = genGRF(inputs(idx));
        offset = (padding - 1) / 2; % Shift source to be centered at 0.5
        X_src = X_src - offset; Y_src = Y_src - offset; Z_src = Z_src - offset;
        F_interp = griddedInterpolant(X_src, Y_src, Z_src, GRF_src, 'linear', 'nearest');
        GRF_i = F_interp(X_q, Y_q, Z_q) - levelset_i;
        
        X_min = limits(1,1);X_max = limits(1,2);
        Y_min = limits(2,1);Y_max = limits(2,2);
        Z_min = limits(3,1);Z_max = limits(3,2);
    
        % Direct Indexing into Global Matrix
        graded_GRF(X_min:X_max,Y_min:Y_max,Z_min:Z_max) = ...
            graded_GRF(X_min:X_max,Y_min:Y_max,Z_min:Z_max) + ...
            weights(idx) * GRF_i;
    
        
    end
    
    fprintf('Generated Gaussian Random Field...\n');
    %% INITIAL PARAMETERS
    
    L = 1.0;                % Side length of RVE (mm)
    
    % Solid phase properties (will be used in later steps)
    E_solid = 1;            
    nu_solid = 0.3;
    
    
  
    %% 1. Surface Extraction (Isosurface)
    % Define the level-set (0 for solid/void boundary)
    isoLevel = 0;
    
    % Extract the surface triangulated mesh
    [F, V] = isosurface(X,Y,Z,graded_GRF, isoLevel);
    C=zeros(size(F,1),1);
    %Compute isocaps
    [fc,vc] = isocaps(X,Y,Z,graded_GRF,isoLevel,'enclose','below');
    
    if ~isempty(fc)
        nc=patchNormal(fc,vc);
        cc=zeros(size(fc,1),1);
        cc(nc(:,1)<-0.5)=1;
        cc(nc(:,1)>0.5)=2;
        cc(nc(:,2)<-0.5)=3;
        cc(nc(:,2)>0.5)=4;
        cc(nc(:,3)<-0.5)=5;
        cc(nc(:,3)>0.5)=6;
    
        %Join sets
        [F,V,C]=joinElementSets({F,fc},{V,vc},{C,cc});
    end
    
    %% Merge nodes and clean-up mesh 
    
    %Merge nodes
    [F,V]=mergeVertices(F,V); 
    
    %Check for unique faces
    [~,indUni,~]=unique(sort(F,2),'rows');
    F=F(indUni,:); %Keep unique faces
    C=C(indUni);
    
    %Remove collapsed faces
    [F,logicKeep]=patchRemoveCollapsed(F); 
    C=C(logicKeep);
    
    %Remove unused points
    [F,V]=patchCleanUnused(F,V);
    
  
    % 
    %% Plot settings
    fontSize=20;
    faceAlpha1=0.8;
    markerSize=40;
    lineWidth1=3;
    lineWidth2=4;
    markerSize1=25;
    % 

    % cFigure; hold on;
    % title(num2str(loop),'FontSize',fontSize);
    % gpatch(F,V,C,'k',1);
    % 
    % % plotV(V(indKeep,:),'k.','MarkerSize',markerSize1);
    % axisGeom(gca,fontSize);
    % colormap gjet; icolorbar;
    % camlight headlight;
    % drawnow;


    %% Control parameters
   
    % Path names
    defaultFolder = fileparts(fileparts(mfilename('fullpath')));
    savePath=fullfile(defaultFolder,'data','temp');
    %%
    % Defining file names
    febioFebFileNamePart='tempModel';
    febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
    febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
    febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
    febioLogFileName_stress_prin=[febioFebFileNamePart,'_stress_prin_out.txt']; %Log file name for exporting stress
    febioLogFileName_stress_full=[febioFebFileNamePart,'_stress_full_out.txt']; %Log file name for exporting stress
    % 
    % porousGeometryCase=3;
    % 
    sampleSize=1; %Height of the sample
    pointSpacing=sampleSize/(5*k);
    tolDir=pointSpacing; %Tolerance for detecting sides after remeshing
    % 
    overSampleRatio=2;
    numStepsLevelset=ceil(overSampleRatio.*(sampleSize./pointSpacing)); %Number of voxel steps across period for image data (roughly number of points on mesh period)
    % 
    %Define applied displacement
    appliedStrain=0.001; %Linear strain (Only used to compute applied stretch)
    
    loadingOption='compression'; % or 'tension' or 'shear'
    switch loadingOption
        case 'compression'
            stretchLoad=1-appliedStrain; %The applied stretch for uniaxial loading
            displacementMagnitude=(stretchLoad*sampleSize)-sampleSize; %The displacement magnitude
        case 'tension'
            stretchLoad=1+appliedStrain; %The applied stretch for uniaxial loading
            displacementMagnitude=(stretchLoad*sampleSize)-sampleSize; %The displacement magnitude
        case 'shear'
            stretchLoad=1+appliedStrain; %The applied stretch for uniaxial loading
            displacementMagnitude=(stretchLoad*sampleSize)-sampleSize; %The displacement magnitude
    end
    % 
    %Material parameter set
    E_youngs=1; %Youngs modulus
    nu=0.3; %Poisson's ratio
    % mu=E_youngs/3;
    % 
    %FEA control settings
    numTimeSteps=1; %Number of time steps desired
    max_refs=15; %Max reforms
    max_ups=0; %Set to zero to use full-Newton iterations
    opt_iter=6; %Optimum number of iterations
    max_retries=5; %Maximum number of retires
    dtmin=(1/numTimeSteps)/100; %Minimum time step size
    dtmax=1/numTimeSteps; %Maximum time step size
    
    runMode='external';% 'internal' or 'external'
    
    
    %% Using grouping to keep only largest group
    
    groupOptStruct.outputType='label';
    [G,~,groupSize]=tesgroup(F,groupOptStruct); %Group connected faces
    [~,indKeep]=max(groupSize); %Index of largest group
    
    %Keep only largest group
    F=F(G==indKeep,:); %Trim faces
    C=C(G==indKeep,:); %Trim color data
    [F,V]=patchCleanUnused(F,V); %Remove unused nodes
    
    
    
    %% Remesh using geomgram
    
    optionStruct.pointSpacing=pointSpacing;
    % optionStruct.max_dist=0;
    [F,V]=ggremesh(F,V,optionStruct);
    C=zeros(size(F,1),1);
    
    
      
    %% Tetrahedral meshing using tetgen (see also |runTetGen|)
    
    % Create tetgen input structure
    inputStruct.stringOpt='-pq2AaY';
    inputStruct.Faces=F;
    inputStruct.Nodes=V;
    inputStruct.holePoints=[];
    inputStruct.faceBoundaryMarker=C; %Face boundary markers
    inputStruct.regionPoints=getInnerPoint(F,V); %region points
    inputStruct.regionA=5*tetVolMeanEst(F,V);
    inputStruct.minRegionMarker=2; %Minimum region marker
    
    % Mesh model using tetrahedral elements using tetGen
    [meshOutput]=runTetGen(inputStruct); %Run tetGen
    
    % Access model element and patch data
    Fb=meshOutput.facesBoundary;
    Cb=meshOutput.boundaryMarker;
    V=meshOutput.nodes;
    CE=meshOutput.elementMaterialID;
    E=meshOutput.elements;
    
    %% Visualizing mesh using |meshView|, see also |anim8|
    
    % meshView(meshOutput);
    
    %% Defining node labels
    
    %% Defining node sets for KUBC
    % Instead of overwriting labels, we create distinct sets for all 6 faces
    nodeSets = struct();
    nodeSets.X_min = find(V(:,1) < tolDir);
    nodeSets.X_max = find(V(:,1) > (sampleSize - tolDir));
    nodeSets.Y_min = find(V(:,2) < tolDir);
    nodeSets.Y_max = find(V(:,2) > (sampleSize - tolDir));
    nodeSets.Z_min = find(V(:,3) < tolDir);
    nodeSets.Z_max = find(V(:,3) > (sampleSize - tolDir));
    
    % Combine all boundary nodes into one master array. 
    % This is crucial for constraining lateral Poisson expansion everywhere.
    allBoundaryNodes = unique([nodeSets.X_min; nodeSets.X_max; ...
                               nodeSets.Y_min; nodeSets.Y_max; ...
                               nodeSets.Z_min; nodeSets.Z_max]);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AUTOMATED HOMOGENIZATION LOOP (6 LOAD CASES)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 1. Initialization
    targetStrain = 0.001; 
    C_homogenized = zeros(6,6); % Initialize the 6x6 stiffness matrix
    elemVol = tetVol(E, V);     % Calculate volumes once
    totalVol = sum(elemVol);
    sampleSize = max(V(:,1)) - min(V(:,1)); % Assuming cubic RVE
    
    % Voigt Labels for documentation
    labels = {'XX', 'YY', 'ZZ', 'XY', 'YZ', 'XZ'};
    febio_exe_path = 'C:\Program Files\FEBioStudio\bin\febio4.exe'; % Ensure this is correct
    
    % 1. Tell MATLAB's system path where the solver is
    setenv('FEBIO_PATH', febio_exe_path);
    loadingOption = 'compression';
    switchVariable = -1;
    dofLabels = {'x', 'y', 'z'};
    
    for k = 1:6
        fprintf('\n--- STARTING LOAD CASE %d/6: %s ---\n', k, labels{k});
        % 0. FORCE MATLAB TO CLOSE ALL DANGLING FILE HANDLES
        fclose('all'); 
        
        % 1. CREATE UNIQUE FILE NAMES FOR THIS SPECIFIC LOAD CASE
        % This prevents FEBio or Windows from locking the file between loops
        currentModelName = sprintf('%s_LC%d', febioFebFileNamePart, k);
        febioFebFileName = fullfile(savePath, [currentModelName, '.feb']);
        febioLogFileName = [currentModelName, '.txt'];
        febioLogFileName_disp = [currentModelName, '_disp_out.txt'];
        febioLogFileName_stress_prin = [currentModelName, '_stress_prin_out.txt'];
        febioLogFileName_stress_full_current = [currentModelName, '_stress_full_out.txt'];
    
        % 1. Default to normal strain (Compression/Tension)
        loadingOption = 'compression';
        switchVariable = -1;
        
        % 2. Switch to shear properties for k = 4, 5, 6
        if k >= 4
            loadingOption = 'shear';
            switchVariable = 1;
        end
    
        % 3. YOUR ORIGINAL LOGIC: Determine the principal axis (1=X, 2=Y, 3=Z)
        loadAxis = mod(k-1, 3) + 1;

        % Dynamically determine the prescribed and support faces based on load axis
        switch loadAxis
            case 1
                bcSupportList = nodeSets.X_min;
                bcPrescribeList = nodeSets.X_max;
            case 2
                bcSupportList = nodeSets.Y_min;
                bcPrescribeList = nodeSets.Y_max;
            case 3
                bcSupportList = nodeSets.Z_min;
                bcPrescribeList = nodeSets.Z_max;
        end
    
        % Reset febio_spec for each run to avoid boundary condition bleed
        febio_spec = febioStructTemplate; 
        febio_spec.ATTR.version = '4.0';
        febio_spec.Module.ATTR.type = 'solid';
        %Control section
        febio_spec.Control.analysis='STATIC';
        febio_spec.Control.time_steps=numTimeSteps;
        febio_spec.Control.step_size=1/numTimeSteps;
        febio_spec.Control.solver.max_refs=max_refs;
        febio_spec.Control.time_stepper.dtmin=dtmin;
        febio_spec.Control.time_stepper.dtmax=dtmax; 
        febio_spec.Control.time_stepper.max_retries=max_retries;
        febio_spec.Control.time_stepper.opt_iter=opt_iter;
        
        % Ensure material E and nu are set as before
        
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
    
        % -> NodeSets
        febio_spec.Mesh.NodeSet{1}.ATTR.name = 'bcSupportList';
        febio_spec.Mesh.NodeSet{1}.VAL = mrow(bcSupportList);
        
        febio_spec.Mesh.NodeSet{2}.ATTR.name = 'bcPrescribeList';
        febio_spec.Mesh.NodeSet{2}.VAL = mrow(bcPrescribeList);
    
        febio_spec.Mesh.NodeSet{3}.ATTR.name = 'allBoundaryNodes';
        febio_spec.Mesh.NodeSet{3}.VAL = mrow(allBoundaryNodes);
           
        %MeshDomains section
        febio_spec.MeshDomains.SolidDomain.ATTR.name=partName1;
        febio_spec.MeshDomains.SolidDomain.ATTR.mat=materialName1;
    
        %% 2. Apply Directional Boundary Conditions
        % Determine which axis we are loading
        % k=1: X, k=2: Y, k=3: Z
        
        
        % 3. Determine the Prescribed Axis
        if k <= 3
            % Compression: Move in the same direction as the face normal
            prescribedAxis = loadAxis;
        else
            % Shear: Move the face in a TRANSVERSE direction (e.g., X-face moves in Y)
            prescribedAxis = mod(loadAxis, 3) + 1; 
        end
        
        %Output section
        % -> log file
        febio_spec.Output.logfile.ATTR.file=febioLogFileName;
        
        febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
        febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
        febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
    
        febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress_prin;
        febio_spec.Output.logfile.element_data{1}.ATTR.data='s1;s2;s3';
        febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
        
        % FIX: Use the unique variable for this load case, but NO fullfile()
        febio_spec.Output.logfile.element_data{2}.ATTR.file=febioLogFileName_stress_full_current;
        febio_spec.Output.logfile.element_data{2}.ATTR.data='sx;sy;sz;sxy;syz;sxz';
        febio_spec.Output.logfile.element_data{2}.ATTR.delim=',';
        
    
        %% 2. Apply Directional Boundary Conditions
        if k <= 3
            % COMPRESSION / TENSION (Normal Strain KUBC)
            prescribedAxis = loadAxis;
            stretchLoad = 1 + switchVariable * appliedStrain;
            displacementMagnitude = (stretchLoad*sampleSize)-sampleSize;
            
            % -> bc{1}: Fix the support face in the load direction ONLY
            febio_spec.Boundary.bc{1}.ATTR.name = 'Fixed_Support';
            febio_spec.Boundary.bc{1}.ATTR.type = 'zero displacement';
            febio_spec.Boundary.bc{1}.ATTR.node_set = 'bcSupportList';
            febio_spec.Boundary.bc{1}.x_dof = double(prescribedAxis == 1);
            febio_spec.Boundary.bc{1}.y_dof = double(prescribedAxis == 2);
            febio_spec.Boundary.bc{1}.z_dof = double(prescribedAxis == 3);
    
            % -> bc{2}: Prescribe the loaded face in the load direction ONLY
            febio_spec.Boundary.bc{2}.ATTR.name = 'Prescribed_Motion';
            febio_spec.Boundary.bc{2}.ATTR.type = 'prescribed displacement';
            febio_spec.Boundary.bc{2}.ATTR.node_set = 'bcPrescribeList';
            febio_spec.Boundary.bc{2}.dof = dofLabels{prescribedAxis}; 
            febio_spec.Boundary.bc{2}.value.ATTR.lc = 1;
            febio_spec.Boundary.bc{2}.value.VAL = displacementMagnitude;
            febio_spec.Boundary.bc{2}.relative = 0;
            
            % -> bc{3}: Constrain lateral Poisson expansion on ALL boundaries
            febio_spec.Boundary.bc{3}.ATTR.name = 'Lateral_Constraint_KUBC';
            febio_spec.Boundary.bc{3}.ATTR.type = 'zero displacement';
            febio_spec.Boundary.bc{3}.ATTR.node_set = 'allBoundaryNodes';
            febio_spec.Boundary.bc{3}.x_dof = double(prescribedAxis ~= 1);
            febio_spec.Boundary.bc{3}.y_dof = double(prescribedAxis ~= 2);
            febio_spec.Boundary.bc{3}.z_dof = double(prescribedAxis ~= 3);
        
        else
            % SHEAR (Kinematic Uniform Boundary Conditions)
            % For simple shear, boundary nodes must be displaced strictly proportional 
            % to their coordinate along the load axis to enforce uniform strain.
            
            prescribedAxis = mod(loadAxis, 3) + 1; 
            coords = {'X', 'Y', 'Z'};
            loadCoord = coords{loadAxis}; % 'X', 'Y', or 'Z'
            strainVal = switchVariable * appliedStrain; 
            
            % Construct the mathematical string for linear displacement (e.g., '0.001*X')
            mathString = sprintf('%g*%s', strainVal, loadCoord);
            
            % -> bc{1}: Prescribe the linear mathematical function to the target DOF 
            % for ALL exterior nodes to enforce the linear shear profile.
            febio_spec.Boundary.bc{1}.ATTR.name = 'Shear_Displacement_KUBC';
            febio_spec.Boundary.bc{1}.ATTR.type = 'prescribed displacement';
            febio_spec.Boundary.bc{1}.ATTR.node_set = 'allBoundaryNodes';
            febio_spec.Boundary.bc{1}.dof = dofLabels{prescribedAxis};
            febio_spec.Boundary.bc{1}.value.ATTR.type = 'math'; % Signals FEBio to evaluate spatially
            febio_spec.Boundary.bc{1}.value.ATTR.lc = 1;
            febio_spec.Boundary.bc{1}.value.VAL = mathString;
            febio_spec.Boundary.bc{1}.relative = 0;
            
            % -> bc{2}: Strictly constrain the remaining two orthogonal DOFs 
            % on all boundary nodes to prevent secondary bending or expansion.
            febio_spec.Boundary.bc{2}.ATTR.name = 'Shear_Orthogonal_Constraint';
            febio_spec.Boundary.bc{2}.ATTR.type = 'zero displacement';
            febio_spec.Boundary.bc{2}.ATTR.node_set = 'allBoundaryNodes';
            febio_spec.Boundary.bc{2}.x_dof = double(prescribedAxis ~= 1);
            febio_spec.Boundary.bc{2}.y_dof = double(prescribedAxis ~= 2);
            febio_spec.Boundary.bc{2}.z_dof = double(prescribedAxis ~= 3);
        end

        %LoadData section
        % -> load_controller
        febio_spec.LoadData.load_controller{1}.ATTR.name='LC1';
        febio_spec.LoadData.load_controller{1}.ATTR.id=1;
        febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
        febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
        febio_spec.LoadData.load_controller{1}.extend='CONSTANT';
        febio_spec.LoadData.load_controller{1}.points.pt.VAL=[0 0; 1 1];
        
        % %Output section
        % % -> log file
        % febio_spec.Output.logfile.ATTR.file=febioLogFileName;
        % febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
        % febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
        % febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
        % 
        % febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress_prin;
        % febio_spec.Output.logfile.element_data{1}.ATTR.data='s1;s2;s3';
        % febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
        % 
        % febio_spec.Output.logfile.element_data{2}.ATTR.file=append(febioLogFileName_stress_full,num2str(k));
        % febio_spec.Output.logfile.element_data{2}.ATTR.data='sx;sy;sz;sxy;syz;sxz';
        % febio_spec.Output.logfile.element_data{2}.ATTR.delim=',';
        
        febio_spec.Output.plotfile.compression=0;
    
        
        febioStruct2xml(febio_spec, febioFebFileName);
        
        febioAnalysis.run_filename=febioFebFileName; %The input file name
        febioAnalysis.run_logname=febioLogFileName; %The name for the log file
        
        % CRITICAL: Point directly to the .exe here to bypass the GUI
        febioAnalysis.FEBioPath = febio_exe_path;
        febioAnalysis.disp_on=1; %Display information on the command window
        febioAnalysis.runMode=runMode;
        febioAnalysis.maxLogCheckTime = 300;
        
        [runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!
    
        %% 4. Data Extraction & Matrix Assembly
        if runFlag == 1
            % Update this line to read the new unique file
            data = importFEBio_logfile(fullfile(savePath, febioLogFileName_stress_full_current), 0, 1);
            stressData = data.data(:, :, end); % Last time step

            % Calculate volume-averaged stress (Voigt vector)
            avgStress = zeros(6,1);
            for i = 1:6
                avgStress(i) = sum(stressData(:, i) .* elemVol) / sampleSize^3;
            end
            
            % Assign column to Stiffness Matrix: C_ij = sigma_i / epsilon_j
            C_homogenized(:, k) = avgStress / (switchVariable * appliedStrain);
        else
            error('Load case %s failed to converge!', labels{k});
        end

        %% UPDATE PROGRESS BAR
        simCounter = simCounter + 1;
        
        % Check if the user manually closed the waitbar to avoid crashing
        if isvalid(hWait)
            waitbar(simCounter / totalSimulations, hWait, ...
                    sprintf('Processing Seed %d, Weight %.1f | Completed: %d / %d', ...
                    currentSeed, weight1, simCounter, totalSimulations));
        end
    end
    
    results(:,:,loop+1) = C_homogenized;

    %% 5. CSV FAIL-SAFE LOGGING
    % Force perfect symmetry to extract the 9 independent Voigt components
    C_sym = (C_homogenized + C_homogenized') / 2;
    
    fid = fopen(csvFileName, 'a');
    fprintf(fid, '%d,%d,%.2f,%.2f,%.1f,%.1f,%.1f,%.2f,%.1f,%.1f,%.1f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n', ...
        currentSeed, loop, weight1, ...
        rho1, theta1(1), theta1(2), theta1(3), ...
        rho2, theta2(1), theta2(2), theta2(3), ...
        C_sym(1,1), C_sym(1,2), C_sym(1,3), ...
        C_sym(2,2), C_sym(2,3), C_sym(3,3), ...
        C_sym(4,4), C_sym(5,5), C_sym(6,6));
    fclose(fid);
    
    fprintf('--> Loop %d safely appended to logging ledger.\n', loop);
    
    % Memory management: flush heavy variables before the next interpolation
    clear graded_GRF X Y Z F V M is_inside_domain D C cc fc Fb Cb CE E;
    java.lang.System.gc();
end

%%
figure()
plot(squeeze(results(1,1,:)))

% --- THE INITIALIZATION & PURGE FIX ---
% Clear any lingering variables to prevent network bleed
load dlnetI

% 1. Load the normalization parameters
load('normalizationParams.mat', 'normParams');

nLoops = 11;
C_arrays = zeros(nLoops, 9);
inputArrays = zeros(nLoops, 4);

for loop = 1:nLoops
    % Extract the 6x6 apparent stiffness matrix computed during this iteration
    C = results(:, :, loop);
    
    % 2. ENFORCE SYMMETRY
    % FEA homogenization often yields slight numerical asymmetry. 
    % We force perfect symmetry before extracting the orthotropic components.
    C_sym = (C + C') / 2;
    
    % Extract the 9 independent Voigt components from the symmetric matrix
    raw_C_vector = [C_sym(1,1), C_sym(1,2), C_sym(1,3), ...
                    C_sym(2,2), C_sym(2,3), C_sym(3,3), ...
                    C_sym(4,4), C_sym(5,5), C_sym(6,6)];
                
    % 3. BASE MATERIAL SCALING
    raw_C_vector = raw_C_vector / E_youngs; 
    
    % 4. NEURAL NETWORK NORMALIZATION (Simple Min-Max)
    C_min = normParams.output.min; 
    C_max = normParams.output.max;
    normalized_C_vector = (raw_C_vector - C_min) ./ (C_max - C_min);
    
    % Clamp to [0, 1] bounds to handle minor FEA floating-point overshoots
    normalized_C_vector = min(max(normalized_C_vector, 0), 1);
    
    C_arrays(loop, :) = normalized_C_vector;
    
    % 5. EVALUATE INVERSE NEURAL NETWORK
    % Feed the 9x1 normalized stiffness profile vector into the network
    nn_output = forward(dlnetI, dlarray(C_arrays(loop, :)', 'CB'));
    nn_output = extractdata(nn_output)'; % Pull vector out of the dlarray wrapper
    
    % 6. DENORMALIZE INVERSE OUTPUTS TO PHYSICAL SPECS (Simple Min-Max)
    in_min = normParams.input.min;
    in_max = normParams.input.max;
    
    % Vectorized denormalization for all 4 components at once:
    % Index 1: Density (rho)
    % Indices 2-4: Conical half angles (theta)
    inputArrays(loop, :) = nn_output .* (in_max - in_min) + in_min;
end

figure()
plot(squeeze(inputArrays(:,1)), 'LineWidth', 2)
title('Homogenized Apparent Relative Density \rho Across Interface Loop');
ylabel('Density \rho');
xlabel('Loop Iteration');
grid on;

% Verification output
fprintf('\n--- Verification Sample (Loop Iteration %d) ---\n', nLoops);
disp('Normalized Stiffness Profile Vector Input:');
disp(C_arrays(end, :));
disp('Recovered Physical Parameters [Density, Theta1, Theta2, Theta3]:');
disp(inputArrays(end, :));

%% Plotting Conical Half-Angles (\theta) Across the Interface Loop
figure('Name', 'Spinodoid Cone Angles - Combined', 'Color', 'w');
hold on;

% inputArrays columns: 1=rho, 2=theta1, 3=theta2, 4=theta3
plot(1:nLoops, inputArrays(:, 2), '-o', 'LineWidth', 2, 'DisplayName', '\theta_1 (X-axis)');
plot(1:nLoops, inputArrays(:, 3), '-s', 'LineWidth', 2, 'DisplayName', '\theta_2 (Y-axis)');
plot(1:nLoops, inputArrays(:, 4), '-^', 'LineWidth', 2, 'DisplayName', '\theta_3 (Z-axis)');

hold off;
title('Recovered Conical Half-Angles (\theta) Across Interface');
ylabel('Angle (Degrees)');
xlabel('Loop Iteration');
ylim([0 90]); % Physical limits of the spinodoid cone angles
legend('Location', 'best');
grid on;

end

%% CLOSE PROGRESS BAR
if isvalid(hWait)
    close(hWait);
end

%% =========================================================================
% POST-PROCESSING: CSV INGESTION, AVERAGING, & HYPOTHESIS TESTING
% =========================================================================
fprintf('\n=== READING FAIL-SAFE LEDGER ===\n');

% Load Network and Normalization Parameters
load dlnetI
load dlnetF
load('normalizationParams.mat', 'normParams');

% 1. Read the CSV Data
defaultFolder = fileparts(mfilename('fullpath'));
csvPath = fullfile(defaultFolder, 'data', 'temp');

if ~isfile(csvPath)
    error('CSV file not found! Check the path: %s', csvPath);
end
data = readtable(csvPath);

% Initial parameters
E_youngs = 1; % Material properties used for base scaling
nLoops = 11;
seeds = [1, 2, 3];

% Preallocate E1 tracking arrays
E1_seeds = zeros(length(seeds), nLoops);
E1_avg = zeros(1, nLoops);
inputArrays_avg = zeros(nLoops, 4);

% 2. Process Data Loop by Loop
for l = 0:(nLoops-1)
    % Extract data for this specific loop iteration across all seeds
    loopData = data(data.LoopIndex == l, :);
    
    % Initialize a zero matrix for averaging the stiffness
    C_sum = zeros(6,6);
    
    % --- A. Calculate E1 for Individual Seeds ---
    for s = 1:length(seeds)
        seedRow = loopData(loopData.Seed == seeds(s), :);
        
        if ~isempty(seedRow)
            % Rebuild the symmetric matrix from the Voigt components
            C_s = zeros(6,6);
            C_s(1,1)=seedRow.C11; C_s(1,2)=seedRow.C12; C_s(1,3)=seedRow.C13;
            C_s(2,1)=seedRow.C12; C_s(2,2)=seedRow.C22; C_s(2,3)=seedRow.C23;
            C_s(3,1)=seedRow.C13; C_s(3,2)=seedRow.C23; C_s(3,3)=seedRow.C33;
            C_s(4,4)=seedRow.C44; C_s(5,5)=seedRow.C55; C_s(6,6)=seedRow.C66;
            
            C_s_sym = (C_s + C_s') / 2; 
            S_s = inv(C_s_sym);         
            E1_seeds(s, l+1) = 1 / S_s(1,1);
            
            % Accumulate for averaging
            C_sum = C_sum + C_s_sym;
        end
    end

    % --- B. Calculate E1 for the Averaged Stiffness ---
    C_avg_sym = C_sum / length(seeds);
    S_avg = inv(C_avg_sym);
    E1_avg(l+1) = 1 / S_avg(1,1);
    
    % --- C. Feed Averaged Stiffness to Inverse NN ---
    raw_C_vector = [C_avg_sym(1,1), C_avg_sym(1,2), C_avg_sym(1,3), ...
                    C_avg_sym(2,2), C_avg_sym(2,3), C_avg_sym(3,3), ...
                    C_avg_sym(4,4), C_avg_sym(5,5), C_avg_sym(6,6)];
                
    raw_C_vector = raw_C_vector / E_youngs; 
    
    normalized_C_vector = (raw_C_vector - normParams.output.min) ./ (normParams.output.max - normParams.output.min);
    normalized_C_vector = min(max(normalized_C_vector, 0), 1);
    
    nn_output = forward(dlnetI, dlarray(normalized_C_vector', 'CB'));
    nn_output = extractdata(nn_output)'; 
    
    inputArrays_avg(l+1, :) = nn_output .* (normParams.input.max - normParams.input.min) + normParams.input.min;
end

% --- D. Evaluate the Hypothesis via Forward NN ---
hyp_raw = [0.6, 15, 15, 15]; % [rho, theta1, theta2, theta3]
hyp_norm = (hyp_raw - normParams.input.min) ./ (normParams.input.max - normParams.input.min);

dl_hyp = dlarray(single(hyp_norm'), 'CB');
prop_hyp_norm = extractdata(forward(dlnetF, dl_hyp))';
prop_hyp = prop_hyp_norm .* (normParams.output.max - normParams.output.min) + normParams.output.min;

% Rebuild C and get E1 for the hypothesis
C_hyp = zeros(6,6);
C_hyp(1,1) = prop_hyp(1); C_hyp(1,2) = prop_hyp(2); C_hyp(1,3) = prop_hyp(3);
C_hyp(2,1) = prop_hyp(2); C_hyp(2,2) = prop_hyp(4); C_hyp(2,3) = prop_hyp(5);
C_hyp(3,1) = prop_hyp(3); C_hyp(3,2) = prop_hyp(5); C_hyp(3,3) = prop_hyp(6);
C_hyp(4,4) = prop_hyp(7); C_hyp(5,5) = prop_hyp(8); C_hyp(6,6) = prop_hyp(9);

S_hyp = inv(C_hyp);
E1_hyp = 1 / S_hyp(1,1);

%% 3. PLOTTING
% Plot 1: E1 Modulus Trajectory vs Hypothesis
figure('Name', 'E1 Modulus Tracking', 'Color', 'w', 'Position', [100, 100, 800, 500]);
hold on;
plot(1:nLoops, E1_seeds(1,:), 'Color', [0.7 0.7 0.7], 'LineStyle', '--', 'DisplayName', 'Seed 1');
plot(1:nLoops, E1_seeds(2,:), 'Color', [0.7 0.7 0.7], 'LineStyle', ':', 'DisplayName', 'Seed 2');
plot(1:nLoops, E1_seeds(3,:), 'Color', [0.7 0.7 0.7], 'LineStyle', '-.', 'DisplayName', 'Seed 3');

plot(1:nLoops, E1_avg, 'k-o', 'LineWidth', 2.5, 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Averaged E_1');
yline(E1_hyp, 'r-', 'LineWidth', 2, 'DisplayName', 'Hypothesis [15,15,15]');

title('Directional Young''s Modulus (E_1) Across Transition');
xlabel('Loop Iteration');
ylabel('Stiffness E_1 (MPa)');
legend('Location', 'best');
grid on; hold off;

% Plot 2: Inverse NN Prediction based on Averaged Stiffness
figure('Name', 'Inverse NN - Averaged', 'Color', 'w', 'Position', [950, 100, 800, 500]);
hold on;
plot(1:nLoops, inputArrays_avg(:, 2), '-o', 'LineWidth', 2, 'DisplayName', '\theta_1 (X-axis)');
plot(1:nLoops, inputArrays_avg(:, 3), '-s', 'LineWidth', 2, 'DisplayName', '\theta_2 (Y-axis)');
plot(1:nLoops, inputArrays_avg(:, 4), '-^', 'LineWidth', 2, 'DisplayName', '\theta_3 (Z-axis)');

title('Averaged Recovered Conical Half-Angles (\theta)');
ylabel('Angle (Degrees)');
xlabel('Loop Iteration');
ylim([0 90]);
legend('Location', 'best');
grid on; hold off;

function [X_q,Y_q,Z_q,weights_i,limits] = calcWeights(center,res,resolution,neighborhoodRadius,kappa)
    % 1. Get Element Center
    % 2. Define Local Bounding Box (Logic moved from calcWeights)
    % We only process voxels within ±1.5 units of the center
    cx = center(1);cy = center(2);cz=center(3);
    mins = max(1, floor((center - neighborhoodRadius) * res) + 1)';
    maxs = min(resolution, ceil((center + neighborhoodRadius) * res))';
    limits = [mins,maxs];
    
    % 3. Create Local Coordinate Grid
    % These points align perfectly with the global grid but only cover the neighborhood
    xv_local = (mins(1) - 0.5 : maxs(1) - 0.5) / res;
    yv_local = (mins(2) - 0.5 : maxs(2) - 0.5) / res;
    zv_local = (mins(3) - 0.5 : maxs(3) - 0.5) / res;
    [X_loc, Y_loc, Z_loc] = ndgrid(xv_local, yv_local, zv_local);

    
    % Sampling coordinates: Global pos - (Center - 0.5) maps center to 0.5
    [X_q, Y_q, Z_q] = ndgrid(xv_local - center(1) + 0.5, yv_local - center(2) + 0.5, zv_local - center(3) + 0.5);
    % 5. Calculate Weights for the local grid
    % Calculate Distance Deltas w.r.t. element center
    dX = X_loc - cx; dY = Y_loc - cy; dZ = Z_loc - cz;
    weights_i = exp(- ((2*abs(dX)).^kappa + (2*abs(dY)).^kappa + (2*abs(dZ)).^kappa));
end

function [GRF, X, Y, Z, levelset] = genGRF(inputStruct)
    % Parse input
    
    %Create default structure
    defaultInputStruct.isocap=true; % option to cap the isosurface
    defaultInputStruct.domainSize=1; % domain size: scalar for cube, 
                                     % or [a,b,c] for cuboid of size [a,b,c]
    defaultInputStruct.resolution=100; % resolution for sampling GRF scalar for cube
                                       % or [ra,rb,rc] for cuboid with ra*rb*rc points
    defaultInputStruct.waveNumber=15*pi; % GRF wave number
    defaultInputStruct.numWaves=1000; % number of waves in GRF
    defaultInputStruct.relativeDensity=0.5; % relative density: between [0.3,1]
    defaultInputStruct.thetas=[15 15 0]; % conical half angles (in degrees) 
    %                                     along xyz axes for controlling the 
    %                                     anisotropy. Note: each entry must be 
    %                                     either 0 or between [15,90] degrees.
    defaultInputStruct.R = eye(3); % Rotate the GRF, R must be SO(3)
    defaultInputStruct.ignoreChecks = false; % Ignore checks on parameters if true (not advised)
    defaultInputStruct.trimDomainFunction = @doNothing; % A function handle handle (preceeded by @)
    %                                                     that takes in coordinates
    %                                                     x,y,z and returns true if
    %                                                     the coordinates lie inside
    %                                                     a desirable domain and
    %                                                     false if outside. The
    %                                                     default function
    %                                                     doNothing() is defined at
    %                                                     the end of this file. For
    %                                                     example, if you want to
    %                                                     generate a spinodoid
    %                                                     sample inside a (1/8th)
    %                                                     unit sphere centered at
    %                                                     (0,0,0), then
    %                                                     trimDomainFunction(x,y,z)
    %                                                     should output true if 
    %                                                     (x.^2+y.^2+z.^2 <= 1^2).
    defaultInputStruct.patchDomain.F=[];
    defaultInputStruct.patchDomain.V=[];
    
    %Complete input with default if incomplete
    [inputStruct]=structComplete(inputStruct,defaultInputStruct,1); %Complement provided with default if missing or empty
    
    %Get parameters from input structure
    isocap = inputStruct.isocap; % option to cap the isosurface
    domainSize = inputStruct.domainSize; % domain size
    resolution = inputStruct.resolution; % resolution for sampling GRF
    waveNumber = inputStruct.waveNumber; % GRF wave number
    numWaves = inputStruct.numWaves; % number of waves in GRF
    relativeDensity = inputStruct.relativeDensity; % relative density: [0.3,1]
    thetas= inputStruct.thetas; % conical half angles (in degrees)
    R = inputStruct.R; % Rotate the GRF, R must be SO(3)
    ignoreChecks = inputStruct.ignoreChecks; %Ignore checks on parameter values
    trimDomainFunction = inputStruct.trimDomainFunction; %trimDomainFunction handle
    F_patchDomain=inputStruct.patchDomain.F;
    V_patchDomain=inputStruct.patchDomain.V;
    
    %% Input checks
    if(ignoreChecks)
        warning(['Ignoring all checks on parameter values. ',...
            'Unreasonable parameters may give wrong topologies.']);
    else
        if((relativeDensity<0.0) || (relativeDensity>1.0))
            error('relativeDensity must be between [0.3,1]')
        end
        if((relativeDensity<0.3))
            warning('Relative density too low, it may produce discontinuous domains')
        end
        if((any(thetas<0)) || (any(thetas>90)))
            error('thetas must be either 0 or between [15,90] degrees')
        end
        for i=1:3
            if((thetas(i)>0) && (thetas(i)<15))
                warning('thetas must be either 0 or between [15,90] degrees to ensure continuous domains')
            end
        end
    
        if(size(R,1)~=size(R,2))
            error('Rotation matrix is not square')
        end
        if(size(R,1)~=3)
            error('Rotation matrix must be 3x3 in size')
        end
        if(abs(det(R)-1)>1e-8)
            error('Rotation matrix: det(R)~=1')
        end
        if(norm(R'*R-eye(3))>1e-8)
            error('Rotation matrix is not orthogonal')
        end
        
        if(any(size(resolution) ~= size(domainSize)))
            error('domainSize and resolution must be of same size')
        end
        
        if((length(domainSize) == 1) || (length(domainSize) == 3))
        else
            error('domainSize must be scalar or 1x3 vector')
        end
    end
    
    %% Define rotated axes
    axes1 = [1,0,0];
    axes2 = [0,1,0];
    axes3 = [0,0,1];
    
    axes1 = (R*axes1')';
    axes2 = (R*axes2')';
    axes3 = (R*axes3')';
    
    %% Generate wave directions for GRF
    %array of all wave directions
    
    waveDirections = zeros(numWaves,3);
    for i=1:numWaves
        flag = true; %keep trying until candidate wave vector is found
        while(flag)
            % generate isotropically sampled candidate wave
            candidate = randn(1,3);
            candidate = candidate/norm(candidate);
            % check for allowed wave vector directions
            % angle along first axis
            angle1 = min(...
                acosd(dot(candidate,axes1)),...
                acosd(dot(candidate,-axes1)));
            % angle along second axis
            angle2 = min(...
                acosd(dot(candidate,axes2)),...
                acosd(dot(candidate,-axes2)));
            % angle along third axis
            angle3 = min(...
                acosd(dot(candidate,axes3)),...
                acosd(dot(candidate,-axes3)));
            % check
            if(any([angle1,angle2,angle3]<thetas))
                flag = false;
                break;
            end
        end
        waveDirections(i,:)=candidate;
    end
    
    %% Generate wave phase angles for GRF
    wavePhases = rand_angle([numWaves,1]); %2*pi*rand(numWaves,1);
    
    %% Discretize the domain
    if (length(domainSize) == 1)
        discretization = linspace(0,domainSize,resolution);
        [X,Y,Z] = ndgrid(discretization,discretization,discretization);
    else
        discretizationX = linspace(0,domainSize(1),resolution(1));
        discretizationY = linspace(0,domainSize(2),resolution(2));
        discretizationZ = linspace(0,domainSize(3),resolution(3));
        [X,Y,Z] = ndgrid(discretizationX,discretizationY,discretizationZ);
    end
    
    
    %% Evaluate GRF on sampling points
    
    GRF = zeros(size(X));
    for i=1:numWaves
        dotProduct = waveDirections(i,1)*X ...
                   + waveDirections(i,2)*Y ...
                   + waveDirections(i,3)*Z;
                 
        GRF = GRF+sqrt(2/numWaves)*cos(dotProduct*waveNumber+wavePhases(i));
    end
    
    %% Trim domain
    
    % Create trim logic using domain function 
    is_inside_domain = trimDomainFunction(X,Y,Z);
    
    % Create/adjust trim logic using patch data
    if ~isempty(F_patchDomain)
        voxelSize=domainSize./resolution;
        imOrigin=[0 0 0];
        imSiz=size(GRF);
        M=patch2Im(F_patchDomain,V_patchDomain,ones(size(F_patchDomain,1),1),voxelSize,imOrigin,imSiz);
        M=permute(M,[2 1 3]);
        is_inside_domain = is_inside_domain & M==1;
    end
    
    if nnz(~is_inside_domain)>0 %If any domain trimming is needed
        % GRF = GRF + (1-is_inside_domain)*1e8;
    
        % D = bwdist(is_inside_domain,'Euclidean');
        % D=D./2;
        % D(D>1)=1;
        % D=1-D;
        % GRF = GRF.*D;
    
        GRF(~is_inside_domain) = 2*max(GRF(:));
    end
    
    %% Apply levelset
    levelset = sqrt(2)*erfinv(2*relativeDensity-1);
end

function [var] = doNothing(x,y,z)
    var = true(size(x));
end

%%
%
% <<gibbVerySmall.gif>>
%
% _*GIBBON*_
% <www.gibboncode.org>
%
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>

%%
% _*GIBBON footer text*_
%
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
%
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
%
% Copyright (C) 2006-2020 Kevin Mattheus Moerman
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
