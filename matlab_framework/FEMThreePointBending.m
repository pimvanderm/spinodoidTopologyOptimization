tr = stlread('TopOpt_three_point_bending_half_88x20x20_v03_k15_2026-06-03_153353.stl');
% 1. RENDERING SETUP
settings.nelx = 88; settings.nely = 20; settings.nelz = 20;
settings.domainSize = [settings.nelx, settings.nely, settings.nelz]; % 1 unit per element

numberOfWaves = 15;
settings.numberOfWaves = numberOfWaves;
settings.res = 5 * numberOfWaves / settings.nelz; 
settings.k_wave = numberOfWaves * 1 / settings.nelz * 2 * pi;
settings.writeSTL = 0;

settings.problemType = 'three_point_bending_half';
settings.volfrac = 0.3;




clear inputStruct
clear meshStruct
clear optionStruct

%% Plot settings
fontSize=20;
faceAlpha1=0.8;
markerSize=40;
lineWidth1=3;
lineWidth2=4;
markerSize1=25;

%Material parameter set
E_youngs=2800; %Youngs modulus in GPa
nu=0.33; %Poisson's ratio
% mu=E_youngs/3;

%FEA control settings
numTimeSteps=1; %Number of time steps desired
max_refs=15; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=6; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size

runMode='external';% 'internal' or 'external'

%% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');

%%

% ---> CRITICAL FIX: Create the folder if it does not exist <---
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

%%

% Defining file names
febioFebFileNamePart='3PB';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_stress_prin=[febioFebFileNamePart,'_stress_prin_out.txt']; %Log file name for exporting stress
febioLogFileName_stress_full=[febioFebFileNamePart,'_stress_full_out.txt']; %Log file name for exporting stress

%% Remesh using geomgram
pointSpacing=20/75;
optionStruct =struct();
optionStruct.pointSpacing = pointSpacing;
overSampleRatio=1;
%numStepsLevelset=ceil(overSampleRatio.*(sampleSize./pointSpacing)); %Number of voxel steps across period for image data (roughly number of points on mesh period)


% optionStruct.max_dist=0;
f = tr.ConnectivityList;
v = tr.Points;

[F,VPreMesh]=ggremesh(f,v,optionStruct);

C=zeros(size(F,1),1);

%%
% Visualizing geometry

cFigure; hold on;
title('Geogram remeshed','FontSize',fontSize);

gpatch(F,VPreMesh,'w','k',1);

axisGeom(gca,fontSize);
camlight headlight;
inner_pt = getInnerPoint(F,VPreMesh);
plot3(inner_pt(1), inner_pt(2), inner_pt(3), 'r.', 'MarkerSize', 50);
drawnow;


%% Tetrahedral meshing using tetgen (see also |runTetGen|)

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


% Create tetgen input structure
inputStruct.stringOpt='-pq2AY';
inputStruct.Faces=F;
inputStruct.Nodes=VPreMesh;
inputStruct.holePoints=[];
inputStruct.faceBoundaryMarker=C; %Face boundary markers
inputStruct.regionPoints=innerPoint; %region points
inputStruct.regionA=10*tetVolMeanEst(F,VPreMesh);
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

%% Coordinates of the BCs
tolDir = 1;
C_vertex = zeros(size(V,1),1);
dimensions = max(v);
xmin = min(v(:,1));

% Keep your existing physical locations for the loading/supports
constraint1 = [(dimensions(1)+xmin)/2-2*dimensions(3) 0];
constraint2 = [(dimensions(1)+xmin)/2+2*dimensions(3) 0];
load1 = [(dimensions(1)+xmin)/2 dimensions(3)];

logic1 = vecnorm(V(:,[1 3])-constraint1,2,2)<tolDir; % Bottom Left Roller Support
logic2 = vecnorm(V(:,[1 3])-constraint2,2,2)<tolDir; % Bottom Right Roller Support
logic3 = vecnorm(V(:,[1 3])-load1,2,2)<tolDir;       % Top Center Loading Line

C_vertex(logic1) = 1;
C_vertex(logic2) = 1;
C_vertex(logic3) = 2;

% Define boundary condition lists based on C_vertex
bcSupportList   = find(C_vertex == 1); % Both roller support lines combined
bcPrescribeList = find(C_vertex == 2); % The entire top loading line

% --- NEW CRITICAL FIXES FOR RIGID BODY MOTION STABILIZATION ---
% 1. Find the exact mid-point of the top loading line along the Y-axis to lock Y-translation
y_center = (max(V(:,2)) + min(V(:,2))) / 2;
[~, center_idx_within_load] = min(abs(V(bcPrescribeList, 2) - y_center));
bcCenterNode = bcPrescribeList(center_idx_within_load); % Single node at the center of the beam
%% 
% Visualizing boundary conditions. Markers plotted on the semi-transparent
% model denote the nodes in the various boundary condition lists.

hf=cFigure;
title('Boundary conditions','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(Fb,V,'w','k',1);

hl(1)=plotV(V(bcSupportList,:),'r.','MarkerSize',markerSize);
hl(2)=plotV(V(bcPrescribeList,:),'k.','MarkerSize',markerSize);

legend(hl,{'BC full support','BC z prescribe'});

axisGeom(gca,fontSize);
camlight headlight;
drawnow;


%% Setup FEBio

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


%%

if runFlag==1 %i.e. a succesful run
    
    %%
    % Importing nodal displacements from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),0,1);
    
    %Access data
    N_disp_mat=dataStruct.data; %Displacement
    timeVec=dataStruct.time; %Time
    

    % 1. Calculate the bounding height of your mesh (Assuming Z is the vertical axis)
    % If your part is vertical along the Y-axis, change V(:,3) to V(:,2)
    meshHeight = max(V(:,3)) - min(V(:,3)); 
    
    % 2. Find the absolute maximum displacement magnitude in the final time step
    maxDisp = max(sqrt(sum(N_disp_mat(:,:,end).^2, 2)));
    
    % 3. Calculate the multiplier needed to make the max displacement = 1/10th of the height
    scaleFactor = (meshHeight / 4) / maxDisp;
    
    % 4. Create the scaled deformed coordinate set
    V_DEF = (N_disp_mat * scaleFactor) + repmat(V,[1 1 size(N_disp_mat,3)]);

    %%
    % Importing element stress from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_stress_prin),0,1);
    
    %Access data
    E_stress_mat=dataStruct.data;
    
    E_stress_mat(isnan(E_stress_mat))=0;
    
    %Compute Von Mises
    S_vm_ND = sqrt( 0.5*((E_stress_mat(:,1,:)-E_stress_mat(:,2,:)).^2 + (E_stress_mat(:,2,:)-E_stress_mat(:,3,:)).^2 + (E_stress_mat(:,1,:)-E_stress_mat(:,3,:)).^2));
    
    %Compute volume weighted mean Von Mises stress
    
    elemVol=tetVol(E,V); %Element volumes
    totalVol=sum(elemVol); %Total volume
    S_vm_mean=squeeze(sum(S_vm_ND.*repmat(elemVol,[1 1 size(S_vm_ND,3)]),1)./totalVol); %Mean Von Mises stress
    
    %%
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations
    
    DN_magnitude=sqrt(sum(N_disp_mat(:,:,end).^2,2)); %Current displacement magnitude
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('Displacement magnitude [mm]','Interpreter','Latex')
    hp=gpatch(Fb,V_DEF(:,:,end),DN_magnitude,'k',1); %Add graphics object to animate
    hp.FaceColor='interp';
    
    axisGeom(gca,fontSize);
    colormap(gjet(250)); colorbar;
    caxis([0 max(DN_magnitude)]);
    axis(axisLim(V_DEF)); %Set axis limits statically
    camlight headlight;
    
    % Set up animation features
    animStruct.Time=timeVec; %The time vector
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments
        DN_magnitude=sqrt(sum(N_disp_mat(:,:,qt).^2,2)); %Current displacement magnitude
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),DN_magnitude}; %Property values for to set in order to animate
    end
    anim8(hf,animStruct); %Initiate animation feature
    drawnow;
    
    %%
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations
    
    [CV]=faceToVertexMeasure(E,V,S_vm_ND(:,:,end));
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('$\sigma_{1}$ [MPa]','Interpreter','Latex')
    hp=gpatch(Fb,V_DEF(:,:,end),CV,'k',1); %Add graphics object to animate
    
    hp.FaceColor='interp';
    
    axisGeom(gca,fontSize);
    colormap(gjet(250)); colorbar;
    caxis([min(S_vm_ND(:)) max(S_vm_ND(:))]/4);
    axis(axisLim(V_DEF)); %Set axis limits statically
    camlight headlight;
    
    % Set up animation features
    animStruct.Time=timeVec; %The time vector
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments
        
        [CV]=faceToVertexMeasure(E,V,S_vm_ND(:,:,qt));
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),CV}; %Property values for to set in order to animate
    end
    anim8(hf,animStruct); %Initiate animation feature
    drawnow;
    
end

%% =====================================================================
% NEW POST-PROCESSING: INTERPOLATE FEBIO DISPLACEMENTS TO TOPOPT GRID NODES
% =====================================================================
disp('Mapping high-fidelity FEBio displacements to TopOpt grid...');

% 1. Reconstruct the exact TopOpt grid node coordinates in physical space (mm)
% Grid element dimensions match your setup: 44 x 10 x 10
% Physical macro-beam dimensions are: 88 x 20 x 20 mm
nelx_to = 44; 
nely_to = 10; 
nelz_to = 10;

L_scale = 20 / nelz_to; % Scaling factor based on span length (mm per element length)

nnodx = nelx_to + 1;
nnody = nely_to + 1;
nnodz = nelz_to + 1;
nnode_to = nnodx * nnody * nnodz;

% Generate the orderly mesh grid of nodes matching your TopOpt ordering
% Elements vary fastest in z, then y, then x (matches setupFEM ndgrid ordering)
[GridZ, GridY, GridX] = ndgrid(0:nelz_to, 0:nely_to, 0:nelx_to);

% Scale grid indices to actual physical millimeter coordinates
V_topopt_nodes = [GridX(:), GridY(:), GridZ(:)] * L_scale; 

% 2. Extract final step displacement components from FEBio results
% N_disp_mat format: [Nodes x Dimensions(3) x TimeSteps]
U_febio_final = N_disp_mat(:, :, end); 

% 3. Calculate centroids of the FEBio tetrahedral elements for quick spatial lookup
% E is the [nele_tet x 4] connectivity matrix from TetGen
tet_centroids = zeros(size(E, 1), 3);
for i = 1:4
    tet_centroids = tet_centroids + V(E(:, i), :);
end
tet_centroids = tet_centroids / 4;

% Build a k-d tree for ultra-fast spatial querying of element centroids
Mdl = triangulation(E, V);

% Initialize the TopOpt-formatted full structural displacement array (3*N x 1)
% Matching the exact DOF structure of columnOrientation: [Ux1; Uy1; Uz1; Ux2; Uy2; Uz2...]
U_interpolated = zeros(3 * nnode_to, 1);

% 4. Loop over every node in the TopOpt grid and interpolate displacement
% We find the closest solid tetrahedral elements to prevent void-space dropouts
num_neighbors_to_average = 5; 

for n = 1:nnode_to
    pt = V_topopt_nodes(n, :);

    % Find the closest tetrahedral node from the high-fidelity mesh directly
    % using native point-location nearestNeighbor mapping
    idx_node = nearestNeighbor(Mdl, pt);

    % Extract the displacement vector from that specific node point location
    disp_mapped = U_febio_final(idx_node, :);

    % Slot components back into the 3*N mapping matrix format (Voigt interleaved)
    U_interpolated(3*n - 2) = disp_mapped(1); % Ux
    U_interpolated(3*n - 1) = disp_mapped(2); % Uy
    U_interpolated(3*n)     = disp_mapped(3); % Uz
end

disp('Mapping complete. Handing over vector to visualizeFEM.');

% 5. Package a mock FEM structure matching columnOrientation specs for visualization
mock_FEM.nelx = nelx_to;
mock_FEM.nely = nely_to;
mock_FEM.nelz = nelz_to;
mock_FEM.nele = nelx_to * nely_to * nelz_to;
mock_FEM.ndof = 3 * nnode_to;

% --- DIRECT CALCULATION OF TOPOPT EDOFMAT ---
% Reconstruct the index mapping configuration matching setupFEM lines 105-135
nodenrs_mock = reshape(1:nnode_to, nnodz, nnody, nnodx);
anchorNodes_mock = nodenrs_mock(1:nelz_to, 1:nely_to, 1:nelx_to);
anchorNodes_mock = anchorNodes_mock(:); % Flatten to column vector

% Strides matching the node arrangement layout matrix
sz_m = 1;           
sy_m = nnodz;       
sx_m = nnodz * nnody; 

% Hex node structural offsets map sequence matching your solver setup
nodeOffsets_mock = [0, sx_m, sy_m+sx_m, sy_m, sz_m, sz_m+sx_m, sz_m+sy_m+sx_m, sz_m+sy_m];

% Allocate and calculate the mock degrees of freedom matrix
mock_edofMat = zeros(mock_FEM.nele, 24);
for n = 1:8
    globalNode = anchorNodes_mock + nodeOffsets_mock(n);
    mock_edofMat(:, 3*n-2) = 3*globalNode - 2;     % x-DOF index mapping link
    mock_edofMat(:, 3*n-1) = 3*globalNode - 1;     % y-DOF index mapping link
    mock_edofMat(:, 3*n  ) = 3*globalNode;          % z-DOF index mapping link
end

% Store the matrix inside the structure handled by the rendering plot loops
mock_FEM.edofMat = mock_edofMat;
mock_FEM.fixedDofs = [];        
mock_FEM.F = sparse(mock_FEM.ndof, 1);

% Call your original visualization function to inspect FEBio results on the grid topology
% scaled to 1 / max displacement magnitude to ensure a visible rendering profile
%visualizeFEM(mock_FEM, U_interpolated, chi_final, 1/max(abs(U_interpolated)));
%sgtitle('High-Fidelity FEBio Displacement Field Mapped back to TopOpt Voxel Mesh Grid');

%% =====================================================================
% NEW POST-PROCESSING: PLOT DISPLACEMENT DIFFERENCE FIELD (DEVIATION)
% =====================================================================
disp('Calculating spatial displacement deviation between TopOpt and FEBio...');

% 1. Ensure the raw TopOpt displacement vector is physically scaled (mm)
% If U_final hasn't been scaled by L_scale yet, do it here to ensure a 1:1 match
L_scale_factor = 0.5; % 2 mm/element scale factor
U_topopt_scaled = U_final * L_scale_factor; 

% 2. Calculate the component-wise error vector at every single degree of freedom
% Error = High-Fidelity Tet FEM - Homogenized Hex TopOpt
U_error_field = U_interpolated - U_topopt_scaled;

% 3. Calculate the absolute Euclidean displacement error magnitude at each node
% Number of nodes matches your grid specs (45 * 11 * 11)
num_nodes_total = length(U_error_field) / 3;
error_magnitude_per_node = zeros(num_nodes_total, 1);

for n = 1:num_nodes_total
    ux_err = U_error_field(3*n - 2);
    uy_err = U_error_field(3*n - 1);
    uz_err = U_error_field(3*n);

    % Vector magnitude of the error tracking mismatch
    error_magnitude_per_node(n) = sqrt(ux_err^2 + uy_err^2 + uz_err^2);
end

% Print global mismatch benchmarks to the command window
fprintf('\n==================================================\n');
fprintf('DISPLACEMENT DISCREPANCY ANALYSIS (FEBio vs TopOpt)\n');
fprintf('==================================================\n');
fprintf('Maximum Nodal Displacement Deviation:   %0.4f mm\n', max(error_magnitude_per_node));
fprintf('Mean Field Displacement Deviation:      %0.4f mm\n', mean(error_magnitude_per_node));
fprintf('==================================================\n');

% 4. Package the mismatch data back into a structured Voigt array vector
% We reuse the interleaved Z-component field space to force the 
% visualization function color map engine to render the scalar error field directly.
U_error_visualization = zeros(size(U_error_field));
for n = 1:num_nodes_total
    % Keep the real directional deformation tracking vectors intact so 
    % the beam model bends correctly on screen...
    U_error_visualization(3*n-2) = U_interpolated(3*n-2); 
    U_error_visualization(3*n-1) = U_interpolated(3*n-1);

    % ...but override the Z-displacement channel with our scalar error magnitude.
    % This forces the visualization mesh color plots to display deviation directly.
    U_error_visualization(3*n) = error_magnitude_per_node(n);
end

% 5. Execute rendering via your original layout engine
% Using an absolute baseline scale factor of 1 to ensure standard frame tracking bounds
visualizeFEM(mock_FEM, U_topopt_scaled, U_interpolated, U_error_field, chi_final, 20)

% % Adjust color bar and title fields to reflect deviation metrics instead of normal tracking
% h_fig = gcf;
% h_axes = findobj(h_fig, 'Type', 'axes');
% for a = 1:length(h_axes)
%     h_cb = findobj(h_fig, 'Type', 'ColorBar');
%     if ~isempty(h_cb)
%         h_cb.Label.String = 'Absolute Displacement Mismatch [mm]';
%     end
% end
% 
% colormap(h_fig, 'hot'); % Switch to a 'hot' heatmap to accentuate high-error zones
sgtitle('Spatial Deviation Map: Unstructured Tet Mesh vs Homogenized Voxel Grid');

function visualizeFEM(FEM, U_topopt, U_febio, U_error, chi, scale)
% VISUALIZEFEM Visualize TopOpt, FEBio, and Error fields side-by-side
% Signature: (FEM, U_topopt, U_febio, U_error, chi, scale)

if nargin < 6
    scale = 1.0;
end
if nargin < 5
    chi = [];
end

% -------------------------------------------------------------------------
% EXTRACT NODE COORDINATES
% -------------------------------------------------------------------------
nnodx = FEM.nelx + 1;
nnody = FEM.nely + 1;
nnodz = FEM.nelz + 1;
nnode = nnodx * nnody * nnodz;

[IZ, IY, IX] = ndgrid(1:nnodz, 1:nnody, 1:nnodx);
X0 = IX(:) - 1;   
Y0 = IY(:) - 1;   
Z0 = IZ(:) - 1;   

% Extract and calculate nodal displacement magnitudes for all three fields
Umag_topopt = sqrt(U_topopt(1:3:end).^2 + U_topopt(2:3:end).^2 + U_topopt(3:3:end).^2);
Umag_febio  = sqrt(U_febio(1:3:end).^2  + U_febio(2:3:end).^2  + U_febio(3:3:end).^2);
Umag_error  = sqrt(U_error(1:3:end).^2  + U_error(2:3:end).^2  + U_error(3:3:end).^2);

% Calculate deformed coordinate sets
X1_to = X0 + scale * U_topopt(1:3:end);
Y1_to = Y0 + scale * U_topopt(2:3:end);
Z1_to = Z0 + scale * U_topopt(3:3:end);

X1_fb = X0 + scale * U_febio(1:3:end);
Y1_fb = Y0 + scale * U_febio(2:3:end);
Z1_fb = Z0 + scale * U_febio(3:3:end);

% -------------------------------------------------------------------------
% IDENTIFY BOUNDARY CONDITION NODES & LOADS
% -------------------------------------------------------------------------
fixedNodes = unique(ceil(FEM.fixedDofs / 3));
[loadDofs, ~, loadVals] = find(FEM.F);
loadNodes = ceil(loadDofs / 3);
loadDirs  = mod(loadDofs - 1, 3) + 1;   
arrowScale = 0.4 * max(FEM.nelx,FEM.nely);

% -------------------------------------------------------------------------
% SETUP FIGURE
% -------------------------------------------------------------------------
figure('Name', 'Displacement Comparison', 'Color', 'white', 'Position', [50, 150, 1800, 500]);

% =========================================================================
% SUBPLOT 1: TopOpt Displacement
% =========================================================================
ax1 = subplot(1, 3, 1);
hold on; axis equal; grid on; box on;
title('TopOpt Homogenized Displacement', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('x'); ylabel('y'); zlabel('z');

drawMeshEdgesColored(FEM, X1_to, Y1_to, Z1_to, Umag_topopt, chi);

% Replot BCs
scatter3(X1_to(fixedNodes), Y1_to(fixedNodes), Z1_to(fixedNodes), 40, 'bs', '.');
for i = 1:length(loadDofs)
    n = loadNodes(i); d = loadDirs(i); dv = zeros(1,3); dv(d) = sign(loadVals(i)) * arrowScale;
    quiver3(X1_to(n)-dv(1), Y1_to(n)-dv(2), Z1_to(n)-dv(3), dv(1), dv(2), dv(3), 0, 'r', 'LineWidth', 2.5, 'MaxHeadSize', 0.5);
end

cb1 = colorbar; cb1.Label.String = '|U| TopOpt [mm]';
colormap(ax1, 'turbo');
clim(ax1, [0, max(Umag_topopt)+eps]); 
view(ax1, getGoodView(FEM));

% =========================================================================
% SUBPLOT 2: FEBio Displacement
% =========================================================================
ax2 = subplot(1, 3, 2);
hold on; axis equal; grid on; box on;
title('FEBio Tetrahedral Displacement', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('x'); ylabel('y'); zlabel('z');

drawMeshEdgesColored(FEM, X1_fb, Y1_fb, Z1_fb, Umag_febio, chi);

% Replot BCs
scatter3(X1_fb(fixedNodes), Y1_fb(fixedNodes), Z1_fb(fixedNodes), 40, 'bs', '.');
for i = 1:length(loadDofs)
    n = loadNodes(i); d = loadDirs(i); dv = zeros(1,3); dv(d) = sign(loadVals(i)) * arrowScale;
    quiver3(X1_fb(n)-dv(1), Y1_fb(n)-dv(2), Z1_fb(n)-dv(3), dv(1), dv(2), dv(3), 0, 'r', 'LineWidth', 2.5, 'MaxHeadSize', 0.5);
end

cb2 = colorbar; cb2.Label.String = '|U| FEBio [mm]';
colormap(ax2, 'turbo');
clim(ax2, [0, max(Umag_febio)+eps]); 
view(ax2, getGoodView(FEM));

% =========================================================================
% SUBPLOT 3: Absolute Error (Deviation)
% =========================================================================
ax3 = subplot(1, 3, 3);
hold on; axis equal; grid on; box on;
title('Absolute Spatial Deviation', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('x'); ylabel('y'); zlabel('z');

% Render the error colors overlaid on the FEBio deformed mesh structure
drawMeshEdgesColored(FEM, X1_fb, Y1_fb, Z1_fb, Umag_error, chi);

cb3 = colorbar; cb3.Label.String = 'Mismatch Magnitude [mm]';
colormap(ax3, 'hot');
disp(max(Umag_error))
clim(ax3, [0, 2*mean(Umag_error)+eps]); 
view(ax3, getGoodView(FEM));

% Link camera movements across all three plots so they rotate together
Link = linkprop([ax1, ax2, ax3], {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', Link);

end

% =========================================================================
% HELPER: DRAW UNCOLORED EDGES
% =========================================================================
function drawMeshEdges(FEM, X, Y, Z, color, linewidth, chi)
    if nargin >= 7 && ~isempty(chi)
        active_elems = find(chi(:,1) > 0.1);
    else
        active_elems = (1:FEM.nele)';
    end
    num_active = length(active_elems);
    if num_active == 0, return; end

    nodes = (FEM.edofMat(active_elems, 1:3:end) + 2) / 3;
    local_edges = [1 2; 2 3; 3 4; 4 1; 5 6; 6 7; 7 8; 8 5; 1 5; 2 6; 3 7; 4 8];

    X = X(:); Y = Y(:); Z = Z(:);
    num_edges = 12 * num_active;
    
    X_lines = NaN(3, num_edges);
    Y_lines = NaN(3, num_edges);
    Z_lines = NaN(3, num_edges);

    for k = 1:12
        idx_start = (k-1)*num_active + 1;
        idx_end   = k*num_active;

        n1 = nodes(:, local_edges(k,1));
        n2 = nodes(:, local_edges(k,2));

        X_lines(1, idx_start:idx_end) = X(n1)';
        X_lines(2, idx_start:idx_end) = X(n2)';
        Y_lines(1, idx_start:idx_end) = Y(n1)';
        Y_lines(2, idx_start:idx_end) = Y(n2)';
        Z_lines(1, idx_start:idx_end) = Z(n1)';
        Z_lines(2, idx_start:idx_end) = Z(n2)';
    end

    plot3(X_lines(:), Y_lines(:), Z_lines(:), 'Color', color, 'LineWidth', linewidth);
end

% =========================================================================
% HELPER: DRAW COLORED EDGES
% =========================================================================
function drawMeshEdgesColored(FEM, X, Y, Z, Umag, chi)
    if nargin >= 6 && ~isempty(chi)
        active_elems = find(chi(:,1) > 0.1);
    else
        active_elems = (1:FEM.nele)';
    end
    num_active = length(active_elems);
    if num_active == 0, return; end

    nodes = (FEM.edofMat(active_elems, 1:3:end) + 2) / 3;
    local_edges = [1 2; 2 3; 3 4; 4 1; 5 6; 6 7; 7 8; 8 5; 1 5; 2 6; 3 7; 4 8];

    avgU = mean(Umag(nodes), 2); % <--- FIXED: Now correctly references node subset

    num_edges = 12 * num_active;
    X_lines = zeros(2, num_edges);
    Y_lines = zeros(2, num_edges);
    Z_lines = zeros(2, num_edges);
    C_lines = zeros(2, num_edges);

    X = X(:); Y = Y(:); Z = Z(:);

    for k = 1:12
        idx_start = (k-1)*num_active + 1;
        idx_end   = k*num_active;

        n1 = nodes(:, local_edges(k,1));
        n2 = nodes(:, local_edges(k,2));

        X_lines(1, idx_start:idx_end) = X(n1)';
        X_lines(2, idx_start:idx_end) = X(n2)';
        Y_lines(1, idx_start:idx_end) = Y(n1)';
        Y_lines(2, idx_start:idx_end) = Y(n2)';
        Z_lines(1, idx_start:idx_end) = Z(n1)';
        Z_lines(2, idx_start:idx_end) = Z(n2)';

        C_lines(1, idx_start:idx_end) = avgU';
        C_lines(2, idx_start:idx_end) = avgU';
    end

    surface(X_lines, Y_lines, Z_lines, C_lines, ...
            'FaceColor', 'none', ...
            'EdgeColor', 'flat', ...
            'MeshStyle', 'column', ... % Prevents the stray diagonal lines
            'LineWidth', 1.5);
end

function v = getGoodView(FEM)

if FEM.nelz == 1
    % Quasi-2D: view from slightly above the XY plane
    v = [0, 0];   % straight-on 2D view (azimuth=0, elevation=90)
else
    % Full 3D: isometric-style view
    v = [35, 25];
end
end

