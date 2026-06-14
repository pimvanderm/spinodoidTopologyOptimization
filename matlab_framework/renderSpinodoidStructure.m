% =========================================================================
% SPINODOID RENDERING
% RBF Blending of discrete, explicitly defined elements
% =========================================================================




% % 1. RENDERING SETUP
% settings.nelx = params.nelx; settings.nely = params.nely; settings.nelz = params.nelz;
% settings.domainSize = [settings.nelx, settings.nely, settings.nelz]; % 1 unit per element
% 
% numberOfWaves = 2;
% settings.numberOfWaves = numberOfWaves;
% settings.res = 5 * numberOfWaves / settings.nelz; 
% settings.k_wave = numberOfWaves * 1 / settings.nelz * 2 * pi;
% settings.writeSTL = 0;
% 
% settings.problemType = params.problemType;
% settings.volfrac = params.volfrac;
% 
% tr = renderSpinodoidStructure(chi_final,settings);

function tr = renderSpinodoidStructure(chi_final,settings)
tic;
% Unpack settings to variables
nelx = settings.nelx; nely = settings.nely; nelz = settings.nelz;
nEle = nelx*nely*nelz;
domainSize = settings.domainSize;
numberOfWaves = settings.numberOfWaves;
res = settings.res;
resolution = domainSize * res;
k_wave = settings.k_wave;
kappa = settings.kappa; % Adjust to make transition sharper or smoother
if settings.graphs == 1
h = waitbar(0, 'Analyzing Spinodoids...'); % Initialize
end

% Load TopOpt output into specific array structures
rhos = chi_final(:,1);
thetas = chi_final(:,2:4);
eulers = chi_final(:,5:7);

% Generate massive shared wave pool
numGlobalWaves = 1000;
globalWaveDirs = randn(numGlobalWaves, 3);
globalWaveDirs = globalWaveDirs ./ vecnorm(globalWaveDirs, 2, 2);
globalPhases = rand(numGlobalWaves, 1) * 2 * pi;

% Pre-calculate global angles for speed
angles1 = min(acosd(globalWaveDirs(:,1)), acosd(-globalWaveDirs(:,1)));
angles2 = min(acosd(globalWaveDirs(:,2)), acosd(-globalWaveDirs(:,2)));
angles3 = min(acosd(globalWaveDirs(:,3)), acosd(-globalWaveDirs(:,3)));

% =========================================================================
% 2. EXPLICIT ELEMENT DEFINITIONS (Your chi_final stand-in)
% =========================================================================
% We define 4 distinct elements. 
% Centers are at 0.5, 1.5, etc. (Centered in each 1x1x1 bounding box)

% Generate all element centers in one line
[Z, Y, X] = ndgrid(1:nelz, 1:nely, 1:nelx);

% Element centers (subtract 0.5 for center)
cx = X(:) - 0.5; cy = Y(:) - 0.5; cz = Z(:) - 0.5;
centers = [cx, cy, cz];

% =========================================================================
% 3. ELEMENT ASSEMBLY AND RBF BLENDING
% =========================================================================
assembled_GRF = zeros(resolution(1), resolution(2), resolution(3));
weight_sum    = zeros(resolution(1), resolution(2), resolution(3));

fprintf('Starting generation...\n');
start_time = toc;

% Loop over our discrete elements
for e = 1:nEle
    elements(e) = struct('cx', centers(e,1), 'cy', centers(e,2), 'cz', centers(e,3), ...
                     'rho', rhos(e), 'theta', thetas(e,:), 'euler', eulers(e,:));
    
    % Extract parameters for this specific element
    cx = elements(e).cx; cy = elements(e).cy; cz = elements(e).cz;
    rho   = elements(e).rho;
    theta = elements(e).theta;
    euler = elements(e).euler;
    
    % Determine Levelset
    if rho ~= 0
        levelset_local = sqrt(2) * erfinv(2 * rho - 1);
    else 
        levelset_local = -5;
    end
    
    % Calculate local evaluation box (Padding = 1.5 units to allow blending)
    pad = 1.5;
    x_idx = max(1, floor((cx - pad)*res)) : min(resolution(1), ceil((cx + pad)*res));
    y_idx = max(1, floor((cy - pad)*res)) : min(resolution(2), ceil((cy + pad)*res));
    z_idx = max(1, floor((cz - pad)*res)) : min(resolution(3), ceil((cz + pad)*res));
    
    [X_loc, Y_loc, Z_loc] = ndgrid(x_idx/res - 0.5/res, ...
                                   y_idx/res - 0.5/res, ...
                                   z_idx/res - 0.5/res);
    
    % Gaussian RBF Weights based on distance from THIS element's center
    dist_sq = (X_loc - cx).^2 + (Y_loc - cy).^2 + (Z_loc - cz).^2;
    
    W = exp(-kappa * dist_sq);
    
    % Rotate the Global Wave Pool into the local element frame
    a = deg2rad(euler(1)); b = deg2rad(euler(2)); c = deg2rad(euler(3));
    Rx = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)];
    Ry = [cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)];
    Rz = [cos(c) -sin(c) 0; sin(c) cos(c) 0; 0 0 1];
    R_mat = Rz * Ry * Rx;
    
    localWaveDirs = globalWaveDirs * R_mat;
    
    ang1 = min(acosd(localWaveDirs(:,1)), acosd(-localWaveDirs(:,1)));
    ang2 = min(acosd(localWaveDirs(:,2)), acosd(-localWaveDirs(:,2)));
    ang3 = min(acosd(localWaveDirs(:,3)), acosd(-localWaveDirs(:,3)));
    
    % Geometric filter based on THIS element's theta
    allowed_idx = (ang1 < theta(1)) | (ang2 < theta(2)) | (ang3 < theta(3));
              
    % Evaluate GRF using the UNROTATED global vectors
    activeDirs   = globalWaveDirs(allowed_idx, :); 
    activePhases = globalPhases(allowed_idx);
    numActive    = size(activeDirs, 1);
    
    % 6. VECTORIZED GRF EVALUATION
    % Flatten the 3D local coordinate grids into an (N_voxels x 3) matrix
    coords_flat = [X_loc(:), Y_loc(:), Z_loc(:)];
    
    % Matrix Multiplication: (N_voxels x 3) * (3 x W) = (N_voxels x W)
    % This calculates the dot product for every voxel and every wave simultaneously!
    dotProducts = coords_flat * activeDirs';
    
    % Evaluate the cosine field
    % activePhases' is (1 x W). MATLAB's implicit expansion broadcasts it to every voxel row.
    wave_matrix = sqrt(2/numActive) * cos(dotProducts * k_wave + activePhases');
    
    % Sum across the wave dimension (dimension 2) to collapse it into a single scalar field
    element_GRF_flat = sum(wave_matrix, 2);
    
    % Reshape back into the 3D voxel grid format
    element_GRF = reshape(element_GRF_flat, size(X_loc));
    
    % --- THE DENSITY TRICK ---
    element_GRF = element_GRF - levelset_local;
    
    % 7. ACCUMULATE
    assembled_GRF(x_idx, y_idx, z_idx) = assembled_GRF(x_idx, y_idx, z_idx) + (W .* element_GRF);
    weight_sum(x_idx, y_idx, z_idx)    = weight_sum(x_idx, y_idx, z_idx) + W;
    if settings.graphs ==1
    % Update waitbar
    waitbar(e/nEle, h, sprintf('Progress: %d%% (%d/%d)', floor(e/nEle*100), e, nEle));
    end
end
if settings.graphs == 1
delete(h);
end

fprintf('Rendered Gaussian Random Field. Time taken: %.2f seconds.\n', toc-start_time);

% =========================================================================
% 4. GLOBAL NORMALIZATION AND RENDERING
% =========================================================================
weight_sum(weight_sum == 0) = 1; 
assembled_GRF = assembled_GRF ./ weight_sum;

[Xg, Yg, Zg] = ndgrid((1:resolution(1))/res, (1:resolution(2))/res, (1:resolution(3))/res);

[f, v] = isosurface(Xg, Yg, Zg, assembled_GRF, 0);
[fc, vc] = isocaps(Xg, Yg, Zg, assembled_GRF, 0, 'enclose', 'below');

% Scale structure to 80x20x20 mm
v = v * 20 / nelz;
vc = vc * 20 / nelz;

%%
if settings.graphs == 1
figure('Name', '2x2 Discrete Spinodoid Assembly', 'Color', 'w', 'Position', [200, 200, 800, 600]);






patch('Faces', f, 'Vertices', v, 'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'none');
patch('Faces', fc, 'Vertices', vc, 'FaceColor', [0.9 0.5 0.1], 'EdgeColor', 'none');

axis equal;% axis([0 2*domainSize(1) 0 2*domainSize(2) 0 2*domainSize(3)]);
view(15, 20); % Tilted 3D view to see all 4 quadrants clearly
camlight headlight; lighting gouraud; material dull;
grid on; box on;
title('2x2 Explicit Element Assembly', 'FontSize', 16, 'FontWeight', 'bold');
end
%%



% Join sets
[f,v]=joinElementSets({f,fc},{v,vc});
tr = triangulation(f,v);
tStamp = datestr(now, 'yyyy-mm-dd_HHMMSS');
fileName = sprintf('TopOpt_%s_%.0fx%.0fx%.0f_v0%.0f_k%.0f_%s.stl', ...
                    settings.problemType,settings.nelx,settings.nely,settings.nelz,round(10*settings.volfrac),numberOfWaves,tStamp);
if settings.writeSTL == 1
stlwrite(TR,fileName)
end

end