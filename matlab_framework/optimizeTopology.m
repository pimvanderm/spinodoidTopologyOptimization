%% TEST MAIN OPTIMIZATION LOOP
% clear; clc; close all;

% %% EXAMPLE FEM SCRIPT
% % Parameters
% params.nelx = 12;
% params.nely = 5;
% params.nelz = 5;
% params.problemType = 'three_point_bending_half';
% 
% % Material and NN
% params.rho_min = 0.3;
% params.rho_max = 0.6;
% params.theta_min = 15;
% params.lambda1 = 50;
% params.lambda2 = 50;
% params.rho_nn_min = 0;
% params.rho_nn_max = 1;
% params.theta_nn_min = 0;
% params.theta_nn_max = 90;
% params.C_scale = 2800;
% params.dlnetF_path = 'dlnetF.mat';
% params.eps = 1e-3;
% 
% % Optimization
% params.maxiter = 1000;
% params.tol = 5e-4;
% params.stagnation_tol = 1e-3;
% params.volfrac = 0.3;
% params.move = 0.15;
% params.beta = 5;
% params.eta = 0.5;
% params.r_filter = 1.5;
% params.alpha_max = 180;
% 
% nEle = params.nelx*params.nely*params.nelz;
% params.nele = nEle;
% 
% % Initial design
% params.rho_init = params.volfrac;
% params.theta_init = [90 90 90];
% params.alpha_init = 0;
% 
% load('dlnetF');
% FEM = setupFEM(params.nelx,params.nely,params.nelz,params.problemType);
% 
% nele = FEM.nele;
% maxiter = params.maxiter;
% tol = params.tol;
% 
% % Initialize design variables
% chi_raw = [params.rho_init .* ones(nele, 1), ...
%            params.theta_init .* ones(nele, 3), ...
%            params.alpha_init * ones(nele, 3)];
% [chi, chi_nn] = applyTransformations(chi_raw, params);
% chi_norm = normalizeForNN(chi_nn, params);
% C_all = predictStiffness(chi_norm, dlnetF, params);
% C_rot_all = rotateStiffness(C_all, chi);
% 
% Ke_all = zeros(24, 24, nele);
% for e = 1:nele
%     Ke_all(:,:,e) = elementStiffness(C_rot_all(:,:,e));
% end
% 
% K = assembleGlobalK(FEM, Ke_all);
% U = solveFEM(K, FEM);
% 
% 
% [c, dC_dKe_scale] = computeCompliance(U, FEM, Ke_all);
% 
% visualizeFEM(FEM, U, chi, 1/max(max(U)))
% %% Topology Optimization Parameters
% % Parameters
% params.nelx = 44;
% params.nely = 20;
% params.nelz = 20;
% params.problemType = 'three_point_bending_half';
% 
% % Material and NN
% params.rho_min = 0.3;
% params.rho_max = 0.6;
% params.theta_min = 15;
% params.lambda1 = 50;
% params.lambda2 = 50;
% params.rho_nn_min = 0;
% params.rho_nn_max = 1;
% params.theta_nn_min = 0;
% params.theta_nn_max = 90;
% params.C_scale = 2800;
% params.dlnetF_path = 'dlnetF.mat';
% params.eps = 1e-3;
% 
% % Optimization
% params.maxiter = 1000;
% params.tol = 5e-4;
% params.stagnation_tol = 5e-3;
% params.volfrac = 0.3;
% params.move = 0.15;
% params.beta = 5;
% params.eta = 0.5;
% params.r_filter = 1.5;
% params.alpha_max = 180;
% 
% nEle = params.nelx*params.nely*params.nelz;
% params.nele = nEle;
% 
% % Initial design
% params.rho_init = params.volfrac;
% params.theta_init = [0 15 15];
% params.alpha_init = 0;
% 
% %%
% % Run optimization
% [chi_raw_history, chi_history,c_history, U_final] = optimizeTopology(params);
% %%
% 
% 
% finalIter = find(all(chi_raw_history(:,:,1) ~= 0, 2), 1, 'last');
% chi_final = squeeze(chi_raw_history(finalIter,:,:));
% chi_final = snapToValidRanges(chi_final,params);
% %%
% if params.problemType == 'three_point_bending_half'
%     [chi_final, params] = mirrorHalfStructure(chi_final,params);
% end
% %%
% topOptResults = struct();
% topOptResults.params = params;
% topOptResults.rawHistory = chi_raw_history;
% topOptResults.chi_final = chi_final;
% topOptResults.cHistory = c_history;
% %%
% tStamp = datestr(now, 'yyyy-mm-dd_HHMM');
% fileName = sprintf('TopOpt_%s_%.0fx%.0fx%.0f_v0%.0f_%s', ...
%                     params.problemType,params.nelx,params.nely,params.nelz,round(10*params.volfrac),tStamp);
% save(fileName,'-struct','topOptResults');
% %%
% % Plot convergence
% figure;
% subplot(2,1,1);
% plot(c_history, 'b-', 'LineWidth', 2);
% xlabel('Iteration');
% ylabel('Compliance');
% title('Convergence History');
% grid on;
% 
% subplot(2,1,2);
% change = abs(diff(c_history)) ./ c_history(1:end-1);
% semilogy(change, 'r-', 'LineWidth', 2);
% xlabel('Iteration');
% ylabel('Relative Change');
% title('Convergence Rate');
% grid on;


function [chi_raw_history,chi_history, compliance_history, U_final] = optimizeTopology(params)
% OPTIMIZETOPOLOGY Main topology optimization loop
%
% INPUTS:
%   params - Structure with all parameters (mesh, material, optimization)
%
% OUTPUTS:
%   chi_raw_history     - (maxiter x nele x 7) design variable history
%   compliance_history  - (maxiter x 1) compliance history

% Initialize
FEM = setupFEM(params.nelx, params.nely, params.nelz, params.problemType);

load(params.dlnetF_path, 'dlnetF');
[H, Hs] = buildFilter(FEM, params);

nele = FEM.nele;
maxiter = params.maxiter;
tol = params.tol;

% Initialize design variables
chi_raw = [params.rho_init .* ones(nele, 1), ...
           params.theta_init .* ones(nele, 3), ...
           params.alpha_init * ones(nele, 3)];


% Preallocate history
chi_raw_history = zeros(maxiter, nele, 7);
chi_history = zeros(maxiter,nele,7);
compliance_history = zeros(maxiter, 1);
window_size = 5;  % Number of iterations to check for stagnation
stagnation_tol = params.stagnation_tol;  % Tolerance for stagnation detection

fprintf('\n=== TOPOLOGY OPTIMIZATION ===\n');
fprintf('Elements: %d, Max iterations: %d\n\n', nele, maxiter);
fprintf('%4s %12s %12s %12s %12s %10s %10s %10s %8s\n', 'Iter', 'Compliance', 'Change', 'Av Change', 'Volume', 'Max dRho', 'Max dTheta', 'Max dAng', 'Time(s)');
fprintf('%s\n', repmat('-', 1, 88)); % Extended the dashed line to fit the new columns

c_old = inf;

% Initialize figure handle (before loop)
if params.graphs == 1
progress_fig = [];
hWaitbar = waitbar(0, 'Iteration 1', 'Name', 'Solving problem','CreateCancelBtn','delete(gcbf)');
end

for iter = 1:maxiter
    
    tic;
    
    % =====================================================================
    % FORWARD PASS
    % =====================================================================
    
    [chi, chi_nn] = applyTransformations(chi_raw, params);
    chi_norm = normalizeForNN(chi_nn, params);
    C_all = predictStiffness(chi_norm, dlnetF, params);
    C_rot_all = rotateStiffness(C_all, chi);
    
    Ke_all = zeros(24, 24, nele);
    for e = 1:nele
        Ke_all(:,:,e) = elementStiffness(C_rot_all(:,:,e));
    end
    
    K = assembleGlobalK(FEM, Ke_all);
    U = solveFEM(K, FEM);
    
    
    [c, dC_dKe_scale] = computeCompliance(U, FEM, Ke_all);
    
    % =====================================================================
    % SENSITIVITY ANALYSIS
    % =====================================================================
    dc_dchi_raw = computeSensitivities(U, FEM, Ke_all, C_all, chi, chi_raw, dlnetF, params);

    % =====================================================================
    % SENSITIVITY FILTERING
    % =====================================================================
    rho_current = chi_raw(:,1);
    
    % 1. Density: Standard Sigmund Filter
    dc_drho = dc_dchi_raw(:,1);
    dc_dchi_raw(:,1) = (H * (rho_current .* dc_drho)) ./ (Hs .* max(1e-3, rho_current));

    % 2. All Angles (Columns 2 through 7): Density-Weighted Filter
    % We weight the angle sensitivities by density so void elements 
    % (which have meaningless gradients) don't pollute the structural load paths!
    for k = 2:7
        dc_dangle = dc_dchi_raw(:,k);
        dc_dchi_raw(:,k) = (H * (rho_current .* dc_dangle)) ./ (Hs .* max(1e-3, rho_current));
    end
    
    % =====================================================================
    % DESIGN UPDATE
    % =====================================================================
    if iter == 1
        params.moveSize = 5*ones(1,3)./max(dc_dchi_raw(:,5:7));
    end
    if iter == 1
        dc_dchi_raw_old = dc_dchi_raw;
    else
        % Blend 50% of the current gradient with 50% of the previous
        % This mathematically destroys 2-period checkerboard limit cycles
        dc_dchi_raw = 0.5 * dc_dchi_raw + 0.5 * dc_dchi_raw_old;
        dc_dchi_raw_old = dc_dchi_raw; % Save for next iteration
    end
    
    if iter == 1
        
        params.grad_scale_theta = max(max(abs(dc_dchi_raw(:,2:4))));
        
    end
    [chi_raw, lmid] = updateDesignVars(chi_raw, dc_dchi_raw, params,FEM);
    [chi,~] = applyTransformations(chi_raw,params);
    
    % =====================================================================
    % CONVERGENCE CHECK
    % =====================================================================
    change = abs(c - c_old) / c;
    vol = mean(chi(:,1));
    
    
    % Store history
    chi_raw_history(iter, :, :) = chi_raw;
    chi_history(iter,:,:) = chi;
    compliance_history(iter) = c;
    % Compute moving average change (to detect stagnation despite oscillations)
    if iter >= 2*window_size
        c_window = compliance_history(iter - window_size + 1 : iter);
        movingAvgChange = (mean(compliance_history(iter - window_size + 1 : iter))-mean(compliance_history(iter - 2*window_size + 1 : iter - window_size)))/mean(compliance_history(iter - window_size + 1 : iter));
    else
        movingAvgChange = inf;  % Not enough data yet
    end

    % =====================================================================
    % CALCULATE MAX CHANGES (Density and Euler Angles)
    % =====================================================================
    if iter == 1
        max_drho = 0.0;
        max_dtheta = 0.0;
        max_dalpha = 0.0;
    else
        % Extract previous state
        prev_chi_raw = squeeze(chi_raw_history(iter-1, :, :));
        
        % 1. Max change in density
        max_drho = max(abs(chi_raw(:,1) - prev_chi_raw(:,1)));
        
        % 2. Max change in Euler angles (Cols 5-7), accounting for wrapping
        max_dtheta = max(max(abs(chi_raw(:,2:4) - prev_chi_raw(:,2:4))));
        
        % 3. Max change in Euler angles (Cols 5-7), accounting for wrapping
        d_angle = abs(chi_raw(:,5:7) - prev_chi_raw(:,5:7));
        
        % NOTE: If you are wrapping to [-180, 180], change the 180 below to 360
        d_angle_wrapped = min(d_angle, 180 - d_angle); 
        
        max_dalpha = max(d_angle_wrapped(:));
    end

    iter_time = toc;
    fprintf('%4d %12.6f %12.6e %12.6e %12.6f %10.5f %10.3f %10.3f %8.2f\n', ...
            iter, c, change, movingAvgChange, vol, max_drho, max_dtheta , max_dalpha, iter_time);

    % Check convergence
    if max_drho < 0.005 && max_dtheta < 0.5 && max_dalpha < 0.5 && iter > 10
        fprintf('\nConverged at iteration %d (change = %.6e)\n', iter, change);
        chi_raw_history = chi_raw_history(1:iter, :, :);
        chi_history = chi_history(1:iter,:,:);
        compliance_history = compliance_history(1:iter);
        break;
    end
    
    % Criterion 2: Stagnation (oscillating within small window)
    if iter >= window_size && abs(movingAvgChange) < stagnation_tol
        fprintf('\nConverged due to oscillation at iteration %d (change = %.6e)\n', iter, change);
        chi_raw_history = chi_raw_history(1:iter, :, :);
        chi_history = chi_history(1:iter,:,:);
        compliance_history = compliance_history(1:iter);
        break;
    end
    
    c_old = c;
    % =====================================================================
    % VISUALIZATION (every N iterations, non-intrusive)
    % =====================================================================
    if params.graphs == 1
    progress_fig = plotOptimizationProgress(iter, chi, FEM, params, progress_fig);
    

    drawnow;
    if ~ishandle(hWaitbar)
        % Stop the if cancel button was pressed
        disp('Stopped by user');
        delete(findall(0,'type','figure','tag','TMWWaitbar'))
        break;
    else
        % Update the wait bar
        waitbar(iter/maxiter,hWaitbar, ['Iteration ' num2str(iter)]);
    end
    end

    
    
end

fprintf('\n=== OPTIMIZATION COMPLETE ===\n');
fprintf('Final compliance: %.6f\n', c);
fprintf('Final volume: %.6f (target: %.6f)\n', vol, params.volfrac);
delete(findall(0,'type','figure','tag','TMWWaitbar'))
% Final visualization
if params.graphs ==1
visualizeFEM(FEM, U, chi, 1/max(max(U)));
title(sprintf('Final Design (Iter %d, C = %.4f)', iter, c));
end

U_final = U;

end


%%
function [FEM] = setupFEM(nelx, nely, nelz, problemType)
% SETUPFEM Initialize finite element mesh, boundary conditions, and loads

% Initialize structure
FEM.nelx = nelx;
FEM.nely = nely;
FEM.nelz = nelz;
FEM.nele = nelx * nely * nelz;

% Total number of nodes and DOFs
nnodx = nelx + 1;
nnody = nely + 1;
nnodz = nelz + 1;
nnode = nnodx * nnody * nnodz;
FEM.ndof = 3 * nnode;

% BUILD ELEMENT CONNECTIVITY (edofMat)
% Hex element: 8 nodes × 3 DOFs = 24 DOFs per element
% Node numbering:
%   Bottom face (z-): 1(0,0),2(1,0),3(1,1),4(0,1)
%   Top face (z+): 5(0,0),6(1,0),7(1,1),8(0,1)

nodenrs = reshape(1:nnode, nnodz, nnody, nnodx);
% For each element (ix, iy, iz), its corner node at (ix, iy, iz) is:
% nodenrs(iz, iy, ix)
% The 8 corners are at offsets (+0 or +1) in each direction

% Build list of "anchor" nodes — the (iz=1,iy=1,ix=1) corner of each element
% Element ordering: z varies fastest, then y, then x
anchorNodes = nodenrs(1:nelz, 1:nely, 1:nelx); % nelz x nely x nelx
anchorNodes = anchorNodes(:);                   % nele x 1, column vector

% Strides in node numbering (how many nodes to step in each direction)
sz = 1;           % stride in z direction
sy = nnodz;       % stride in y direction  
sx = nnodz*nnody; % stride in x direction

% The 8 nodes of each element, expressed as offsets from the anchor node
% Convention (matches standard hex element numbering):
%   Node 1: (iz,   iy,   ix  ) -> offset 0
%   Node 2: (iz,   iy,   ix+1) -> offset sx
%   Node 3: (iz,   iy+1, ix+1) -> offset sy + sx
%   Node 4: (iz,   iy+1, ix  ) -> offset sy
%   Node 5: (iz+1, iy,   ix  ) -> offset sz
%   Node 6: (iz+1, iy,   ix+1) -> offset sz + sx
%   Node 7: (iz+1, iy+1, ix+1) -> offset sz + sy + sx
%   Node 8: (iz+1, iy+1, ix  ) -> offset sz + sy
nodeOffsets = [0, sx, sy+sx, sy, sz, sz+sx, sz+sy+sx, sz+sy];

% Build edofMat: nele x 24
% For each element, list the 3 DOFs of each of its 8 nodes
edofMat = zeros(FEM.nele, 24);
for n = 1:8
    globalNode = anchorNodes + nodeOffsets(n); % nele x 1
    edofMat(:, 3*n-2) = 3*globalNode - 2;     % x-DOF
    edofMat(:, 3*n-1) = 3*globalNode - 1;     % y-DOF
    edofMat(:, 3*n  ) = 3*globalNode;          % z-DOF
end
FEM.edofMat = edofMat;

% Prepare sparse assembly indices
iK = kron(edofMat, ones(24,1))';
jK = kron(edofMat, ones(1,24))';
FEM.iK = iK(:);
FEM.jK = jK(:);

% BOUNDARY CONDITIONS
switch lower(problemType)
    case 'cantilever'
        % Fixed: left face (x=0)
        fixedNodes = find(nodenrs(:, :, 1));
        fixedDofs = [3*fixedNodes(:)-2; 3*fixedNodes(:)-1; 3*fixedNodes(:)];
        
        % Load: upward at bottom-right corner
        loadNode = [(nelx)*(nely+1)*(nelz+1)+1:nelz+1: (nelx+1)*(nely+1)*(nelz+1)];

        F = sparse(3*loadNode, 1, 1, FEM.ndof, 1);
    
    case 'compression'
        % Fixed: left face (x=0)
        fixedNodes = find(nodenrs(:, :, 1));
        fixedDofs = [3*fixedNodes(:)-2; 3*fixedNodes(:)-1; 3*fixedNodes(:)];
        
        % Load: rightward at right face (x = L)
        loadNode = [(nelx)*(nely+1)*(nelz+1)+1:1:(nelx+1)*(nely+1)*(nelz+1)];

        F = sparse(3*loadNode-2, 1, 1, FEM.ndof, 1);
    
    case 'cantilever_halfway'
        % Fixed: left face (x=0)
        fixedNodes = find(nodenrs(:, :, 1));
        fixedDofs = [3*fixedNodes(:)-2; 3*fixedNodes(:)-1; 3*fixedNodes(:)];
        
        % Load: upward halfway along the right face (x = L)
        loadNode = [(nelx)*(nely+1)*(nelz+1)+ceil(nelz/2)+1 (nelx)*(nely+1)*(nelz+1)+(nelz+1)+ceil(nelz/2)+1];
        F = sparse(3*loadNode, 1, 1, FEM.ndof, 1);
    
    case 'cantilever_compression'
        % Fixed: left face (x=0)
        fixedNodes = find(nodenrs(:, :, 1));
        fixedDofs = [3*fixedNodes(:)-2; 3*fixedNodes(:)-1; 3*fixedNodes(:)];
        
        % Load: upward halfway along the right face (x = L)
        loadNode = [(nelx)*(nely+1)*(nelz+1)+ceil(nelz/2)+1 (nelx)*(nely+1)*(nelz+1)+(nelz+1)+ceil(nelz/2)+1];
        F = sparse(3*loadNode-2, 1, -1, FEM.ndof, 1);    

    case 'three_point_bending'
        % Physical dimensions tracking: Find space offset relative to center
        space = nelx/2 - 2*nelz;
        
        % 1. Roller Supports: Bottom left and bottom right lines
        % Only constrain movement in the Z-direction
        supportNodes = [nodenrs(1, :, 1+space), nodenrs(1, :, nnodx-space)];
        dofs_rollers_z = 3 * supportNodes(:); % Z-direction DOFs
        
        % 2. Top Loading Line: Located exactly halfway along the length
        % Constrain this entire line in the X-direction to prevent sliding left/right
        loadNodes = nodenrs(nelz+1, :, nelx/2 + 1);
        loadNodes = loadNodes(:);
        dofs_load_x = 3 * loadNodes - 2; % X-direction DOFs
        
        % 3. Middle Stabilizer Node: Find the center-most node along the loading line
        % Constrain a single node in the Y-direction to stop rigid-body axial telescoping
        mid_y_idx = round(nely/2) + 1;
        centerLoadNode = nodenrs(nelz+1, mid_y_idx, nelx/2 + 1);
        dofs_center_y = 3 * centerLoadNode - 1; % Y-direction DOF
        
        % Combine all unique fixed DOFs cleanly
        fixedDofs = unique([dofs_rollers_z; dofs_load_x; dofs_center_y]);
        
        % 4. Apply Force: Introduce a total downward load distributed over the load line
        totalLoad = -50; % Negative for downward bending force
        F = sparse(3*loadNodes, 1, totalLoad/length(loadNodes), FEM.ndof, 1);

    case 'three_point_bending_half'
        % 1. Symmetry Plane: x = 0 (first index in the x-direction)
        % Constrain all nodes on this face in the X-direction to enforce symmetry
        symNodes = nodenrs(:, :, 1);
        dofs_sym_x = 3 * symNodes(:) - 2; % X-direction DOFs
        
        % 2. Roller Support: located at x = 2*nelz + 1 on the bottom face (z = 1)
        % Constrain the line in the Z-direction
        supportNodes = nodenrs(nelz + 1, :, 1);
        dofs_roller_z = 3 * supportNodes(:); % Z-direction DOFs
        
        % 3. Middle Stabilizer Node: Center of the loading line (x = 1, z = nelz+1)
        % Constrain a single node in the Y-direction to prevent rigid body telescoping
        mid_y_idx = round(nely/2) + 1;
        centerLoadNode = nodenrs(nelz+1, mid_y_idx, 1);
        dofs_center_y = 3 * centerLoadNode - 1; % Y-direction DOF
        
        % Combine all unique fixed DOFs cleanly
        fixedDofs = unique([dofs_sym_x; dofs_roller_z; dofs_center_y]);
        
        % 4. Apply Force: 25 N distributed over the loading line (x = 1, z = nelz+1)
        loadNodes = nodenrs(1, :, 2*nelz + 1);
        loadNodes = loadNodes(:);
        totalLoad = 25; % Negative for downward bending force
        
        F = sparse(3*loadNodes, 1, totalLoad/length(loadNodes), FEM.ndof, 1);

    case 'ge_bracket'
        % 1. Extract Node Coordinates (0-indexed to match element grid)
        [IZ, IY, IX] = ndgrid(1:nnodz, 1:nnody, 1:nnodx);
        Xn = IX(:) - 1;
        Yn = IY(:) - 1;
        Zn = IZ(:) - 1;
        
        % 2. Geometric Parameters (Matching the mask dimensions)
        t_base = 0.20 * nelz;
        base_hole_r = 0.08 * nelx;
        x_holes = [0.15, 0.85] * nelx;
        y_holes = [0.15, 0.85] * nely;
        
        pin_x = 0.80 * nelx;
        pin_z = 0.80 * nelz;
        pin_hole_r = 0.12 * nelz;

        % 3. FIXED DOFS (The 4 base holes)
        fixedNodes = [];
        for i = 1:2
            for j = 1:2
                % Distance squared from hole center in XY plane
                r2 = (Xn - x_holes(i)).^2 + (Yn - y_holes(j)).^2;
                
                % Grab nodes strictly inside or on the boundary of the hole, 
                % up to the thickness of the base plate
                hole_nodes = find((r2 <= base_hole_r^2) & (Zn <= t_base));
                fixedNodes = [fixedNodes; hole_nodes];
            end
        end
        fixedNodes = unique(fixedNodes);
        
        % Fix in all 3 directions (x, y, z)
        fixedDofs = [3*fixedNodes(:)-2; 3*fixedNodes(:)-1; 3*fixedNodes(:)];
        
        % 4. LOADED DOFS (Main pin hole)
        % Distance squared from pin center in XZ plane
        r2_pin = (Xn - pin_x).^2 + (Zn - pin_z).^2;
        
        % Grab nodes making up the pin hole
        loadNodes = find(r2_pin <= pin_hole_r^2);
        
        % Apply load: Let's do a classic downward/shear load combination. 
        % Or purely vertical (-Z). We distribute the total load equally across all nodes.
        totalLoad = -100;
        
        % Applying load in the Z-direction (3rd DOF)
        F = sparse(3*loadNodes, 1, totalLoad/length(loadNodes), FEM.ndof, 1);
        
    otherwise
        error('Unknown problem type: %s', problemType);
end



FEM.fixedDofs = fixedDofs;
FEM.freeDofs = setdiff(1:FEM.ndof, fixedDofs);
FEM.F = F;

fprintf('FEM Setup: %dx%dx%d = %d elements, %d DOFs\n', ...
    nelx, nely, nelz, FEM.nele, FEM.ndof);
end

function [Ke] = elementStiffness(C)
% ELEMENTSTIFFNESS Compute element stiffness matrix for 8-node hexahedral element
%
% INPUTS:
%   C - 6x6 material stiffness matrix (Voigt notation)
%
% OUTPUTS:
%   Ke - 24x24 element stiffness matrix
%
% Uses 2x2x2 Gauss integration

% Gauss points and weights (2-point quadrature)
gp = [-1 1] / sqrt(3);
gw = [1 1];

Ke = zeros(24, 24);

% Loop over Gauss points
for i = 1:2
    for j = 1:2
        for k = 1:2
            xi = gp(i);
            eta = gp(j);
            zeta = gp(k);
            w = gw(i) * gw(j) * gw(k);
            
            % Shape function derivatives in natural coordinates
            dN = 1/8 * [
                -(1-eta)*(1-zeta),  (1-eta)*(1-zeta), (1+eta)*(1-zeta), -(1+eta)*(1-zeta), ...
                -(1-eta)*(1+zeta),  (1-eta)*(1+zeta), (1+eta)*(1+zeta), -(1+eta)*(1+zeta);
                -(1-xi)*(1-zeta), -(1+xi)*(1-zeta),  (1+xi)*(1-zeta),   (1-xi)*(1-zeta), ...
                -(1-xi)*(1+zeta), -(1+xi)*(1+zeta),  (1+xi)*(1+zeta),   (1-xi)*(1+zeta);
                -(1-xi)*(1-eta), -(1+xi)*(1-eta), -(1+xi)*(1+eta), -(1-xi)*(1+eta), ...
                 (1-xi)*(1-eta),  (1+xi)*(1-eta),  (1+xi)*(1+eta),  (1-xi)*(1+eta)
                ];
            
            % Jacobian for unit element (1x1x1)
            J = dN * [0 1 1 0 0 1 1 0; 0 0 1 1 0 0 1 1; 0 0 0 0 1 1 1 1]';
            
            % Derivatives in physical coordinates
            dNdx = J \ dN;
            
            % B-matrix (6 x 24) for strain-displacement
            
            B = zeros(6, 24);
            for n = 1:8
                B(:, 3*n-2:3*n) = [
                    dNdx(1,n),         0,         0;  % 11 (xx)
                            0, dNdx(2,n),         0;  % 22 (yy)
                            0,         0, dNdx(3,n);  % 33 (zz)
                            0, dNdx(3,n), dNdx(2,n);  % 23 (yz) <-- FIXED
                    dNdx(3,n),         0, dNdx(1,n);  % 13 (xz) <-- FIXED
                    dNdx(2,n), dNdx(1,n),         0   % 12 (xy) <-- FIXED
                    ];
            end
            
            % Add to stiffness matrix
            Ke = Ke + w * B' * C * B * det(J);
        end
    end
end

end

function [C] = isotropicStiffness(E, nu)
% ISOTROPICSTIFFNESS Create 6x6 isotropic stiffness matrix
%
% INPUTS:
%   E  - Young's modulus
%   nu - Poisson's ratio
%
% OUTPUTS:
%   C - 6x6 stiffness matrix (Voigt notation)

lambda = E * nu / ((1 + nu) * (1 - 2*nu));
mu = E / (2 * (1 + nu));

C = [
    lambda + 2*mu, lambda, lambda, 0, 0, 0;
    lambda, lambda + 2*mu, lambda, 0, 0, 0;
    lambda, lambda, lambda + 2*mu, 0, 0, 0;
    0, 0, 0, mu, 0, 0;
    0, 0, 0, 0, mu, 0;
    0, 0, 0, 0, 0, mu
    ];

end

function [K] = assembleGlobalK(FEM, Ke_all)
% ASSEMBLEGLOBALK Assemble global stiffness matrix from element matrices
%
% INPUTS:
%   FEM    - FEM structure from setupFEM (needs .iK, .jK, .ndof, .nele)
%   Ke_all - 24x24xnele array of element stiffness matrices
%
% OUTPUTS:
%   K - ndof x ndof sparse global stiffness matrix

% Flatten all element stiffness matrices into a single vector
% Each Ke is 24x24 = 576 entries; total length = 576 * nele
sK = reshape(Ke_all, 24*24, FEM.nele);  % (576 x nele)
sK = sK(:);                              % flatten to column vector

% Assemble into sparse matrix
K = sparse(FEM.iK, FEM.jK, sK, FEM.ndof, FEM.ndof);

% Symmetrize to remove floating point asymmetry
K = (K + K') / 2;

end

function [U] = solveFEM(K, FEM)
% SOLVEFEM Solve the FEM system KU = F with boundary conditions
%
% INPUTS:
%   K   - Global stiffness matrix (sparse, ndof x ndof)
%   FEM - FEM structure with .F, .freeDofs, .fixedDofs, .ndof
%
% OUTPUTS:
%   U - Full displacement vector (ndof x 1)

% Initialize displacement vector
U = zeros(FEM.ndof, 1);

% Extract free DOF submatrices
% We only solve for free DOFs; fixed DOFs remain zero
K_free = K(FEM.freeDofs, FEM.freeDofs);
F_free = FEM.F(FEM.freeDofs);

% Solve using MATLAB's sparse direct solver
% The backslash operator automatically selects an appropriate method
% (typically Cholesky for symmetric positive definite K)
U(FEM.freeDofs) = K_free \ F_free;

% Verify boundary conditions are satisfied
assert(max(abs(U(FEM.fixedDofs))) < 1e-12, ...
    'Fixed DOF displacements are non-zero — check boundary conditions');

% Check for non-finite values (indicates singular K)
if any(~isfinite(U))
    error('Displacement contains Inf or NaN — K may be singular. Check BCs.');
end

end

function [compliance, dC_dKe_scale] = computeCompliance(U, FEM, Ke_all)
% COMPUTECOMPLIANCE Compute compliance and element-wise sensitivities
%
% INPUTS:
%   U      - Displacement vector (ndof x 1)
%   FEM    - FEM structure with .edofMat, .F, .nele
%   Ke_all - 24x24xnele element stiffness matrices
%
% OUTPUTS:
%   compliance     - Scalar total compliance c = F^T * U
%   dC_dKe_scale   - (nele x 1) sensitivity of c w.r.t. a uniform scaling
%                    of each element's stiffness matrix
%                    i.e., dc/d(scale_e) = -Ue' * Ke * Ue
%
% NOTE on sensitivities:
%   The analytical sensitivity of compliance w.r.t. element stiffness is:
%       dc/dKe = -Ue * Ue'   (a 24x24 outer product)
%   For the chain rule through the NN, we need dc/d(design_vars).
%   We return the scalar dc/d(scale_e) = -Ue'*Ke*Ue here, which is the
%   element strain energy. The full chain rule is:
%       dc/d(chi_e) = (dc/dCe) : (dCe/d(chi_e))
%   where dCe/d(chi_e) comes from automatic differentiation of the NN.

% Total compliance: c = F^T * U (equivalent to U^T * K * U for self-adjoint)
compliance = full(FEM.F' * U);

% Element-wise sensitivities
dC_dKe_scale = zeros(FEM.nele, 1);

for e = 1:FEM.nele
    % Extract element displacement vector (24 x 1)
    Ue = U(FEM.edofMat(e, :));
    
    % Element strain energy (= -dc/d(scale_e))
    % This is the sensitivity of compliance w.r.t. scaling element e's stiffness
    dC_dKe_scale(e) = -Ue' * Ke_all(:,:,e) * Ue;
end

% % Sanity check: all sensitivities should be non-positive
% % (scaling up any element can only reduce or maintain compliance)
% if any(dC_dKe_scale > 1e-10)
%     warning('Positive compliance sensitivity detected — check element stiffness matrices');
% end

end

function visualizeFEM(FEM, U, chi, scale)
% VISUALIZEFEM Visualize FEM solution with boundary conditions and loads
% Signature updated to (FEM, U, chi, scale)

if nargin < 4
    scale = 1;
end
if nargin < 3
    chi = [];
end

% -------------------------------------------------------------------------
% EXTRACT NODE COORDINATES
% -------------------------------------------------------------------------
nnodx = FEM.nelx + 1;
nnody = FEM.nely + 1;
nnodz = FEM.nelz + 1;
nnode = nnodx * nnody * nnodz;

nodenrs = reshape(1:nnode, nnodz, nnody, nnodx);
[IZ, IY, IX] = ndgrid(1:nnodz, 1:nnody, 1:nnodx);
X0 = IX(:) - 1;   
Y0 = IY(:) - 1;   
Z0 = IZ(:) - 1;   

% Extract nodal displacements
Ux = U(1:3:end);   
Uy = U(2:3:end);   
Uz = U(3:3:end);   

% Deformed coordinates
X1 = X0 + scale * Ux;
Y1 = Y0 + scale * Uy;
Z1 = Z0 + scale * Uz;
Umag = sqrt(Ux.^2 + Uy.^2 + Uz.^2);

% --- SMART FILTERING (Identify active nodes for scatter3) ---
if ~isempty(chi)
    solid_elems = find(chi(:,1) > 0.1);
    solid_nodes = unique(ceil(FEM.edofMat(solid_elems, :) / 3));
else
    solid_nodes = 1:nnode;
end

% -------------------------------------------------------------------------
% IDENTIFY BOUNDARY CONDITION NODES & LOADS
% -------------------------------------------------------------------------
fixedNodes = unique(ceil(FEM.fixedDofs / 3));
[loadDofs, ~, loadVals] = find(FEM.F);
loadNodes = ceil(loadDofs / 3);
loadDirs  = mod(loadDofs - 1, 3) + 1;   
dirNames  = {'x', 'y', 'z'};
arrowScale = 0.4 * max(FEM.nelx,FEM.nely);

% -------------------------------------------------------------------------
% SETUP FIGURE
% -------------------------------------------------------------------------
figure('Name', 'FEM Solution', 'Color', 'white', 'Position', [100, 100, 1400, 500]);

% =========================================================================
% SUBPLOT 1: Undeformed mesh with boundary conditions
% =========================================================================
ax1 = subplot(1, 2, 1);
hold on; axis equal; grid on; box on;
title('Boundary Conditions', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('x'); ylabel('y'); zlabel('z');

% Draw undeformed mesh edges (Filtered by chi)
drawMeshEdges(FEM, X0, Y0, Z0, [0.7 0.7 0.7], 0.5);

% Highlight fixed nodes
scatter3(X0(fixedNodes), Y0(fixedNodes), Z0(fixedNodes), 60, 'bs', 'filled', 'DisplayName', 'Fixed DOFs');

% Draw load arrows
for i = 1:length(loadDofs)
    n = loadNodes(i);
    d = loadDirs(i);
    dv = zeros(1,3);
    dv(d) = sign(loadVals(i)) * arrowScale;
    if i == 1, visFlag = 'on'; else, visFlag = 'off'; end
    quiver3(X0(n) - dv(1), Y0(n) - dv(2), Z0(n) - dv(3), dv(1), dv(2), dv(3), 0, ...
            'r', 'LineWidth', 2.5, 'MaxHeadSize', 0.5, 'DisplayName', sprintf('Load (%s)', dirNames{d}), 'HandleVisibility', visFlag);
end

view(getGoodView(FEM));
colormap(ax1, 'gray');

% =========================================================================
% SUBPLOT 2: Deformed mesh colored by displacement magnitude
% =========================================================================
ax2 = subplot(1, 2, 2);
hold on; axis equal; grid on; box on;
title(sprintf('Deformed Mesh (scale \\times%g)', scale), 'FontSize', 12, 'FontWeight', 'bold');
xlabel('x'); ylabel('y'); zlabel('z');

% Draw deformed mesh colored by Umag (Filtered by chi)
drawMeshEdgesColored(FEM, X1, Y1, Z1, Umag, chi);

% Replot BCs on deformed view
scatter3(X1(fixedNodes), Y1(fixedNodes), Z1(fixedNodes), 40, 'bs', '.');
for i = 1:length(loadDofs)
    n = loadNodes(i); d = loadDirs(i); dv = zeros(1,3); dv(d) = sign(loadVals(i)) * arrowScale;
    quiver3(X1(n)-dv(1), Y1(n)-dv(2), Z1(n)-dv(3), dv(1), dv(2), dv(3), 0, 'r', 'LineWidth', 2.5, 'MaxHeadSize', 0.5);
end

cb2 = colorbar;
cb2.Label.String = 'Displacement magnitude';
colormap(ax2, 'turbo');
clim([0, max(Umag)+eps]); 
view(getGoodView(FEM));

% % =========================================================================
% % SUBPLOT 3: Per-node displacement components
% % =========================================================================
% ax3 = subplot(1, 3, 3);
% hold on; axis equal; grid on; box on;
% title('Displacement Magnitude Map', 'FontSize', 12, 'FontWeight', 'bold');
% xlabel('x'); ylabel('y'); zlabel('z');
% 
% % Draw undeformed mesh background (Filtered by chi)
% drawMeshEdges(FEM, X0, Y0, Z0, [0.7 0.7 0.7], 0.2, chi);
% 
% % Scatter plot only the solid nodes
% scatter3(X0(solid_nodes), Y0(solid_nodes), Z0(solid_nodes), 80, Umag(solid_nodes), 'filled');
% 
% % Annotate max displacement node (ensure we only look at solid nodes)
% [maxU, maxIdx] = max(Umag(solid_nodes));
% maxNode = solid_nodes(maxIdx);
% text(X0(maxNode)+0.1, Y0(maxNode)+0.1, Z0(maxNode)+0.1, ...
%      sprintf('max=%.3f', maxU), 'FontSize', 9, 'Color', 'red');
% 
% cb3 = colorbar;
% cb3.Label.String = 'Displacement magnitude';
% colormap(ax3, 'turbo');
% clim([0, max(Umag)+eps]); 
% view(getGoodView(FEM));
% 
% sgtitle(sprintf('FEM Analysis | %dx%dx%d mesh | Compliance = %.4f', ...
%         FEM.nelx, FEM.nely, FEM.nelz, full(FEM.F' * U)), ...
%         'FontSize', 13, 'FontWeight', 'bold');
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

function [chi, chi_nn] = applyTransformations(chi_raw, params)
% APPLYTRANSFORMATIONS Apply sigmoid transformations to raw design variables
%
% INPUTS:
%   chi_raw - (nele x 7) raw optimization variables
%             columns: [rho, theta1, theta2, theta3, alpha, beta, gamma]
%             rho in [0,1] nominally; thetas in degrees; angles in degrees
%   params  - Structure with fields:
%               .rho_min   - minimum density (paper: 0.3)
%               .rho_max   - maximum density (paper: 0.7)
%               .theta_min - minimum cone angle in degrees (paper: 15)
%               .lambda1   - sharpness for rho transformation (paper: 4)
%               .lambda2   - sharpness for theta transformation (paper: 4)
%
% OUTPUTS:
%   chi    - (nele x 7) physically transformed variables
%            [rho_t, theta1_t, theta2_t, theta3_t, alpha, beta, gamma]
%            angles remain in degrees, rho dimensionless
%   chi_nn - (nele x 4) variables formatted for NN input BEFORE normalization
%            [rho_t, theta1_t, theta2_t, theta3_t]
%            (Euler angles are not inputs to the NN)

rho_min   = params.rho_min;
rho_max   = params.rho_max;
theta_min = params.theta_min;
lambda1   = params.lambda1;
lambda2   = params.lambda2;

nele = size(chi_raw, 1);
chi  = chi_raw;   % copy — Euler angles pass through unchanged

% --- DENSITY TRANSFORMATION ---
% Maps rho_raw onto (0, rho_max] with a soft lower bound at rho_min
% As rho_raw -> +inf: rho_t -> rho_max (we scale by rho_max/1)
% As rho_raw -> -inf: rho_t -> 0 (sigmoid kills it)
% The sigmoid is centred at rho_min so the transition occurs there
rho_raw = chi_raw(:, 1);
rho_arg = max(rho_raw,rho_min);
rho_t   = rho_arg ./ (1 + exp(-lambda1 * (rho_raw - rho_min/2)));

% Clamp to valid range [rho_min, rho_max] as a safety net
rho_t = min(rho_t, rho_max);
chi(:, 1) = rho_t;

% --- ANGLE TRANSFORMATIONS (theta1, theta2, theta3) ---
% Maps theta_raw onto [theta_min, 90] with soft lower bound
% max(theta_i, theta_min) ensures the argument is always >= theta_min
% The sigmoid is centred at theta_min/2
for i = 2:4
    theta_raw = chi_raw(:, i);
    theta_arg = max(theta_raw, theta_min);
    theta_t   = theta_arg ./ (1 + exp(-lambda2 * (theta_raw - theta_min/2)));
    
    % Clamp to valid range [theta_min, 90] as safety net
    theta_t = min(theta_t, 90);
    chi(:, i) = theta_t;
end

% --- EULER ANGLES (alpha, beta, gamma) ---
% No transformation needed — these are free variables in [-180, 180]
% They are already in degrees and pass through unchanged
% chi(:, 5:7) = chi_raw(:, 5:7);  (already copied above)

% --- OUTPUT FOR NEURAL NETWORK ---
% NN takes [rho, theta1, theta2, theta3] only — no Euler angles
chi_nn = chi(:, 1:4);

end

function chi_norm = normalizeForNN(chi_nn, params)
% NORMALIZEFORNNN Normalize physical design variables for NN input
%
% INPUTS:
%   chi_nn  - (nele x 4) physical variables [rho, theta1, theta2, theta3]
%   params  - Structure with fields:
%               .rho_nn_min   - min rho in training data (0.1)
%               .rho_nn_max   - max rho in training data (0.9)
%               .theta_nn_min - min theta in training data (0, for theta=0 case)
%               .theta_nn_max - max theta in training data (90)
%
% OUTPUTS:
%   chi_norm - (nele x 4) normalized to [0, 1]

chi_norm = chi_nn;

% Normalize density
chi_norm(:,1) = (chi_nn(:,1) - params.rho_nn_min) / ...
                (params.rho_nn_max - params.rho_nn_min);

% Normalize angles (same bounds for theta1, theta2, theta3)
for i = 2:4
    chi_norm(:,i) = (chi_nn(:,i) - params.theta_nn_min) / ...
                    (params.theta_nn_max - params.theta_nn_min);
end

% % Warn if any value is outside [0,1] — indicates input outside training range
% if any(chi_norm(:) < -0.01) || any(chi_norm(:) > 1.01)
%     warning('normalizeForNN: Some inputs fall outside training data range [0,1]');
%     min(chi_norm(:))
% end

% Clamp to [0,1] as safety net
chi_norm = min(max(chi_norm, -0.125), 1);

end

function C_all = predictStiffness(chi_norm, dlnetF, params)
% PREDICTSTIFFNESS Evaluate forward NN to get orthotropic stiffness tensors
%
% INPUTS:
%   chi_norm - (nele x 4) normalized design variables [rho, t1, t2, t3]
%   dlnetF   - Trained dlnetwork (forward NN)
%   params   - Structure with fields:
%               .C_scale - scalar or (9x1) vector to denormalize NN output
%                          (multiply NN output by this to get true stiffness)
%               The 9 outputs correspond to:
%               [C11, C12, C13, C22, C23, C33, C44, C55, C66]
%
% OUTPUTS:
%   C_all - (6 x 6 x nele) orthotropic stiffness matrices in Voigt notation

nele = size(chi_norm, 1);
C_all = zeros(6, 6, nele);

% Convert to dlarray for NN evaluation
% dlnetF expects input as (features x batch) = (4 x nele)
X_dl = dlarray(chi_norm', 'CB');   % 'CB' = Channel x Batch

% Forward pass through NN
% Output: (9 x nele) — the 9 independent components of orthotropic C
Y_dl = predict(dlnetF, X_dl);

% % Extract and denormalize
% Y = extractdata(Y_dl)';   % nele x 9
% 
% % Denormalize NN output to physical stiffness values
% if isscalar(params.C_scale)
%     Y_phys = Y * params.C_scale;
% else
%     Y_phys = Y .* repmat(params.C_scale(:)', nele, 1);
% end

Y_norm = extractdata(Y_dl)';   % nele x 9

% THE FIX: Load min-max parameters and denormalize correctly
persistent normParams
if isempty(normParams)
    load('normalizationParams.mat', 'normParams');
end

out_min = repmat(normParams.output.min, nele, 1);
out_max = repmat(normParams.output.max, nele, 1);

% Restore true physical stiffness ratios!
Y_phys = (Y_norm .* (out_max - out_min) + out_min) * params.C_scale;

% Assemble 6x6 stiffness matrices
% Voigt notation for orthotropic material:
%   [C11  C12  C13   0    0    0  ]
%   [C12  C22  C23   0    0    0  ]
%   [C13  C23  C33   0    0    0  ]
%   [ 0    0    0   C44   0    0  ]
%   [ 0    0    0    0   C55   0  ]
%   [ 0    0    0    0    0   C66 ]
%
% NN output columns: [C11,C12,C13,C22,C23,C33,C44,C55,C66]
%                      1    2    3   4    5   6   7    8    9

for e = 1:nele
    c = Y_phys(e, :);
    C_all(:,:,e) = [
        c(1), c(2), c(3),    0,    0,    0;
        c(2), c(4), c(5),    0,    0,    0;
        c(3), c(5), c(6),    0,    0,    0;
           0,    0,    0, c(7),    0,    0;
           0,    0,    0,    0, c(8),    0;
           0,    0,    0,    0,    0, c(9)
        ];
end

% =========================================================================
% ASSIGN NEAR-ZERO STIFFNESS TO VOID ELEMENTS
% =========================================================================
is_void = chi_norm(:, 1) < 0.01;
C_void = 1e-4 * eye(6);
if sum(is_void) > 0.5
    C_all(:,:,is_void) = repmat(C_void, [1, 1, sum(is_void)]);
end
end

function T = buildRotationMatrix(alpha, beta, gamma)
% BUILDROTATIONMATRIX Build 6x6 Bond transformation matrix for stiffness rotation
%
% INPUTS:
%   alpha - rotation about Z axis (degrees)
%   beta  - rotation about Y axis (degrees)
%   gamma - rotation about X axis (degrees)
%   Convention: R = Rz(alpha) * Ry(beta) * Rx(gamma)  [ZYX / aerospace]
%
% OUTPUTS:
%   T - 6x6 Bond matrix such that C_rotated = T * C * T'
%
% BACKGROUND:
%   For a Voigt stiffness tensor, the transformation under a rotation R is:
%   C_rotated = T(R) * C * T(R)'
%   where T is the Bond matrix built from the 3x3 rotation matrix R.
%   See: Auld (1973), Acoustic Fields in Solids, or Slawinski (2010).

% Convert to radians
a = alpha * pi / 180;   % Z rotation
b = beta  * pi / 180;   % Y rotation
g = gamma * pi / 180;   % X rotation

% Build 3x3 rotation matrix (ZYX convention)
Rz = [cos(a), -sin(a), 0;
      sin(a),  cos(a), 0;
          0,       0,  1];

Ry = [cos(b),  0, sin(b);
          0,   1,      0;
     -sin(b),  0, cos(b)];

Rx = [1,      0,       0;
      0, cos(g), -sin(g);
      0, sin(g),  cos(g)];

R = Rz * Ry * Rx;   % Full rotation matrix (3x3)

% Extract components for readability
l1=R(1,1); l2=R(1,2); l3=R(1,3);
m1=R(2,1); m2=R(2,2); m3=R(2,3);
n1=R(3,1); n2=R(3,2); n3=R(3,3);

% Build Bond matrix T (6x6)
% Voigt ordering: [11, 22, 33, 23, 13, 12]
% NOTE: Some texts use [11,22,33,12,23,13] — we use the ordering that matches
% the stiffness matrix assembled in predictStiffness:
% [C11,C22,C33,C44(=C2323),C55(=C1313),C66(=C1212)]
% This Bond matrix uses Voigt ordering [1,2,3,4,5,6] = [11,22,33,23,13,12]

T = [
    l1^2,      l2^2,      l3^2,         2*l2*l3,       2*l1*l3,       2*l1*l2;
    m1^2,      m2^2,      m3^2,         2*m2*m3,       2*m1*m3,       2*m1*m2;
    n1^2,      n2^2,      n3^2,         2*n2*n3,       2*n1*n3,       2*n1*n2;
    m1*n1,     m2*n2,     m3*n3,  m2*n3+m3*n2,   m1*n3+m3*n1,   m1*n2+m2*n1;
    l1*n1,     l2*n2,     l3*n3,  l2*n3+l3*n2,   l1*n3+l3*n1,   l1*n2+l2*n1;
    l1*m1,     l2*m2,     l3*m3,  l2*m3+l3*m2,   l1*m3+l3*m1,   l1*m2+l2*m1
    ];

end

function C_rot_all = rotateStiffness(C_all, chi)
% ROTATESTIFFNESS Apply Bond rotation to all element stiffness matrices
%
% INPUTS:
%   C_all - (6 x 6 x nele) unrotated stiffness matrices
%   chi   - (nele x 7) design variables; columns 5,6,7 are alpha,beta,gamma
%
% OUTPUTS:
%   C_rot_all - (6 x 6 x nele) rotated stiffness matrices

nele = size(C_all, 3);
C_rot_all = zeros(6, 6, nele);

for e = 1:nele
    alpha = chi(e, 5);
    beta  = chi(e, 6);
    gamma = chi(e, 7);
    
    T = buildRotationMatrix(alpha, beta, gamma);
    C_rot_all(:,:,e) = T * C_all(:,:,e) * T';
end

end

function dc_dchi_raw = computeSensitivities(U, FEM, Ke_all, C_all, ...
                                             chi, chi_raw, dlnetF, params)
% COMPUTESENSITIVITIES Compute compliance sensitivities w.r.t. raw design variables
%
% Full chain rule:
%   dc/dchi_raw = (dc/dKe) * (dKe/dC_rot) * (dC_rot/dC) *
%                (dC/dchi_norm) * (dchi_norm/dchi) * (dchi/dchi_raw)
%
% INPUTS:
%   U        - (ndof x 1) displacement vector
%   FEM      - FEM structure
%   Ke_all   - (24 x 24 x nele) element stiffness matrices
%   C_all    - (6 x 6 x nele) unrotated stiffness matrices from NN
%   chi      - (nele x 7) transformed design variables
%   chi_raw  - (nele x 7) raw design variables
%   dlnetF   - Forward neural network (dlnetwork)
%   params   - Parameter structure
%
% OUTPUTS:
%   dc_dchi_raw - (nele x 7) sensitivities

nele = FEM.nele;
dc_dchi_raw = zeros(nele, 7);

% =========================================================================
% STEP 1: dc/dKe (scalar per element)
% From computeCompliance: dc/d(scale_e) = -Ue'*Ke*Ue = element strain energy
% This is the sensitivity of compliance w.r.t. uniform scaling of Ke
% =========================================================================
[~, dc_dKe_scale] = computeCompliance(U, FEM, Ke_all);
% dc_dKe_scale is (nele x 1), each entry = -Ue'*Ke*Ue

% =========================================================================
% STEP 2: dKe/dC_rot (analytical — elementStiffness is LINEAR in C)
% Ke = integral(B' * C * B) dV
% So dKe_ij / dC_rot_kl = integral(B'_ik * B_lj) dV
% This means: for a perturbation dC, dKe = integral(B' * dC * B) dV
%
% We need dc/dC_rot for each element:
%   dc/dC_rot_kl = sum_ij (dc/dKe_ij) * (dKe_ij/dC_rot_kl)
%                = Ue' * [dKe/dC_rot_kl] * Ue
%
% More compactly: dc/dC_rot = -(Be*Ue) * (Be*Ue)'  [contracted with B]
% which is the 6x6 "strain outer product" matrix
% We compute this by evaluating elementStiffness with each unit basis matrix
% =========================================================================
dc_dC_rot = zeros(6, 6, nele);

% Gauss points (same as elementStiffness)
gp = [-1 1] / sqrt(3);
gw = [1 1];

for e = 1:nele
    Ue = U(FEM.edofMat(e,:));  % 24x1 element displacement
    
    % Accumulate strain-displacement product B*Ue over Gauss points
    eps_e = zeros(6, 1);   % element strain at centroid (approximate)
    
    for i = 1:2
        for j = 1:2
            for k = 1:2
                xi = gp(i); eta = gp(j); zeta = gp(k);
                w  = gw(i) * gw(j) * gw(k);
                
                [B, detJ] = computeB(xi, eta, zeta);
                
                % Strain at this Gauss point
                eps_gp = B * Ue;   % 6x1
                
                % dc/dC_rot += w * eps_gp * eps_gp' * detJ
                % (negative because dc/dKe_scale is negative)
                dc_dC_rot(:,:,e) = dc_dC_rot(:,:,e) + ...
                    w * detJ * (eps_gp * eps_gp');
            end
        end
    end
    
    % Apply sign: dc/dKe_scale = -Ue'*Ke*Ue, so
    % dc/dC_rot_e = -eps_e * eps_e' (integrated)
    dc_dC_rot(:,:,e) = -dc_dC_rot(:,:,e);
end


% =========================================================================
% STEP 3: dC_rot/dC (Bond matrix — analytical)
% C_rot = T * C * T'
% So dC_rot/dC is the linear map: dC_rot = T * dC * T'
% For the chain rule: dc/dC = T' * (dc/dC_rot) * T
% =========================================================================
dc_dC = zeros(6, 6, nele);

for e = 1:nele
    T = buildRotationMatrix(chi(e,5), chi(e,6), chi(e,7));
    dc_dC(:,:,e) = T' * dc_dC_rot(:,:,e) * T;
end

% =========================================================================
% STEP 4: dC/dchi_norm (automatic differentiation through dlnetF)
% We need the Jacobian of the 9 NN outputs w.r.t. the 4 NN inputs
% Shape: (nele x 9 x 4) — for each element, a 9x4 Jacobian
% We use dlfeval + dlgradient, looping over outputs
% =========================================================================
chi_nn    = chi(:, 1:4);
chi_norm  = normalizeForNN(chi_nn, params);

X_dl      = dlarray(chi_norm', 'CB');   % 4 x nele

% Get Jacobian: dY/dX where Y is (9 x nele), X is (4 x nele)
% We compute this column by column (one per NN output)
dY_dX = zeros(9, 4, nele);   % (outputs x inputs x elements)


% Define function that returns scalar sum of one output across batch
% We differentiate w.r.t. X to get gradient for this output
dY_dX = computeNNJacobianRow(dlnetF, chi_norm, nele);

% % CRITICAL: Scale by C_scale to get dC/dchi_norm from d(NN_output)/dchi_norm
% if isscalar(params.C_scale)
%     dY_dX = params.C_scale * dY_dX;
% else
%     % If C_scale is a vector (one per output component)
%     for i = 1:9
%         dY_dX(i, :, :) = params.C_scale(i) * squeeze(dY_dX(i, :, :));
%     end
% end

% THE FIX: Chain rule for Min-Max denormalization
load('normalizationParams.mat', 'normParams');
output_range = normParams.output.max - normParams.output.min;

for i = 1:9
    dY_dX(i, :, :) = output_range(i) * squeeze(dY_dX(i, :, :));
end

dY_dX = params.C_scale * dY_dX;

% STEP 5: Contract dc/dC with dC/dchi_norm
% dc/dchi_norm(e,k) = sum_{i,j} (dc/dC_ij(e)) * (dC_ij/dchi_norm(e,k))
%
% dC_ij/dchi_norm maps NN outputs [C11,C12,C13,C22,C23,C33,C44,C55,C66]
% to the 6x6 C matrix. We need to account for the assembly pattern.
% =========================================================================
dc_dchi_norm = zeros(nele, 4);

dc_dchi_norm = zeros(nele, 4);

% Index map from NN output index to (row,col) in C matrix
% NN outputs: [C11,C12,C13,C22,C23,C33,C44,C55,C66]
%              1    2    3   4    5   6    7    8    9
C_idx = [1,1; 1,2; 1,3; 2,2; 2,3; 3,3; 4,4; 5,5; 6,6];

for e = 1:nele
    for out_idx = 1:9
        r = C_idx(out_idx, 1);
        c = C_idx(out_idx, 2);
        
        for k = 1:4
            % For diagonal entries (r==c): single contribution
            if r == c
                dc_dchi_norm(e,k) = dc_dchi_norm(e,k) + ...
                    dc_dC(r, c, e) * dY_dX(out_idx, k, e);
            else
                % For off-diagonal: both C(r,c) and C(c,r) get the same NN output
                % So we sum contributions from both positions in dc_dC
                dc_dchi_norm(e,k) = dc_dchi_norm(e,k) + ...
                    (dc_dC(r, c, e) + dc_dC(c, r, e)) * dY_dX(out_idx, k, e);
            end
        end
    end
end


% =========================================================================
% STEP 6: dchi_norm/dchi (analytical — normalizeForNN is linear)
% chi_norm_k = (chi_k - min_k) / (max_k - min_k)
% So dchi_norm_k/dchi_k = 1 / (max_k - min_k)
% =========================================================================
scale_rho   = 1 / (params.rho_nn_max   - params.rho_nn_min);
scale_theta = 1 / (params.theta_nn_max - params.theta_nn_min);
dchi_norm_dchi = diag([scale_rho, scale_theta, scale_theta, scale_theta]);
dc_dchi_phys = dc_dchi_norm * dchi_norm_dchi';   % nele x 4


% =========================================================================
% STEP 7: dchi/dchi_raw (analytical — sigmoid transformations)
% Differentiate applyTransformations w.r.t. raw inputs
% =========================================================================
d_transform = computeTransformationJacobian(chi_raw, params);
% d_transform is (nele x 7) — diagonal of the Jacobian (no cross-terms)




% Assemble final sensitivity
% Columns 1-4: rho and theta (go through NN)
for k = 1:4
    dc_dchi_raw(:,k) = dc_dchi_phys(:,k) .* d_transform(:,k);
end
% Columns 5-7: Euler angles (do NOT go through NN, affect rotation only)
% dc/d(euler) = sum_{i,j} (dc/dC_rot_ij) * (dC_rot_ij/d(euler))
% This requires dT/d(euler) which we compute by finite difference
% (T is cheap to evaluate, and this avoids a complex analytical derivation)
euler_fd_eps = 5.0;   % degrees
for k = 5:7
    for e = 1:nele
        chi_plus  = chi(e,:); chi_plus(k)  = chi_plus(k)  + euler_fd_eps;
        chi_minus = chi(e,:); chi_minus(k) = chi_minus(k) - euler_fd_eps;
        
        T_plus  = buildRotationMatrix(chi_plus(5),  chi_plus(6),  chi_plus(7));
        T_minus = buildRotationMatrix(chi_minus(5), chi_minus(6), chi_minus(7));
        
        dC_rot_deuler = (T_plus * C_all(:,:,e) * T_plus' - ...
                         T_minus * C_all(:,:,e) * T_minus') / (2*euler_fd_eps);
        
        % Contract with dc/dC_rot
        dc_dchi_raw(e,k) = sum(sum(dc_dC_rot(:,:,e) .* dC_rot_deuler));
    end
    % Chain through transformation (Euler angles have no transformation, so
    % dchi/dchi_raw = 1 for Euler angles)
end

end


function [B, detJ] = computeB(xi, eta, zeta)

dN = 1/8 * [
    -(1-eta)*(1-zeta),  (1-eta)*(1-zeta),  (1+eta)*(1-zeta), -(1+eta)*(1-zeta), ...
    -(1-eta)*(1+zeta),  (1-eta)*(1+zeta),  (1+eta)*(1+zeta), -(1+eta)*(1+zeta);
    -(1-xi)*(1-zeta),  -(1+xi)*(1-zeta),   (1+xi)*(1-zeta),   (1-xi)*(1-zeta), ...
    -(1-xi)*(1+zeta),  -(1+xi)*(1+zeta),   (1+xi)*(1+zeta),   (1-xi)*(1+zeta);
    -(1-xi)*(1-eta),   -(1+xi)*(1-eta),    -(1+xi)*(1+eta),   -(1-xi)*(1+eta), ...
     (1-xi)*(1-eta),    (1+xi)*(1-eta),     (1+xi)*(1+eta),    (1-xi)*(1+eta)
    ];

nodeCoords = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
J    = dN * nodeCoords;
detJ = det(J);
dNdx = J \ dN;

B = zeros(6, 24);
for n = 1:8
    B(:, 3*n-2:3*n) = [
        dNdx(1,n),         0,         0;  % 11 (xx)
                0, dNdx(2,n),         0;  % 22 (yy)
                0,         0, dNdx(3,n);  % 33 (zz)
                0, dNdx(3,n), dNdx(2,n);  % 23 (yz) <-- Swapped
        dNdx(3,n),         0, dNdx(1,n);  % 13 (xz) <-- Swapped
        dNdx(2,n), dNdx(1,n),         0   % 12 (xy) <-- Swapped
        ];
end

end


function dY_dX = computeNNJacobianRow(dlnetF, chi_norm, nele)
% COMPUTENNJACOBIAN Compute per-element Jacobian of NN outputs w.r.t. inputs
%
% INPUTS:
%   dlnetF    - Forward neural network (dlnetwork)
%   chi_norm  - (nele x 4) normalized design variables
%   nele      - number of elements
%
% OUTPUTS:
%   dY_dX - (9 x 4 x nele) per-element Jacobian
%           dY_dX(i,k,e) = d(C_component_i)/d(chi_norm_k) for element e


% dY_dX = zeros(9, 4, nele);
% 
% for e = 1:nele
%     x_e = dlarray(chi_norm(e,:)', 'CB');
% 
%     J_e = dlfeval(@nnJacobian, dlnetF, x_e);
% 
%     disp(num2str(e))
%     dY_dX(:,:,e) = extractdata(J_e);
% end



    % A small perturbation step for finite differences
    epsilon = 0.05;
    dY_dX = zeros(9, 4, nele);
    
    % 1. Baseline forward pass for the ENTIRE batch
    X0 = dlarray(chi_norm', 'CB');
    Y0 = extractdata(predict(dlnetF, X0)); % Size: [9 x nele]
    
    % 2. Perturb each of the 4 inputs one by one
    for k = 1:4
        chi_pert = chi_norm;
        chi_pert(:, k) = chi_pert(:, k) + epsilon; % Add epsilon to the k-th input
        
        X_pert = dlarray(chi_pert', 'CB');
        Y_pert = extractdata(predict(dlnetF, X_pert)); % Size: [9 x nele]
        
        % Calculate the gradient for this input across all elements
        dY_dX(:, k, :) = (Y_pert - Y0) / epsilon;
    end
end

function J = nnJacobian(net, X)
    Y = forward(net, X);      % evaluate network (traced)
    X = stripdims(X);         % strip format from input
    Y = stripdims(Y);         % strip format from output
    J = dljacobian(Y, X, 1);  % Jacobian: d(Y)/d(X) along dim 1
end

function [y_sum, grad] = forwardAndGrad(net, X, out_idx)

Y = forward(net, X);          % 9 x nele
y_sum = sum(Y(out_idx, :));   % scalar
grad  = dlgradient(y_sum, X);
end


function d_transform = computeTransformationJacobian(chi_raw, params)

nele      = size(chi_raw, 1);
d_transform = ones(nele, 7);   % Euler angles have derivative = 1

rho_min   = params.rho_min;
rho_max   = params.rho_max;
theta_min = params.theta_min;
lambda1   = params.lambda1;
lambda2   = params.lambda2;

% --- Density transformation derivative ---
% rho_t = rho_max * rho / (1 + exp(-lambda1*(rho - rho_min)))
% Let s = sigmoid(lambda1*(rho - rho_min)) = 1/(1+exp(-lambda1*(rho-rho_min)))
% rho_t = rho_max * rho * s
% d(rho_t)/d(rho) = rho_max * (s + rho * s*(1-s)*lambda1)
rho_raw = chi_raw(:,1);
s1 = 1 ./ (1 + exp(-lambda1 * (rho_raw - rho_min/2)));
d_transform(:,1) = s1 + rho_raw .* s1 .* (1-s1) * lambda1;

% --- Angle transformation derivative ---
% theta_t = max(theta,theta_min) / (1 + exp(-lambda2*(theta - theta_min/2)))
% Let s = sigmoid(lambda2*(theta - theta_min/2))
% For theta > theta_min (the active region):
%   theta_t = theta * s
%   d(theta_t)/d(theta) = s + theta * s*(1-s)*lambda2
for i = 2:4
    theta_raw = chi_raw(:,i);
    s2 = 1 ./ (1 + exp(-lambda2 * (theta_raw - theta_min/2)));
    theta_arg = max(theta_raw, theta_min);
    
    % Subgradient: d/dtheta of max(theta,theta_min) = 1 if theta>theta_min, else 0
    d_max = double(theta_raw >= theta_min);
    
    d_transform(:,i) = d_max .* s2 + theta_arg .* s2 .* (1-s2) * lambda2;
end

end

function [H, Hs] = buildFilter(FEM, params)
% BUILDFILTER Precompute filter weight matrix
%
% INPUTS:
%   FEM    - FEM structure
%   params - Structure with .r_filter (filter radius in element lengths)
%
% OUTPUTS:
%   H  - (nele x nele) sparse weight matrix
%   Hs - (nele x 1) row sums of H (for normalization)

nele      = FEM.nele;
r_filter  = params.r_filter;

% Element centroids
% Element ordering: z-fast, y-mid, x-slow (matches setupFEM)
[IZ, IY, IX] = ndgrid(1:FEM.nelz, 1:FEM.nely, 1:FEM.nelx);
ex = IX(:) - 0.5;   % centroid x-coordinate
ey = IY(:) - 0.5;   % centroid y-coordinate
ez = IZ(:) - 0.5;   % centroid z-coordinate

% Build sparse weight matrix
% H(i,j) = max(0, r_filter - dist(centroid_i, centroid_j))
iH = zeros(nele * (2*ceil(r_filter)+1)^3, 1);
jH = zeros(size(iH));
sH = zeros(size(iH));
idx = 0;

for e1 = 1:nele
    for e2 = 1:nele
        dist = sqrt((ex(e1)-ex(e2))^2 + ...
                    (ey(e1)-ey(e2))^2 + ...
                    (ez(e1)-ez(e2))^2);
        if dist <= r_filter
            idx = idx + 1;
            iH(idx) = e1;
            jH(idx) = e2;
            sH(idx) = r_filter - dist;
        end
    end
end

H  = sparse(iH(1:idx), jH(1:idx), sH(1:idx), nele, nele);
Hs = sum(H, 2);

end


function rho_f = applyFilter(rho, H, Hs)
% APPLYFILTER Apply precomputed density filter
%
% INPUTS:
%   rho - (nele x 1) unfiltered densities
%   H   - (nele x nele) filter weight matrix
%   Hs  - (nele x 1) row sums for normalization
%
% OUTPUTS:
%   rho_f - (nele x 1) filtered densities

rho_f = (H * rho) ./ Hs;

end

function [chi_raw_new, lmid] = updateDesignVars(chi_raw, dc_dchi_raw, params,FEM)
% UPDATEDESIGNVARS Optimality Criteria update for all design variables
%
% The OC update for a variable x with sensitivity dc/dx is:
%   x_new = x * (-dc/dx / lambda)^eta
% clamped to [x-move, x+move] and [x_min, x_max]
% where lambda is a Lagrange multiplier enforcing the volume constraint
% (found by bisection on rho only; other variables update freely)
%
% INPUTS:
%   chi_raw     - (nele x 7) current raw design variables
%   dc_dchi_raw - (nele x 7) compliance sensitivities
%   params      - Structure with fields:
%     .volfrac     - target volume fraction
%     .move        - OC move limit (typical: 0.2)
%     .eta         - OC damping (typical: 0.5)
%     .rho_min/max - density bounds
%     .theta_min   - minimum angle (degrees)
%     .alpha_max   - maximum Euler angle magnitude (degrees)
%
% OUTPUTS:
%   chi_raw_new - (nele x 7) updated raw design variables

nele    = FEM.nele;
move    = params.move;
eta     = params.eta;
volfrac = params.volfrac;
theta_min = params.theta_min;
chi_raw_new = chi_raw;

% =========================================================================
% UPDATE RHO (column 1) WITH VOLUME CONSTRAINT VIA BISECTION
% =========================================================================
rho     = chi_raw(:,1);
dc_drho = dc_dchi_raw(:,1);

% Bisection to find Lagrange multiplier lambda
l1 = 0; l2 = 1e9;
while (l2 - l1) / (l1 + l2 + eps) > 1e-6
    lmid = (l1 + l2) / 2;
    
    % OC update
    B_e      = max(eps, -dc_drho / lmid);
    rho_cand = rho .* B_e.^eta;
    
    % Apply move limits and bounds
    rho_new = max(0, ...
              max(rho - move, ...
              min(params.rho_max, ...
              min(rho + move, rho_cand))));
    
    % Check volume constraint
    % --- NEW CODE: Calculate actual physical volume ---
    rho_arg  = max(rho_new, params.rho_min);
    rho_phys = rho_arg ./ (1 + exp(-params.lambda1 * (rho_new - params.rho_min/2)));
    rho_phys = min(rho_phys, params.rho_max); % Clamp to upper bound
    
    % Check volume constraint against physical reality, not raw variables
    if mean(rho_phys) > volfrac
        l1 = lmid;
    else
        l2 = lmid;
    end
end

chi_raw_new(:,1) = rho_new;

% =========================================================================
% UPDATE THETA1, THETA2, THETA3 (columns 2-4) — unconstrained OC
% No volume constraint on angles; they update freely within bounds
% =========================================================================

% % First, compute unconstrained OC updates for all angles
% theta_unconstrained = zeros(nele, 3);
% thetaMove = 0.05;
% for k = 2:4
%     theta     = chi_raw(:,k);
%     dc_dtheta = dc_dchi_raw(:,k);
% 
%     B_e       = max(eps, -dc_dtheta / (mean(abs(dc_dtheta)) + eps));
%     theta_new = theta .* B_e.^eta;
% 
%     % Apply move limits [0, 90] (temporarily allow anything)
%     theta_new = max(0, max(theta - thetaMove * 90, ...
%                 min(90, min(theta + thetaMove * 90, theta_new))));
% 
%     theta_unconstrained(:, k-1) = theta_new;
% end
% 
% % Now enforce constraint: at most one angle < theta_min per element
% for e = 1:nele
%     theta_e = theta_unconstrained(e, :);  % [theta1, theta2, theta3]
%     is_small = theta_e < theta_min;
%     num_small = sum(is_small);
% 
%     % New Case: only intervene if ALL three want to be small
%     if num_small == 3
%         % Get sensitivities for the three angles
%         dc_angles = dc_dchi_raw(e, 2:4);
% 
%         % We want to find the LEAST influential angle.
%         % In this context, "influence" is the magnitude of the sensitivity.
%         % Since these angles are trying to go to 0, the sensitivities are 
%         % likely positive (meaning decreasing the angle decreases compliance).
%         % We find the index of the minimum sensitivity:
%         [~, idx_least_influential] = min(dc_angles);
% 
%         % Set only that one to theta_min, leaving the other two small
%         theta_unconstrained(e, idx_least_influential) = theta_min;
%     end
% 
% end
% 
% % Write back to chi_raw_new
% chi_raw_new(:, 2:4) = theta_unconstrained;

% =========================================================================
% UPDATE THETA1, THETA2, THETA3 — Gradient Descent (NOT OC!)
% =========================================================================
max_theta_step = 15; % Your desired move limit (0.05 * 90)
theta_cand = zeros(nele, 3);
is_solid = chi_raw(:,1) > 0.1;
for k = 2:4
    theta     = chi_raw(:,k);
    dc_dtheta = dc_dchi_raw(:,k);
    
    % 1. Scale-invariant gradient
    dim_grad = dc_dtheta / (params.grad_scale_theta + eps);
    
    % 2. Additive Gradient Descent Step
    % (Tune this learning rate if it moves too slow/fast. 15.0 is a good start)
    learning_rate_theta = 50.0; 
    raw_step = learning_rate_theta * dim_grad;
    
    % 3. Hard Clip to your 4.5 degree move limit
    clipped_step = max(-max_theta_step, min(max_theta_step, raw_step));
    
    % 4. Apply step ONLY to solid elements (starves cancerous growths)
    theta_new = theta - (clipped_step .* is_solid);
    
    % 5. Strictly apply absolute bounds
    theta_cand(:, k-1) = max(0, min(90, theta_new));
end

% 6. Enforce the Spinodoid Constraint smoothly
for e = 1:nele
    if all(theta_cand(e, :) < theta_min)
        valid_idx = find(chi_raw(e, 2:4) >= theta_min - 1e-5);
        if isempty(valid_idx)
            [~, hold_idx] = max(theta_cand(e, :));
        elseif length(valid_idx) == 1
            hold_idx = valid_idx(1);
        else
            dc_angles = dc_dchi_raw(e, 2:4);
            [~, temp_idx] = min(dc_angles(valid_idx));
            hold_idx = valid_idx(temp_idx);
        end
        
        % Because we used additive steps, holding this at theta_min 
        % guarantees we never violate the 4.5 deg limit!
        theta_cand(e, hold_idx) = theta_min;
    end
end

% Write back to chi_raw_new
chi_raw_new(:, 2:4) = theta_cand;


% =========================================================================
% UPDATE EULER ANGLES — Strain-Energy Anchored Descent
% =========================================================================
max_angle_step = 15; 

% 1. Extract physical scale reference (Mean Density Sensitivity)
% This grounds the math in the actual physical strain energy of the mesh,
% preventing noise amplification even if density updates are disabled.
dc_drho = dc_dchi_raw(:,1);
phys_scale = mean(abs(dc_drho));

% Safety catch for zero-load cases
if phys_scale < 1e-12
    phys_scale = 1e-12;
end

% Define a noise floor (adjust if it freezes too early or spins too long)
noise_threshold = 1e-4; 

for k = 5:7 
    alpha     = chi_raw(:,k);
    dc_dalpha = dc_dchi_raw(:,k);
    
    % Calculate the dimensionless, scale-invariant gradient
    dim_grad = dc_dalpha / phys_scale;
    
    % --- THE FIX: Cartwheel Prevention ---
    % If the gradient ratio is smaller than the noise threshold, it's a ghost.
    % Force the gradient to exactly zero.
    dim_grad(abs(dim_grad) < noise_threshold) = 0.0;
    
    % Apply the learning rate
    learning_rate = 150.0;
    raw_step = learning_rate * dim_grad;
    
    % Hard Clip to prevent explosions
    clipped_step = max(-max_angle_step, min(max_angle_step, raw_step));
    
    % Apply additive step and safely wrap to [-90, 90]
    alpha_new = alpha - clipped_step;
    alpha_new = mod(alpha_new + 90, 180) - 90;
    
    chi_raw_new(:,k) = alpha_new;
end

% for k = 5:7 % (Change to 5:7 to optimize all axes later)
%     alpha     = chi_raw(:,k);
%     dc_dalpha = dc_dchi_raw(:,k);
% 
%     % 2. Calculate the dimensionless, scale-invariant gradient
%     dim_grad = dc_dalpha / phys_scale;
% 
%     % 3. Apply the learning rate
%     % (A factor of ~500 translates the dimensionless ratio into a crisp degree step)
%     learning_rate = 500.0;
%     raw_step = learning_rate * dim_grad;
% 
%     % 4. Hard Clip to prevent explosions
%     clipped_step = max(-max_angle_step, min(max_angle_step, raw_step));
% 
%     % 5. Apply additive step and safely wrap to [-90, 90]
%     alpha_new = alpha - clipped_step;
%     alpha_new = mod(alpha_new + 90, 180) - 90;
% 
%     chi_raw_new(:,k) = alpha_new;
% end

% for k = 6
% 
%     alpha     = chi_raw(:,k);
%     dc_dalpha = dc_dchi_raw(:,k);
%     max(dc_dalpha);
%     % Normalized gradient step
%     alpha_new  = alpha - params.moveSize(k-4)*dc_dalpha;
% 
%     % Wrap to [-180, 180]
%     alpha_new = mod(alpha_new + 90, 180) - 90;
% 
%     chi_raw_new(:,k) = alpha_new;
% end

end

function [rho_proj, dH_drho] = applyHeaviside(rho_f,beta)

% APPLYHEAVISIDE Dilated Heaviside projection for minimum solid length scale
%
% Projects filtered density toward 1 for any rho_f > 0, while preserving
% exact void (rho_f = 0 maps to 0). This enforces a minimum feature size
% on solid regions equal to the filter radius, without constraining voids.
%
% Uses the dilated projection from Wang, Lazarov & Sigmund (2011):
%   H_d(rho_f) = 1 - exp(-beta * rho_f) + rho_f * exp(-beta)
%
% INPUTS:
%   rho_f - (nele x 1) filtered density field
%   beta  - sharpness parameter (start ~1, increase to ~32 via continuation)
%
% OUTPUTS:
%   rho_proj - (nele x 1) projected density
%   dH_drho  - (nele x 1) pointwise derivative d(rho_proj)/d(rho_f)

% Smoothed Heaviside centred at eta
% tanh form is numerically stable for large beta
eta = 0.3;
rho_proj = (tanh(beta * eta) + tanh(beta * (rho_f - eta))) ...
         / (tanh(beta * eta) + tanh(beta * (1 - eta)));

% Derivative
dH_drho  = beta * (1 - tanh(beta * (rho_f - eta)).^2) ...
         / (tanh(beta * eta) + tanh(beta * (1 - eta)));

end



function chi_final = snapToValidRanges(chi, params)
% SNAPTOVALIDRANGES Snap near-boundary values to exact bounds
%
% Handles numerical noise that leaves values just outside valid ranges
% (e.g., rho = 0.29999 → 0.3, theta = 14.9999 → 15)
%
% INPUTS:
%   chi    - (nele x 7) design variables [rho, theta1, theta2, theta3, alpha, beta, gamma]
%   params - Structure with .rho_min, .theta_min
%
% OUTPUTS:
%   chi_clean - (nele x 7) cleaned design variables

chi_final = chi;

snap_tol = params.rho_min/2;  % Tolerance for snapping (adjust if needed)

% -------------------------------------------------------------------------
% DENSITY: Snap to [0, rho_min, rho_max]
% -------------------------------------------------------------------------
rho = chi(:, 1);

% Snap near-zero to exactly zero (void)
is_near_zero = (rho > 0) & (rho < snap_tol);
rho(is_near_zero) = 0;

% Snap near-rho_min to exactly rho_min
is_near_min = (rho > params.rho_min - snap_tol) & (rho < params.rho_min);
rho(is_near_min) = params.rho_min;

% Snap near-rho_max to exactly rho_max
is_near_max = (rho > params.rho_max) & (rho < params.rho_max + snap_tol);
rho(is_near_max) = params.rho_max;

chi_final(:, 1) = rho;

% -------------------------------------------------------------------------
% ANGLES (theta1, theta2, theta3): Snap to [0, theta_min, 90]
% -------------------------------------------------------------------------

snap_tol = params.theta_min/2;

for i = 2:4
    theta = chi(:, i);
    
    % Snap near-zero to exactly zero (columnar)
    is_near_zero = (theta > 0) & (theta < snap_tol);
    theta(is_near_zero) = 0;
    
    % Snap near-theta_min to exactly theta_min
    is_near_min = (theta > params.theta_min - snap_tol) & (theta < params.theta_min);
    theta(is_near_min) = params.theta_min;
    
    % Snap near-90 to exactly 90
    is_near_max = theta > 90;
    theta(is_near_max) = 90;
    
    chi_final(:, i) = theta;
end

% -------------------------------------------------------------------------
% EULER ANGLES: No snapping needed (continuous [-180, 180])
% -------------------------------------------------------------------------
% Optionally wrap to [-180, 180] for consistency
for i = 5:7
    chi_final(:, i) = mod(chi(:, i) + 180, 360) - 180;
end

end

function fig = plotOptimizationProgress(iter, chi, FEM, params, fig)
% PLOTOPTIMIZATIONPROGRESS Fast 3D Voxel-based Density Field Visualizer
% Renders solid element cubes colored by density. Voids are omitted completely.
%
% INPUTS:
%   iter   - Current iteration number
%   chi    - (nele x 7) transformed variables OR single-column density vector
%   FEM    - FEM structure with .nelx, .nely, .nelz
%   params - Parameter structure with mesh metrics
%   fig    - Figure handle for non-intrusive viewport updates

    % Extract structural grid dimensions
    nelx = params.nelx;
    nely = params.nely;
    nelz = params.nelz;

    % =========================================================================
    % CREATE OR REUSE VIEWPORT
    % =========================================================================
    if isempty(fig) || ~isvalid(fig)
        fig = figure('Position', [100, 100, 800, 600], ...
                     'Name', 'Voxel Optimization Progress', ...
                     'NumberTitle', 'off', ...
                     'Color', 'white');
    else
        set(0, 'CurrentFigure', fig); % Safe update without forcing window focus
        clf(fig);
    end

    % =========================================================================
    % FILTER STRUTTER VS. VOID MATTER
    % =========================================================================
    if size(chi, 2) >= 7
        rho = chi(:, 1); 
    else
        rho = chi(:); 
    end
    
    % Only construct geometric shapes for active matter to maximize framerates
    void_threshold = 0.1; 
    solid_idx = find(rho > void_threshold);
    N = length(solid_idx);

    if N == 0
        % Text notice if the design domain is completely empty early on
        annotation('textarrow', [0.5 0.5], [0.5 0.5], 'String', 'Initializing Field...', ...
            'HeadStyle', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        drawnow;
        return;
    end

    % Unravel 1D indices to spatial grid matrix index arrays
    [iz, iy, ix] = ind2sub([nelz, nely, nelx], solid_idx);

    % =========================================================================
    % COMPILE GLOBAL PATCH OBJECT (Vectorized Matrix Allocation)
    % =========================================================================
    % Base template definition for a single 1x1x1 solid voxel element
    v_base = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
    f_base = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];

    % Expand coordinate arrays flatly across the entire batch
    V = zeros(N * 8, 3);
    F = zeros(N * 6, 4);

    for i = 1:N
        offset = [ix(i)-1, iy(i)-1, iz(i)-1];
        V((i-1)*8+1 : i*8, :) = v_base + offset;
        F((i-1)*6+1 : i*6, :) = f_base + (i-1)*8;
    end

    % =========================================================================
    % RENDER PIPELINE & IMAGE COLOR STYLING
    % =========================================================================
    % Duplicate element density scalar across all 6 faces of each respective cube
    face_densities = repelem(rho(solid_idx), 6);

    hold on;
    grid on;
    box on;
    axis equal tight;
    view(35, 25); % Standard isometric viewport view

    % Generate the high-performance patch draw call
    patch('Vertices', V, 'Faces', F, ...
          'FaceVertexCData', face_densities, ...
          'FaceColor', 'flat', 'EdgeColor', 'none');

    % Set limits matching your true domain boundaries
    xlim([0, nelx]);
    ylim([0, nely]);
    zlim([0, nelz]);

    xlabel('X'); ylabel('Y'); zlabel('Z');
    set(gca, 'ZDir', 'normal');

    % Apply colormap styling matching your workspace layout
    colormap(jet); 
    clim([0, 1]); 
    
    cb = colorbar;
    cb.Label.String = 'Relative Density \rho';
    cb.Label.FontWeight = 'bold';

    % Clean flat illumination to highlight voxel features without calculation drag
    camlight headlight;
    lighting flat;

    title(sprintf('Iteration %d | Active Voxels: %d | Volume: %.3f', ...
          iter, N, mean(rho)), 'FontSize', 12, 'FontWeight', 'bold');

    % Flash graphics pipeline immediately without blocking optimization processes
    drawnow limitrate;
end

%%

