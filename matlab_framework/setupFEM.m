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