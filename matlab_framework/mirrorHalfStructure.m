function [chi_final, params] = mirrorHalfStructure(chi_final,params)
% --- Mirroring Script for Symmetric Half-Model ---
% Assumes chi_final (nele x 7) and params are in your workspace

% 1. Extract dimensions
nelx = params.nelx;
nely = params.nely;
nelz = params.nelz;
num_vars = size(chi_final, 2); % Should be 7

% 2. Preallocate the full mirrored array (Double the elements along X)
% New size will be (2 * nele) x 7
chi_full = zeros(2 * params.nele, num_vars);

% 3. Loop through each design variable column to preserve the layout
for v = 1:num_vars
    % Reshape column vector to 3D grid (Z-fast, Y-mid, X-slow)
    var_3D = reshape(chi_final(:, v), nelz, nely, nelx);
    
    % Mirror along the X-axis (Dimension 3)
    % flip(..., 3) reverses the elements along the X-axis
    var_mirrored = flip(var_3D, 3);
    if v == 6 || v == 7
        var_mirrored = -var_mirrored; 
    end
    
    % Stack them together along the X-axis: [Mirrored, Original]
    % This places the mirror plane (x = 0) perfectly at the seam interface
    var_full_3D = cat(3, var_mirrored, var_3D);
    
    % Flatten back into a 1D column vector using column-major ordering (:)
    chi_full(:, v) = var_full_3D(:);
end
chi_final = chi_full;
% 4. Update parameter structure to reflect the full structure size
params.nelx = 2 * nelx;
params.nele = 2 * params.nele;

fprintf('Structure successfully mirrored along x=0.\n');
fprintf('New grid size: %d x %d x %d (%d elements total).\n', ...
        params.nelx, params.nely, params.nelz, params.nele);

end