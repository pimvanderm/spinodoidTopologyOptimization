

%% Optional: Run your interactive viewer or plotting functions on the full model
interactive3DTopOptViewer(chi_final, params);
% plotMacroscopicSpinodoidSummary(chi_full, params_full);




%%
nelx = params.nelx;
nely = params.nely;
nelz = params.nelz;
rho_final = squeeze(reshape(chi_final(:,1),nelz,nely,nelx));
% 1. Create the grid and a counter
grid = zeros(2*nelx+1, 2*nelz+1);
count = zeros(2*nelx+1, 2*nelz+1);

for i = 1:nelx
    for j = 1:nelz
        % We identify the 4 corners of the current element
        idx_i = [2*i-1, 2*i+1];
        idx_j = [2*j-1, 2*j+1];
        
        % Add the value to the grid
        % (Removed the 0.25, we will normalize later)
        grid(idx_i, idx_j) = grid(idx_i, idx_j) + rho_final(j,i);
        
        % Keep track of how many elements hit these vertices
        count(idx_i, idx_j) = count(idx_i, idx_j) + 1;
    end
end

% 2. Normalize: Divide by the number of contributors
% This automatically handles corners (1), edges (2), and interior (4)
grid = grid ./ count;

% 3. Extract the vertex grid as you did before
grid = grid(1:2:end, 1:2:end);
figure(5)
[X,Y] = ndgrid(0:nelx,0:nelz);%

h = pcolor(X,Y,grid);
shading interp

% Create a colormap where the first 10 colors are white
k=25;
custom_map = turbo(256-k); 
white_padding = ones(k, 3); 
final_map = [white_padding; custom_map];
colormap(final_map);
clim([0.29,0.7])
%colormap(turbo)
%%
clear v
% --- Setup Video Writer ---
v = VideoWriter('TopOpt_volfrac_20x10', 'MPEG-4');
v.FrameRate = 10; % Frames per second
open(v);

% --- Prepare Figure ---
figure()
figure('Color', 'w');
colormap(flipud(gray)); % 1 (solid) = black, 0 (void) = white


% Assume nelx, nely, nelz are defined in your workspace
% If nely = 1 (2D case), we squeeze it for plotting
num_iterations = size(chi_history, 1);
rho_column = chi_raw_history(:,1);
for it = 1:num_iterations
    % 1. Extract the first parameter (density) for all elements at this iteration
    rho_column = chi_history(it, :, 1);
    
    % 2. Reshape according to your provided logic
    % Result is (nelz x nely x nelx)
    rho_3D = reshape(rho_column, nelz, nely, nelx);
    
    % 3. Create the 2D slice for plotting (Squeezing out the Y-dimension)
    rho_2D = squeeze(rho_3D);
    
    % 4. Render the image
    % We use imagesc for the classic SIMP "pixel" look
    imagesc(rho_2D);
    set(gca, 'YDir', 'normal'); % Ensure z-axis points up
    clim([0 1]); % Force the scale to stay consistent 0 to 1
    axis equal tight;
    
    % Add info tags
    title(['Iteration: ', num2str(it)]);
    xlabel('x'); ylabel('z');
    drawnow; % Force MATLAB to render the frame
    
    % 5. Write to video
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% --- Finalize ---
close(v);
fprintf('Video saved as Topology_Optimization.mp4\n');


%%
fprintf('\n=== GENERATING VISUALIZATION ===\n');
plotOptimizedDesign(chi_final, params);

function plotOptimizedDesign(chi_final, params)
% PLOTOPTIMIZEDDESIGN Create publication-quality visualization of TopOpt results
% Similar to Figure from Zheng et al. 2021 paper
%
% INPUTS:
%   chi_final - (nele x 7) final design variables [rho, theta1-3, alpha, beta, gamma]
%   params    - Structure with .nelx, .nely, .nelz

nelx = params.nelx;
nely = params.nely;
nelz = params.nelz;
nele = nelx * nely * nelz;

% Verify input size
assert(size(chi_final, 1) == nele, ...
    'chi_final has %d elements, expected %d', size(chi_final, 1), nele);

chi_phys = chi_final;
% Reshape design variables into 3D grids for visualization
% Remember: element ordering is z-fast, y-mid, x-slow
rho     = reshape(chi_phys(:,1), nelz, nely, nelx);
theta1  = reshape(chi_phys(:,2), nelz, nely, nelx);
theta2  = reshape(chi_phys(:,3), nelz, nely, nelx);
theta3  = reshape(chi_phys(:,4), nelz, nely, nelx);
beta    = reshape(chi_phys(:,6), nelz, nely, nelx);  % Rotation around y-axis (k=6)

% Create void mask (threshold for "empty")
void_threshold = 0.01;
is_void = rho < void_threshold;

% Create figure
fig = figure('Position', [100, 100, 1400, 900], 'Color', 'white');

% =========================================================================
% SUBPLOT 1: Density
% =========================================================================
subplot(2, 3, 1);
plotSliceWithMask(rho, is_void, nelx, nely, nelz, 'Optimized Relative Density ρ');
colormap(gca, 'jet');
cbar = colorbar;
cbar.Label.String = 'ρ';
cbar.Label.FontSize = 12;
clim([0.3, 0.7]);

% =========================================================================
% SUBPLOT 2: Theta1 (Cone angle 1)
% =========================================================================
subplot(2, 3, 2);
plotSliceWithMask(theta1, is_void, nelx, nely, nelz, 'Optimized Spinodoid Parameter θ₁');
colormap(gca, 'turbo');
cbar = colorbar;
cbar.Label.String = 'θ₁ [°]';
cbar.Label.FontSize = 12;
clim([0, 90]);

% =========================================================================
% SUBPLOT 3: Theta2 (Cone angle 2)
% =========================================================================
subplot(2, 3, 3);
plotSliceWithMask(theta2, is_void, nelx, nely, nelz, 'Optimized Spinodoid Parameter θ₂');
colormap(gca, 'turbo');
cbar = colorbar;
cbar.Label.String = 'θ₂ [°]';
cbar.Label.FontSize = 12;
clim([0, 90]);

% =========================================================================
% SUBPLOT 4: Theta3 (Cone angle 3)
% =========================================================================
subplot(2, 3, 4);
plotSliceWithMask(theta3, is_void, nelx, nely, nelz, 'Optimized Spinodoid Parameter θ₃');
colormap(gca, 'turbo');
cbar = colorbar;
cbar.Label.String = 'θ₃ [°]';
cbar.Label.FontSize = 12;
clim([0, 90]);

% =========================================================================
% SUBPLOT 5: Beta (Rotation around y-axis)
% =========================================================================
subplot(2, 3, 5);
plotSliceWithMask(beta, is_void, nelx, nely, nelz, 'Rotation Angle β (y-axis)');
colormap(gca, 'hsv');  % Circular colormap for angles
cbar = colorbar;
cbar.Label.String = 'β [°]';
cbar.Label.FontSize = 12;
clim([-90, 90]);

% =========================================================================
% SUBPLOT 6: Summary Statistics
% =========================================================================
% subplot(2, 3, 6);
% axis off;
% 
% % Compute statistics
% rho_vec = chi_phys(:,1);
% theta1_vec = chi_phys(:,2);
% theta2_vec = chi_phys(:,3);
% theta3_vec = chi_phys(:,4);
% 
% n_void = sum(rho_vec < void_threshold);
% n_solid = sum(rho_vec > params.rho_min);
% n_columnar = sum((theta1_vec < 1) | (theta2_vec < 1) | (theta3_vec < 1));
% 
% stats_text = sprintf([...
%     '\\bf{Optimization Summary}\n\n' ...
%     '\\rm{Mesh: %d × %d × %d elements}\n' ...
%     'Total elements: %d\n\n' ...
%     '\\bf{Density Distribution:}\n' ...
%     '\\rm{Void (ρ < %.2f): %d (%.1f%%)}\n' ...
%     'Solid (ρ > %.2f): %d (%.1f%%)}\n' ...
%     'Volume fraction: %.3f\n\n' ...
%     '\\bf{Microstructure:}\n' ...
%     '\\rm{Columnar elements: %d (%.1f%%)}\n' ...
%     'Mean θ₁: %.1f°\n' ...
%     'Mean θ₂: %.1f°\n' ...
%     'Mean θ₃: %.1f°\n\n' ...
%     '\\bf{Rotation:}\n' ...
%     '\\rm{Mean β: %.1f°}\n' ...
%     'Std β: %.1f°}\n'], ...
%     nelx, nely, nelz, nele, ...
%     void_threshold, n_void, 100*n_void/nele, ...
%     params.rho_min, n_solid, 100*n_solid/nele, ...
%     mean(rho_vec), ...
%     n_columnar, 100*n_columnar/nele, ...
%     mean(theta1_vec(theta1_vec > 0)), ...
%     mean(theta2_vec(theta2_vec > 0)), ...
%     mean(theta3_vec(theta3_vec > 0)), ...
%     mean(beta(:)), std(beta(:)));
% 
% text(0.1, 0.5, stats_text, 'FontSize', 10, ...
%     'VerticalAlignment', 'middle', 'Interpreter', 'tex');

% Overall title
sgtitle('Spinodoid Topology Optimization Results', ...
    'FontSize', 16, 'FontWeight', 'bold');

end


function plotSliceWithMask(data3D, void_mask, nelx, nely, nelz, titleStr)
% PLOTSLICEWITHMASK Plot 3D field with void regions masked (not rendered)
%
% For quasi-2D in x-z plane (nely = 1): shows x-z slice
% For full 3D: shows x-z slice at mid-y

% Debug: check input dimensions
expected_size = [nelz, nely, nelx];
if ~isequal(size(data3D), expected_size)
    error('plotSliceWithMask: data3D has size [%s], expected [%s]', ...
        mat2str(size(data3D)), mat2str(expected_size));
end

if nely == 1
    % data2D is (nelz x nelx)
    data2D = squeeze(data3D);      
    mask2D = squeeze(void_mask);   
    
    data_masked = data2D;
    data_masked(mask2D) = NaN;
    X = 0.5:nelx+0.5;
    Z = 0.5:nelz+0.5;
    % imagesc plots row indices on y-axis and col indices on x-axis
    % It maps the matrix [z, x] correctly to the screen
    h = imagesc(X,Z,data_masked);
    
    % This handles the NaN regions properly
    % 'AlphaData' controls transparency: 0 for NaN, 1 for valid data
    set(h, 'AlphaData', ~isnan(data_masked));
    
    % Set background color for the 'holes'
    set(gca, 'Color', [0.95 0.95 0.95]);
    
    % Standard orientation: TopOpt usually has z=1 at the bottom
    % imagesc defaults to 'ij' (y-axis pointing down). 
    % We flip it so z increases upwards.
    set(gca, 'YDir', 'normal'); 
    axis equal tight;
    xlabel('x');
    ylabel('z');
else
    % Full 3D: show x-z slice at mid-y
    mid_y = ceil(nely / 2);
    data_slice = squeeze(data3D(:, mid_y, :));      % (nelz x nelx)
    mask_slice = squeeze(void_mask(:, mid_y, :));   % (nelz x nelx)
    
    % Mask void regions
    data_masked = data_slice;
    data_masked(mask_slice) = NaN;
    
    [X, Z] = meshgrid(1:nelx, 1:nelz);
    h = pcolor(X, Z, data_masked);
    shading interp;
    
    set(gca, 'Color', [0.95 0.95 0.95]);
    
    axis equal tight;
    xlabel('x', 'FontSize', 11);
    ylabel('z', 'FontSize', 11);
    title([titleStr, sprintf(' (y = %d/%d)', mid_y, nely)], ...
        'FontSize', 12, 'FontWeight', 'bold');
    
    return;  % Skip the 2D title below
end

title(titleStr, 'FontSize', 12, 'FontWeight', 'bold');

end
%%
function interactive3DTopOptViewer(chi_final, params)
% INTERACTIVE3DTOPOPTVIEWER 3D Interactive Voxel Explorer for Spinodoids

    % 1. Extract Dimensions and Filter Voids
    nelx = params.nelx; nely = params.nely; nelz = params.nelz;
    void_threshold = 0.01;
    
    solid_idx = find(chi_final(:,1) > void_threshold);
    N = length(solid_idx);
    
    if N == 0
        disp('No solid elements found to visualize.');
        return;
    end

    [iz, iy, ix] = ind2sub([nelz, nely, nelx], solid_idx);

    fprintf('Building 3D geometry for %d solid elements...\n', N);
    
    v_base = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
    f_base = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];

    V = zeros(N*8, 3);
    F = zeros(N*6, 4);

    for i = 1:N
        offset = [ix(i)-1, iy(i)-1, iz(i)-1];
        V((i-1)*8+1 : i*8, :) = v_base + offset;
        F((i-1)*6+1 : i*6, :) = f_base + (i-1)*8;
    end

    fig = figure('Name', '3D Spinodoid Explorer', 'Color', 'white', ...
                 'Position', [100, 100, 1000, 800], 'NumberTitle', 'off');
    
    ax = axes('Parent', fig, 'Position', [0.05, 0.2, 0.8, 0.75]);
    
    initial_data = chi_final(solid_idx, 1);
    initial_cdata = repelem(initial_data(:), 6); 
    
    p = patch(ax, 'Vertices', V, 'Faces', F, ...
              'FaceVertexCData', initial_cdata, ... 
              'FaceColor', 'flat', 'EdgeColor', 'none');
          
    axis(ax, 'equal', 'tight', 'off');
    view(ax, 3);
    
    camlight(ax, 'headlight');
    lighting(ax, 'flat');
    material(ax, 'dull');
    
    cb = colorbar(ax);
    cb.FontSize = 11;
    cb.FontWeight = 'bold';

    uicontrol('Style', 'text', 'String', 'Display Variable:', ...
        'Units', 'normalized', 'Position', [0.05 0.08 0.15 0.04], ...
        'BackgroundColor', 'w', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    
    var_dropdown = uicontrol('Style', 'popupmenu', ...
        'String', {'1. Relative Density (ρ)', '2. Cone Angle (θ₁)', '3. Cone Angle (θ₂)', ...
                   '4. Cone Angle (θ₃)', '5. Euler Alpha (α)', '6. Euler Beta (β)', '7. Euler Gamma (γ)'}, ...
        'Units', 'normalized', 'Position', [0.22 0.08 0.25 0.04], ...
        'FontSize', 11, 'Callback', @updatePlot);

    uicontrol('Style', 'text', 'String', 'Front Plane:', ...
        'Units', 'normalized', 'Position', [0.55 0.08 0.1 0.04], ...
        'BackgroundColor', 'w', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    
    % --- THE ULTIMATE SLIDER FIX ---
    % Native 0-to-1 continuous slider. No internal step glitches. 
    if nely > 1
        s_enable = 'on';
    else
        s_enable = 'off';
    end
    
    y_slider = uicontrol('Style', 'slider', 'Min', 0, 'Max', 1, ...
        'Value', 0, 'Units', 'normalized', 'Position', [0.67 0.08 0.25 0.04], ...
        'Enable', s_enable);

    % Add continuous listener for live 3D slicing while dragging!
    addlistener(y_slider, 'ContinuousValueChange', @updatePlot);
    y_slider.Callback = @updatePlot; % Catch standard clicks as well

    updatePlot();

    function updatePlot(~, ~)
        val_idx = var_dropdown.Value;
        
        % Map the 0-1 slider value to the physical y-layer [1 to nely]
        y_front = 1 + round(y_slider.Value * (nely - 1));
        
        current_data = chi_final(solid_idx, val_idx);
        
        % Mask out elements in front of the slider plane
        visible_mask = iy >= y_front;
        current_data(~visible_mask) = NaN;
        
        p.FaceVertexCData = repelem(current_data(:), 6); 
        
        if val_idx == 1
            colormap(ax, 'jet');
            clim(ax, [params.rho_min, params.rho_max]);
            title(ax, sprintf('Relative Density \\rho (Showing Y-Layers %d to %d)', y_front, nely), 'FontSize', 16);
        elseif val_idx >= 2 && val_idx <= 4
            colormap(ax, 'turbo');
            clim(ax, [0, 90]);
            title(ax, sprintf('Spinodoid Cone Angle \\theta_%d (Showing Y-Layers %d to %d)', val_idx-1, y_front, nely), 'FontSize', 16);
        else
            colormap(ax, 'hsv');
            clim(ax, [-90, 90]);
            angle_names = {'\alpha (Z-axis)', '\beta (Y-axis)', '\gamma (X-axis)'};
            title(ax, sprintf('Euler Rotation %s (Showing Y-Layers %d to %d)', angle_names{val_idx-4}, y_front, nely), 'FontSize', 16);
        end
        
        delete(findall(ax, 'Type', 'light'));
        camlight(ax, 'headlight');
        
        % Force graphics pipeline to flush so the UI doesn't lag
        drawnow limitrate;
    end
end


interactive3DTopOptViewer(chi_final, params)

%%
function plotMacroscopicSpinodoidSummary(chi_final, params)
% PLOTMACROSCOPICSPINODOIDSUMMARY Fast glyph-based topology visualizer
% Classifies elements as Lamellar, Columnar, or Cubic and plots rotated proxies.

    fprintf('\n=== GENERATING MACROSCOPIC SUMMARY ===\n');
    
    nelx = params.nelx; nely = params.nely; nelz = params.nelz;
    void_threshold = 0.1; % Only plot structural elements
    
    % Find solid elements
    solid_idx = find(chi_final(:,1) > void_threshold);
    N = length(solid_idx);
    if N == 0, disp('No solid elements to plot.'); return; end
    
    [iz, iy, ix] = ind2sub([nelz, nely, nelx], solid_idx);
    
    % Preallocate massive arrays for the patch object
    % (Over-allocating to 40 faces per element to handle complex proxies)
    max_faces = N * 40; 
    V_all = zeros(max_faces * 4, 3); % 4 vertices per quad face
    F_all = zeros(max_faces, 4);
    C_all = zeros(max_faces, 3);     % RGB colors
    
    v_count = 0; f_count = 0;
    
    % Base Colors (Matching your uploaded image)
    col_lamellar = [0.85, 0.55, 0.55]; % Pinkish
    col_columnar = [0.85, 0.65, 0.15]; % Gold/Yellow
    col_cubic    = [0.35, 0.75, 0.75]; % Teal/Cyan
    col_iso      = [0.60, 0.60, 0.60]; % Grey
    
    for i = 1:N
        e = solid_idx(i);
        rho = chi_final(e, 1);
        thetas = chi_final(e, 2:4);
        angles = chi_final(e, 5:7);
        
        % 1. CLASSIFY THE METAMATERIAL
        % How many angles are effectively "zero" (restricted wave vectors)?
        zero_tol = 5.0; 
        is_zero = thetas < zero_tol;
        num_zeros = sum(is_zero);
        
        % Generate the local proxy geometry (centered at 0,0,0, size 1x1x1)
        if num_zeros >= 2
            % LAMELLAR (1D variation)
            [V_local, F_local] = createLamellarProxy(is_zero, rho);
            elem_color = col_lamellar;
            
        elseif num_zeros == 1
            % COLUMNAR (2D variation)
            [V_local, F_local] = createColumnarProxy(is_zero, rho);
            elem_color = col_columnar;
            
        elseif num_zeros == 0 && all(thetas < 80)
            % CUBIC (3D variation)
            [V_local, F_local] = createCubicProxy(rho);
            elem_color = col_cubic;
            
        else
            % ISOTROPIC / RANDOM (Solid block)
            [V_local, F_local] = createIsotropicProxy(rho);
            elem_color = col_iso;
        end
        
        % 2. ROTATE TO GLOBAL FRAME
        % Use the exact Z-Y-X Euler rotation from your TopOpt code
        a = deg2rad(angles(1)); b = deg2rad(angles(2)); g = deg2rad(angles(3));
        Rz = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
        Ry = [cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)];
        Rx = [1 0 0; 0 cos(g) -sin(g); 0 sin(g) cos(g)];
        R_mat = Rz * Ry * Rx;
        
        V_rot = (R_mat * V_local')';
        
        % 3. TRANSLATE TO ELEMENT POSITION
        offset = [ix(i)-1, iy(i)-1, iz(i)-1];
        V_final = V_rot + offset;
        
        % 4. APPEND TO GLOBAL MESH
        num_v = size(V_final, 1);
        num_f = size(F_local, 1);
        
        V_all(v_count+1 : v_count+num_v, :) = V_final;
        F_all(f_count+1 : f_count+num_f, :) = F_local + v_count;
        C_all(f_count+1 : f_count+num_f, :) = repmat(elem_color, num_f, 1);
        
        v_count = v_count + num_v;
        f_count = f_count + num_f;
    end
    
    % Trim preallocated arrays
    V_all = V_all(1:v_count, :);
    F_all = F_all(1:f_count, :);
    C_all = C_all(1:f_count, :);
    
    % 5. RENDER THE SUMMARY
    fig = figure('Name', 'Macroscopic Spinodoid Summary', 'Color', 'w', 'Position', [100, 100, 1200, 800]);
    ax = axes('Parent', fig);
    
    patch(ax, 'Vertices', V_all, 'Faces', F_all, ...
          'FaceVertexCData', C_all, 'FaceColor', 'flat', 'EdgeColor', 'none');
      
    axis(ax, 'equal', 'tight', 'off');
    view(ax, 3);
    
    % Professional lighting setup
    camlight(ax, 'headlight');
    camlight(ax, 'left');
    lighting(ax, 'flat');
    material(ax, 'dull');
    
    title(ax, 'Macroscopic Spinodoid Orientations', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Add a custom legend
    hold on;
    plot3(NaN, NaN, NaN, 's', 'MarkerSize', 15, 'MarkerFaceColor', col_lamellar, 'MarkerEdgeColor', 'k', 'DisplayName', 'Lamellar');
    plot3(NaN, NaN, NaN, 's', 'MarkerSize', 15, 'MarkerFaceColor', col_columnar, 'MarkerEdgeColor', 'k', 'DisplayName', 'Columnar');
    plot3(NaN, NaN, NaN, 's', 'MarkerSize', 15, 'MarkerFaceColor', col_cubic,    'MarkerEdgeColor', 'k', 'DisplayName', 'Cubic');
    plot3(NaN, NaN, NaN, 's', 'MarkerSize', 15, 'MarkerFaceColor', col_iso,      'MarkerEdgeColor', 'k', 'DisplayName', 'Isotropic');
    legend('Location', 'northeast', 'FontSize', 12);
    
    fprintf('Summary generated successfully.\n');
end

% =========================================================================
% GEOMETRY GENERATORS (Creates 1x1x1 normalized proxy shapes)
% =========================================================================

function [V, F] = createLamellarProxy(is_zero, rho)
    % Creates 3 stacked plates
    thickness = max(0.1, rho * 0.8) / 3; % Scale thickness by density
    gap = (1 - 3*thickness) / 2;
    
    % Determine stack axis (the axis of variation = the non-zero theta)
    stack_axis = find(~is_zero, 1); 
    if isempty(stack_axis), stack_axis = 3; end % Default to Z
    
    V = []; F = [];
    z_start = -0.5;
    for i = 1:3
        [v_box, f_box] = makeBox([-0.5, -0.5, z_start], [0.5, 0.5, z_start + thickness]);
        F = [F; f_box + size(V,1)]; V = [V; v_box];
        z_start = z_start + thickness + gap;
    end
    V = alignLocalAxis(V, stack_axis);
end

function [V, F] = createColumnarProxy(is_zero, rho)
    % Creates 4 pillars in the corners
    w = max(0.1, sqrt(rho) * 0.4); % Width of pillar
    
    % Determine pillar axis (the axis of NO variation = the zero theta)
    col_axis = find(is_zero, 1);
    if isempty(col_axis), col_axis = 3; end
    
    V = []; F = [];
    centers = [-0.25, -0.25; 0.25, -0.25; -0.25, 0.25; 0.25, 0.25];
    for i = 1:4
        cx = centers(i,1); cy = centers(i,2);
        [v_box, f_box] = makeBox([cx-w/2, cy-w/2, -0.5], [cx+w/2, cy+w/2, 0.5]);
        F = [F; f_box + size(V,1)]; V = [V; v_box];
    end
    V = alignLocalAxis(V, col_axis);
end

function [V, F] = createCubicProxy(rho)
    % Creates a 3D cross
    w = max(0.1, rho * 0.6); 
    V = []; F = [];
    % X-bar
    [v_b, f_b] = makeBox([-0.5, -w/2, -w/2], [0.5, w/2, w/2]);
    F = [F; f_b + size(V,1)]; V = [V; v_b];
    % Y-bar
    [v_b, f_b] = makeBox([-w/2, -0.5, -w/2], [w/2, 0.5, w/2]);
    F = [F; f_b + size(V,1)]; V = [V; v_b];
    % Z-bar
    [v_b, f_b] = makeBox([-w/2, -w/2, -0.5], [w/2, w/2, 0.5]);
    F = [F; f_b + size(V,1)]; V = [V; v_b];
end

function [V, F] = createIsotropicProxy(rho)
    % Creates a single scaled cube
    w = max(0.2, rho);
    [V, F] = makeBox([-w/2, -w/2, -w/2], [w/2, w/2, w/2]);
end

% --- Helper: Create a standard 3D Box ---
function [V, F] = makeBox(pmin, pmax)
    x1=pmin(1); y1=pmin(2); z1=pmin(3);
    x2=pmax(1); y2=pmax(2); z2=pmax(3);
    V = [x1 y1 z1; x2 y1 z1; x2 y2 z1; x1 y2 z1; 
         x1 y1 z2; x2 y1 z2; x2 y2 z2; x1 y2 z2];
    F = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
end

% --- Helper: Swap local axes so Z aligns with the target axis ---
function V = alignLocalAxis(V, target_axis)
    if target_axis == 1 % Align Z proxy geometry to X
        V = [V(:,3), V(:,2), -V(:,1)];
    elseif target_axis == 2 % Align Z proxy geometry to Y
        V = [V(:,1), V(:,3), -V(:,2)];
    end
    % If target_axis == 3, it is already aligned to Z
end

plotMacroscopicSpinodoidSummary(chi_final, params)