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