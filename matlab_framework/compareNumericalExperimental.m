%% Comparison of Numerical Neural Network Predictions vs Experimental Data
% Analyzes the effect of wavenumber, density, and cone angles on prediction accuracy
% Based on Kumar et al. spinodoid metamaterials

clear; clc; close all;

%% 1. LOAD EXPERIMENTAL DATA

fprintf('=== LOADING EXPERIMENTAL DATA ===\n');

% Read Excel file
expData = readtable('ExperimentData.xlsx', 'Sheet', 'Sheet2');
youngModulus = 1600;
% The first row contains headers, extract them
headers = expData{1,:};
data = expData(1:end,:);

% Convert to numeric and create properly structured table
expResults = table();
expResults.rho = data{1:end,1};
expResults.wavenumber_over_pi = data{1:end,2};
expResults.theta1 = data{1:end,3};
expResults.theta2 = data{1:end,4};
expResults.theta3 = data{1:end,5};
expResults.E1_exp = data{1:end,6};
expResults.E2_exp = data{1:end,7};
expResults.E3_exp = data{1:end,8};

%Normalize for the Young's Modulus of FormLabs Clear Resin V4
expResults.E1_exp = expResults.E1_exp / 1600;
expResults.E2_exp = expResults.E2_exp / 1600;
expResults.E3_exp = expResults.E3_exp / 1600;

fprintf('Loaded %d experimental data points\n', height(expResults));
fprintf('Unique densities: %s\n', mat2str(unique(expResults.rho)'));
fprintf('Unique wavenumbers/π: %s\n', mat2str(unique(expResults.wavenumber_over_pi)'));
fprintf('Unique θ₁ values: %s\n', mat2str(unique(expResults.theta1)'));
fprintf('\n');

%% 2. LOAD TRAINED NEURAL NETWORKS

fprintf('=== LOADING TRAINED NETWORKS ===\n');

% Load the trained networks (adjust paths as needed)
if exist('fNet.mat', 'file')
    load('fNet.mat', 'fNet');
    fprintf('✓ Forward network loaded\n');
else
    error('fNet.mat not found. Train the network first or adjust the path.');
end

if exist('normalizationParams.mat', 'file')
    load('normalizationParams.mat', 'normParams');
    fprintf('✓ Normalization parameters loaded\n');
else
    warning('normalizationParams.mat not found. Assuming no normalization.');
    normParams.input.min = [0 0 0 0];
    normParams.input.max = [1 90 90 90];
end

fprintf('\n');

%% 3. GENERATE NUMERICAL PREDICTIONS

fprintf('=== GENERATING NUMERICAL PREDICTIONS ===\n');

% Prepare input for neural network
inputParams = [expResults.rho, expResults.theta1, expResults.theta2, expResults.theta3];

% Normalize if normalization parameters exist
if ~isempty(normParams)
    inputMin = normParams.input.min;
    inputMax = normParams.input.max;
    inputParams_norm = (inputParams - inputMin) ./ (inputMax - inputMin);
    fprintf('Input parameters normalized\n');
else
    inputParams_norm = inputParams;
    fprintf('No normalization applied\n');
end

% Predict using forward network
stiffnessPred_norm = predict(fNet, inputParams_norm);

% Denormalize predictions if needed
if ~isempty(normParams) && isfield(normParams, 'output')
    outputMin = normParams.output.min;
    outputMax = normParams.output.max;
    stiffnessPred = stiffnessPred_norm .* (outputMax - outputMin) + outputMin;
    fprintf('Output denormalized\n');
else
    stiffnessPred = stiffnessPred_norm;
end

% --- Existing Setup ---
% Note: Kumar et al. uses wavenumber β = 10π/l for l×l×l elements
% The stiffness components are:
% 1:C1111, 2:C1122, 3:C1133, 4:C2222, 5:C2233, 6:C3333, 7:C2323, 8:C3131, 9:C1212

num_predictions = size(stiffnessPred, 1);

% Preallocate arrays for the true engineering constants
E1_num = zeros(num_predictions, 1);
E2_num = zeros(num_predictions, 1);
E3_num = zeros(num_predictions, 1);

% Loop through each prediction to handle the matrix inversion
for i = 1:num_predictions
    % 1. Extract the components for the current prediction
    c11 = stiffnessPred(i, 1);
    c12 = stiffnessPred(i, 2);
    c13 = stiffnessPred(i, 3);
    c22 = stiffnessPred(i, 4);
    c23 = stiffnessPred(i, 5);
    c33 = stiffnessPred(i, 6);
    c44 = stiffnessPred(i, 7); % C2323 (G23)
    c55 = stiffnessPred(i, 8); % C3131 (G31)
    c66 = stiffnessPred(i, 9); % C1212 (G12)
    
    % 2. Construct the 6x6 Orthotropic Stiffness Matrix (C)
    C = [c11, c12, c13,   0,   0,   0;
         c12, c22, c23,   0,   0,   0;
         c13, c23, c33,   0,   0,   0;
           0,   0,   0, c44,   0,   0;
           0,   0,   0,   0, c55,   0;
           0,   0,   0,   0,   0, c66];
       
    % 3. Invert to get the Compliance Matrix (S)
    S = inv(C);
    
    % 4. Extract true Young's Moduli (E = 1 / S_ii)
    E1_num(i) = 1 / S(1,1);
    E2_num(i) = 1 / S(2,2);
    E3_num(i) = 1 / S(3,3);
end

% Store the exact predictions (replacing the old approximations)
expResults.E1_num = E1_num;
expResults.E2_num = E2_num;
expResults.E3_num = E3_num;

fprintf('Generated %d numerical predictions\n', height(expResults));
fprintf('\n');

%% 4. CALCULATE ERRORS

fprintf('=== CALCULATING ERRORS ===\n');

% Absolute errors
expResults.error_E1_abs = expResults.E1_exp - expResults.E1_num;
expResults.error_E2_abs = expResults.E2_exp - expResults.E2_num;
expResults.error_E3_abs = expResults.E3_exp - expResults.E3_num;

% Relative errors (percentage)
expResults.error_E1_rel = 100 * (expResults.E1_exp - expResults.E1_num) ./ expResults.E1_exp;
expResults.error_E2_rel = 100 * (expResults.E2_exp - expResults.E2_num) ./ expResults.E2_exp;
expResults.error_E3_rel = 100 * (expResults.E3_exp - expResults.E3_num) ./ expResults.E3_exp;

% Overall metrics
fprintf('Overall Error Statistics:\n');
fprintf('E1 - Mean Abs Error: %.2f, Mean Rel Error: %.2f%%\n', ...
    mean(abs(expResults.error_E1_abs)), mean(abs(expResults.error_E1_rel)));
fprintf('E2 - Mean Abs Error: %.2f, Mean Rel Error: %.2f%%\n', ...
    mean(abs(expResults.error_E2_abs)), mean(abs(expResults.error_E2_rel)));
fprintf('E3 - Mean Abs Error: %.2f, Mean Rel Error: %.2f%%\n', ...
    mean(abs(expResults.error_E3_abs)), mean(abs(expResults.error_E3_rel)));

% R² calculations
R2_E1 = 1 - sum(expResults.error_E1_abs.^2) / sum((expResults.E1_exp - mean(expResults.E1_exp)).^2);
R2_E2 = 1 - sum(expResults.error_E2_abs.^2) / sum((expResults.E2_exp - mean(expResults.E2_exp)).^2);
R2_E3 = 1 - sum(expResults.error_E3_abs.^2) / sum((expResults.E3_exp - mean(expResults.E3_exp)).^2);

fprintf('\nR² Scores:\n');
fprintf('E1: R² = %.4f\n', R2_E1);
fprintf('E2: R² = %.4f\n', R2_E2);
fprintf('E3: R² = %.4f\n', R2_E3);
fprintf('\n');

%% 5. VISUALIZATION: OVERALL COMPARISON, COLOURED ATTRIBUTES

figure('Position', [50 50 1600 500], 'Name', "Effective Young's Modulus of ");

% E1 comparison
subplot(2,3,1)
scatter(expResults.E1_exp, expResults.E1_num, 60, expResults.wavenumber_over_pi, 'filled', ...
    'MarkerFaceAlpha', 0.7)
hold on
minVal = min([expResults.E1_exp; expResults.E1_num]);
maxVal = max([expResults.E1_exp; expResults.E1_num]);
plot([minVal maxVal], [minVal maxVal], 'k--', 'LineWidth', 2)
xlabel('Experimental E_1 (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Numerical E_1 (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
title(sprintf('Grouped by wavenumber, R² = %.3f', R2_E1), 'FontSize', 14, 'FontWeight', 'bold')
c = colorbar; 
c.Label.String = 'Wavenumber/\pi';
grid on
axis equal tight
set(gca, 'FontSize', 11)

% E2 comparison
subplot(2,3,4)
scatter(expResults.E2_exp, expResults.E2_num, 60, expResults.wavenumber_over_pi, 'filled', ...
    'MarkerFaceAlpha', 0.7)
hold on
minVal = min([expResults.E2_exp; expResults.E2_num]);
maxVal = max([expResults.E2_exp; expResults.E2_num]);
plot([minVal maxVal], [minVal maxVal], 'k--', 'LineWidth', 2)
xlabel('Experimental E_2 (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Numerical E_2 (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
title(sprintf('Grouped by wavenumber, R² = %.3f', R2_E2), 'FontSize', 14, 'FontWeight', 'bold')
c = colorbar; 
c.Label.String = 'Wavenumber/\pi';
grid on
axis equal tight
set(gca, 'FontSize', 11)


% E1 comparison
subplot(2,3,2)
scatter(expResults.E1_exp, expResults.E1_num, 60, expResults.rho, 'filled', ...
    'MarkerFaceAlpha', 0.7)
hold on
minVal = min([expResults.E1_exp; expResults.E1_num]);
maxVal = max([expResults.E1_exp; expResults.E1_num]);
plot([minVal maxVal], [minVal maxVal], 'k--', 'LineWidth', 2)
xlabel('Experimental E_1 (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Numerical E_1 (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
title('Grouped by density', 'FontSize', 14, 'FontWeight', 'bold')
c = colorbar; 
c.Label.String = 'Relative density';
grid on
axis equal tight
set(gca, 'FontSize', 11)

% E2 comparison
subplot(2,3,5)
scatter(expResults.E2_exp, expResults.E2_num, 60, expResults.rho, 'filled', ...
    'MarkerFaceAlpha', 0.7)
hold on
minVal = min([expResults.E2_exp; expResults.E2_num]);
maxVal = max([expResults.E2_exp; expResults.E2_num]);
plot([minVal maxVal], [minVal maxVal], 'k--', 'LineWidth', 2)
xlabel('Experimental E_2 (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Numerical E_2 (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
title('Grouped by density', 'FontSize', 14, 'FontWeight', 'bold')
c = colorbar; 
c.Label.String = 'Relative density';
grid on
axis equal tight
set(gca, 'FontSize', 11)





% E1 comparison
subplot(2,3,3)
scatter(expResults.E1_exp, expResults.E1_num, 60, expResults.theta1, 'filled', ...
    'MarkerFaceAlpha', 0.7)
hold on
minVal = min([expResults.E1_exp; expResults.E1_num]);
maxVal = max([expResults.E1_exp; expResults.E1_num]);
plot([minVal maxVal], [minVal maxVal], 'k--', 'LineWidth', 2)
xlabel('Experimental E_1 (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Numerical E_1 (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
title('Grouped by \theta_1', 'FontSize', 14, 'FontWeight', 'bold')
c = colorbar; 
c.Label.String = '\theta_1';
grid on
axis equal tight
set(gca, 'FontSize', 11)

% E2 comparison
subplot(2,3,6)
scatter(expResults.E2_exp, expResults.E2_num, 60, expResults.theta1, 'filled', ...
    'MarkerFaceAlpha', 0.7)
hold on
minVal = min([expResults.E2_exp; expResults.E2_num]);
maxVal = max([expResults.E2_exp; expResults.E2_num]);
plot([minVal maxVal], [minVal maxVal], 'k--', 'LineWidth', 2)
xlabel('Experimental E_2 (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Numerical E_2 (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
title('Grouped by \theta_1', 'FontSize', 14, 'FontWeight', 'bold')
c = colorbar; 
c.Label.String = '\theta_1';
grid on
axis equal tight
set(gca, 'FontSize', 11)

sgtitle("Numerical Predictions vs Experimental Measurements of Effective Young's Modulus, grouped by input parameters", ...
    'FontSize', 14, 'FontWeight', 'bold')

%% DATA CLEANING: Remove outliers based on "Factor of 10" rule
% We calculate the ratio first to identify the "noisy" experiments
ratio_check = expResults.E1_exp ./ expResults.E1_num;

% Logical mask: Keep only rows where 0.1 < ratio < 10
% This removes cases where Num and Exp are wildly disconnected
keepIdx = (ratio_check > 0.1) & (ratio_check < 10);

% Create a filtered table for regression
expResultsFiltered = expResults(keepIdx, :);

fprintf('Outlier Removal: Removed %d rows out of %d total.\n', ...
    sum(~keepIdx), height(expResults));

% % ADJUSTING THE NUMERICAL RESULTS (Using Filtered Data)
% % Target: The multiplier needed to correct the NN output
% Ratio_E1 = expResults.E1_exp ./ expResults.E1_num;
% 
% % Predictors: focus only on rho and wavenumber
% X = [expResults.rho, expResults.wavenumber_over_pi];

% ADJUSTING THE NUMERICAL RESULTS (Using Filtered Data)
% Target: The multiplier needed to correct the NN output
Ratio_E1 = expResultsFiltered.E1_exp ./ expResultsFiltered.E1_num;

% Predictors: focus only on rho and wavenumber
X = [expResultsFiltered.rho, expResultsFiltered.wavenumber_over_pi];
%X = expResultsFiltered.wavenumber_over_pi;
% Fit a linear model: Ratio = beta0 + beta1*rho + beta2*wavenumber
mdl_correction = fitlm(X, Ratio_E1, 'VarNames', {'density','wavenumber', 'Ratio'});

% Extract the coefficients
b = mdl_correction.Coefficients.Estimate;
fprintf('\n--- Correction Function Coefficients ---\n');
fprintf('Intercept: %.4f\n', b(1));
%fprintf('Rho Coeff: %.4f\n', b(2));
fprintf('Wavenumber Coeff: %.4f\n', b(2));

% 1. Extract the coefficients from your mdl_correction
beta0 = b(1); % Intercept
beta_rho = b(2); % Rho Coefficient
beta_wave = b(3); % Wavenumber Coefficient

% 2. Calculate the Correction Multiplier (K) for every row in the table
% K is the predicted Ratio (Exp/Num) based on the design parameters
K = beta0 + (beta_rho * expResults.rho) + (beta_wave * expResults.wavenumber_over_pi);
%K = beta0 + beta_wave * expResults.wavenumber_over_pi;
% 3. Apply the correction to generate the 'Adjusted' numerical values
% We multiply the original numerical output by the correction factor
expResults.E1_num_adjust = expResults.E1_num .* K;
expResults.E2_num_adjust = expResults.E2_num .* K;
expResults.E3_num_adjust = expResults.E3_num .* K;

% 4. Safety Clip (Optional)
% Since stiffness cannot be negative, we set any potential negative 
% results (from the linear intercept) to zero.
expResults.E1_num_adjust = max(0, expResults.E1_num_adjust);
expResults.E2_num_adjust = max(0, expResults.E2_num_adjust);
expResults.E3_num_adjust = max(0, expResults.E3_num_adjust);

% 5. Print a quick verification of the first few rows
disp('Verification: Original Num vs Adjusted Num vs Experimental');
disp(expResults(1:5, {'E1_num', 'E1_num_adjust', 'E1_exp'}));

%% 4. CALCULATE NEW R² VALUES
R2_E1_original = calcR2(expResults.E1_exp,expResults.E1_num);
R2_E1_adjusted = calcR2(expResults.E1_exp,expResults.E1_num_adjust);

R2_E2_original = calcR2(expResults.E2_exp,expResults.E2_num);
R2_E2_adjusted = calcR2(expResults.E2_exp,expResults.E2_num_adjust);

fprintf('R² Improvement:\n');
fprintf('E1: %.4f → %.4f\n', R2_E1_original, R2_E1_adjusted);
fprintf('E2: %.4f → %.4f\n', R2_E2_original, R2_E2_adjusted);

%% 5. PLOT: Before and After Correction

figure('Position', [100 100 1600 500],'Name',"Effective Young's Modulus after correction");

% E1
subplot(1,2,1)
hold on
% Before correction
scatter(expResults.E1_exp, expResults.E1_num, 40, [0.8 0.8 0.8], 'filled', 'MarkerFaceAlpha', 0.5)
% After correction
scatter(expResults.E1_exp, expResults.E1_num_adjust, 60, 'r', 'filled', 'MarkerFaceAlpha', 0.7)
% Perfect alignment line
plot([0 max(expResults.E1_num)], [0 max(expResults.E1_num)], 'k--', 'LineWidth', 2)
xlabel('Experimental E_1', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Numerical E_1', 'FontSize', 12, 'FontWeight', 'bold')
legend('Original', 'Adjusted', 'Perfect Fit', 'Location', 'best')
title(sprintf('E_1 Correction, R² = %.3f', R2_E1_adjusted), 'FontSize', 14, 'FontWeight', 'bold')
grid on
axis equal tight

% E2
subplot(1,2,2)
hold on
scatter(expResults.E2_exp, expResults.E2_num, 40, [0.8 0.8 0.8], 'filled', 'MarkerFaceAlpha', 0.5)
scatter(expResults.E2_exp, expResults.E2_num_adjust, 60, 'r', 'filled', 'MarkerFaceAlpha', 0.7)
plot([0 max(expResults.E2_num)], [0 max(expResults.E2_num)], 'k--', 'LineWidth', 2)
xlabel('Experimental E_2', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Numerical E_2', 'FontSize', 12, 'FontWeight', 'bold')
legend('Original', 'Adjusted', 'Perfect Fit', 'Location', 'best')
title(sprintf('E_2 Correction, R² = %.3f', R2_E2_adjusted), 'FontSize', 14, 'FontWeight', 'bold')
grid on
axis equal tight

sgtitle('Wavenumber-Based Correction', 'FontSize', 16, 'FontWeight', 'bold')


%% 6. EXPAND: Introduce the data from the larger wavenumbers and see if it is indeed better.

%% 1. LOAD EXPERIMENTAL DATA

fprintf('=== LOADING EXPERIMENTAL DATASET 2 ===\n');

% Read Excel file
expData2 = readtable('ExperimentData.xlsx', 'Sheet', 'Sheet3');
youngModulus = 1600;
% The first row contains headers, extract them

data = expData2(:,:);
data(1:8,6) = data(1:8,7);
% Convert to numeric and create properly structured table
expResults2 = table();
expResults2.rho = data{1:end,1};
expResults2.wavenumber_over_pi = data{1:end,2};
expResults2.theta1 = data{1:end,3};
expResults2.theta2 = data{1:end,4};
expResults2.theta3 = data{1:end,5};
expResults2.E_exp = data{1:end,6};


%Normalize for the Young's Modulus of FormLabs Clear Resin V4
expResults2.E_exp = expResults2.E_exp / 1600;


fprintf('Loaded %d experimental data points\n', height(expResults2));
fprintf('Unique densities: %s\n', mat2str(unique(expResults2.rho)'));
fprintf('Unique wavenumbers/π: %s\n', mat2str(unique(expResults2.wavenumber_over_pi)'));
fprintf('Unique θ₁ values: %s\n', mat2str(unique(expResults2.theta1)'));
fprintf('\n');



%% 3. GENERATE NUMERICAL PREDICTIONS

fprintf('=== GENERATING NUMERICAL PREDICTIONS ===\n');

% Prepare input for neural network
inputParams2 = [expResults2.rho, expResults2.theta1, expResults2.theta2, expResults2.theta3];

% Normalize if normalization parameters exist
if ~isempty(normParams)
    inputMin = normParams.input.min;
    inputMax = normParams.input.max;
    inputParams2_norm = (inputParams2 - inputMin) ./ (inputMax - inputMin);
    fprintf('Input parameters normalized\n');
else
    inputParams2_norm = inputParams2;
    fprintf('No normalization applied\n');
end

% Predict using forward network
stiffnessPred2_norm = predict(fNet, inputParams2_norm);

% Denormalize predictions if needed
if ~isempty(normParams) && isfield(normParams, 'output')
    outputMin = normParams.output.min;
    outputMax = normParams.output.max;
    stiffnessPred2 = stiffnessPred2_norm .* (outputMax - outputMin) + outputMin;
    fprintf('Output denormalized\n');
else
    stiffnessPred2 = stiffnessPred2_norm;
end

% Note: Kumar et al. uses wavenumber β = 10π/l for l×l×l elements
% The stiffness components are C1111, C1122, C1133, C2222, C2233, C3333, C2323, C3131, C1212
% For orthotropic materials, E1 ≈ 1/S1111 where S = inv(C) (compliance matrix)

% Extract diagonal stiffness components (approximation for Young's moduli)
C1111_num = stiffnessPred2(:,1);
C2222_num = stiffnessPred2(:,4);
C3333_num = stiffnessPred2(:,6);

% Store predictions
expResults2.E_num(9:16) = C1111_num(9:16);
expResults2.E_num(1:8) = C2222_num(1:8);


fprintf('Generated %d numerical predictions\n', height(expResults2));
fprintf('\n');

%% 4. CALCULATE ERRORS

fprintf('=== CALCULATING ERRORS ===\n');

% Absolute errors
expResults2.error_E_abs = expResults2.E_exp - expResults2.E_num;


% Relative errors (percentage)
expResults2.error_E_rel = 100 * (expResults2.E_exp - expResults2.E_num) ./ expResults2.E_exp;

% Overall metrics
fprintf('Overall Error Statistics:\n');
fprintf('E - Mean Abs Error: %.2f, Mean Rel Error: %.2f%%\n', ...
    mean(abs(expResults2.error_E_abs)), mean(abs(expResults2.error_E_rel)));

% R² calculations
R2_E = 1 - sum(expResults2.error_E_abs.^2) / sum((expResults2.E_exp - mean(expResults2.E_exp)).^2);

fprintf('\nR² Scores:\n');
fprintf('E1: R² = %.4f\n', R2_E);
fprintf('\n');

%% 5. VISUALIZATION: OVERALL COMPARISON, COLOURED ATTRIBUTES

figure('Position', [50 50 1600 500], 'Name', 'Overall Comparison');

% E1 comparison
subplot(1,3,1)
scatter(expResults2.E_exp, expResults2.E_num, 60, expResults2.wavenumber_over_pi, 'filled', ...
    'MarkerFaceAlpha', 0.7)
hold on
minVal = min([expResults2.E_exp; expResults2.E_num]);
maxVal = max([expResults2.E_exp; expResults2.E_num]);
plot([0 0.6], [0 0.6], 'k--', 'LineWidth', 2)
xlabel('Experimental E (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Numerical E (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
title(sprintf('Grouped by wavenumber, R² = %.3f', R2_E), 'FontSize', 14, 'FontWeight', 'bold')
c = colorbar; 
c.Label.String = 'Wavenumber/\pi';
grid on
axis equal tight
set(gca, 'FontSize', 11)

% E1 comparison
subplot(1,3,2)
scatter(expResults2.E_exp, expResults2.E_num, 60, expResults2.rho, 'filled', ...
    'MarkerFaceAlpha', 0.7)
hold on
minVal = min([expResults2.E_exp; expResults2.E_num]);
maxVal = max([expResults2.E_exp; expResults2.E_num]);
plot([0 0.6], [0 0.6], 'k--', 'LineWidth', 2)
xlabel('Experimental E (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Numerical E (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
title(sprintf('Grouped by density'), 'FontSize', 14, 'FontWeight', 'bold')
c = colorbar; 
c.Label.String = 'Density';
grid on
axis equal tight
set(gca, 'FontSize', 11)

% E1 comparison
subplot(1,3,3)
scatter(expResults2.E_exp, expResults2.E_num, 60, expResults2.theta1, 'filled', ...
    'MarkerFaceAlpha', 0.7)
hold on
minVal = min([expResults2.E_exp; expResults2.E_num]);
maxVal = max([expResults2.E_exp; expResults2.E_num]);
plot([0 0.6], [0 0.6], 'k--', 'LineWidth', 2)
xlabel('Experimental E (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Numerical E (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
title('Grouped by \theta_1', 'FontSize', 14, 'FontWeight', 'bold')
c = colorbar; 
c.Label.String = '\theta_1';
grid on
axis equal tight
set(gca, 'FontSize', 11)

K = beta0 + (beta_wave * expResults2.wavenumber_over_pi);


% 
% % E1 comparison
% subplot(2,3,2)
% scatter(expResults.E1_exp, expResults.E1_num, 60, expResults.rho, 'filled', ...
%     'MarkerFaceAlpha', 0.7)
% hold on
% minVal = min([expResults.E1_exp; expResults.E1_num]);
% maxVal = max([expResults.E1_exp; expResults.E1_num]);
% plot([minVal maxVal], [minVal maxVal], 'k--', 'LineWidth', 2)
% xlabel('Experimental E_1 (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
% ylabel('Numerical E_1 (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
% title('Grouped by density', 'FontSize', 14, 'FontWeight', 'bold')
% c = colorbar; 
% c.Label.String = 'Relative density';
% grid on
% axis equal tight
% set(gca, 'FontSize', 11)
% 
% % E2 comparison
% subplot(2,3,5)
% scatter(expResults.E2_exp, expResults.E2_num, 60, expResults.rho, 'filled', ...
%     'MarkerFaceAlpha', 0.7)
% hold on
% minVal = min([expResults.E2_exp; expResults.E2_num]);
% maxVal = max([expResults.E2_exp; expResults.E2_num]);
% plot([minVal maxVal], [minVal maxVal], 'k--', 'LineWidth', 2)
% xlabel('Experimental E_2 (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
% ylabel('Numerical E_2 (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
% title('Grouped by density', 'FontSize', 14, 'FontWeight', 'bold')
% c = colorbar; 
% c.Label.String = 'Relative density';
% grid on
% axis equal tight
% set(gca, 'FontSize', 11)
% 
% 
% 
% 
% 
% % E1 comparison
% subplot(2,3,3)
% scatter(expResults.E1_exp, expResults.E1_num, 60, expResults.theta1, 'filled', ...
%     'MarkerFaceAlpha', 0.7)
% hold on
% minVal = min([expResults.E1_exp; expResults.E1_num]);
% maxVal = max([expResults.E1_exp; expResults.E1_num]);
% plot([minVal maxVal], [minVal maxVal], 'k--', 'LineWidth', 2)
% xlabel('Experimental E_1 (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
% ylabel('Numerical E_1 (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
% title('Grouped by \theta_1', 'FontSize', 14, 'FontWeight', 'bold')
% c = colorbar; 
% c.Label.String = '\theta_1';
% grid on
% axis equal tight
% set(gca, 'FontSize', 11)
% 
% % E2 comparison
% subplot(2,3,6)
% scatter(expResults.E2_exp, expResults.E2_num, 60, expResults.theta1, 'filled', ...
%     'MarkerFaceAlpha', 0.7)
% hold on
% minVal = min([expResults.E2_exp; expResults.E2_num]);
% maxVal = max([expResults.E2_exp; expResults.E2_num]);
% plot([minVal maxVal], [minVal maxVal], 'k--', 'LineWidth', 2)
% xlabel('Experimental E_2 (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
% ylabel('Numerical E_2 (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
% title('Grouped by \theta_1', 'FontSize', 14, 'FontWeight', 'bold')
% c = colorbar; 
% c.Label.String = '\theta_1';
% grid on
% axis equal tight
% set(gca, 'FontSize', 11)
% 
% sgtitle('Numerical Predictions vs Experimental Measurements, grouped by input parameters', ...
%     'FontSize', 16, 'FontWeight', 'bold')

%%
% expResults = table();
% expResults.rho = data{1:end,1};
% expResults.wavenumber_over_pi = data{1:end,2};
% expResults.theta1 = data{1:end,3};
% expResults.theta2 = data{1:end,4};
% expResults.theta3 = data{1:end,5};
% expResults.E1_exp = data{1:end,6};
% expResults.E2_exp = data{1:end,7};
% expResults.E3_exp = data{1:end,8};
combinedResults = table();
combinedResults.rho = [expResults.rho; expResults2.rho];
combinedResults.wavenumber_over_pi = [expResults.wavenumber_over_pi; expResults2.wavenumber_over_pi];
combinedResults.theta1 = [expResults.theta1; expResults2.theta1];
combinedResults.theta2 = [expResults.theta2; expResults2.theta2];
combinedResults.theta3 = [expResults.theta3; expResults2.theta3];
combinedResults.E_exp = [expResults.E1_exp; expResults2.E_exp];
combinedResults.E_num = [expResults.E1_num; expResults2.E_num];

%%
figure('Position', [50 50 1600 500], 'Name', 'Overall Comparison');

% E1 comparison
subplot(1,3,1)
scatter(combinedResults.E_exp, combinedResults.E_num, 60, combinedResults.wavenumber_over_pi, 'filled', ...
    'MarkerFaceAlpha', 0.7)
hold on
minVal = min([expResults2.E_exp; expResults2.E_num]);
maxVal = max([expResults2.E_exp; expResults2.E_num]);
plot([0 0.6], [0 0.6], 'k--', 'LineWidth', 2)
xlabel('Experimental E (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Numerical E (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
title(sprintf('Grouped by wavenumber, R² = %.3f', R2_E), 'FontSize', 14, 'FontWeight', 'bold')
c = colorbar; 
c.Label.String = 'Wavenumber/\pi';
grid on
axis equal tight
set(gca, 'FontSize', 11)

% E1 comparison
subplot(1,3,2)
scatter(combinedResults.E_exp, combinedResults.E_num, 60, combinedResults.rho, 'filled', ...
    'MarkerFaceAlpha', 0.7)
hold on
minVal = min([expResults2.E_exp; expResults2.E_num]);
maxVal = max([expResults2.E_exp; expResults2.E_num]);
plot([0 0.6], [0 0.6], 'k--', 'LineWidth', 2)
xlabel('Experimental E (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Numerical E (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
title(sprintf('Grouped by density'), 'FontSize', 14, 'FontWeight', 'bold')
c = colorbar; 
c.Label.String = 'Density';
grid on
axis equal tight
set(gca, 'FontSize', 11)

% E1 comparison
subplot(1,3,3)
scatter(combinedResults.E_exp, combinedResults.E_num, 60, combinedResults.theta1, 'filled', ...
    'MarkerFaceAlpha', 0.7)
hold on
minVal = min([expResults2.E_exp; expResults2.E_num]);
maxVal = max([expResults2.E_exp; expResults2.E_num]);
plot([0 0.6], [0 0.6], 'k--', 'LineWidth', 2)
xlabel('Experimental E (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Numerical E (MPa)', 'FontSize', 12, 'FontWeight', 'bold')
title('Grouped by \theta_1', 'FontSize', 14, 'FontWeight', 'bold')
c = colorbar; 
c.Label.String = '\theta_1';
grid on
axis equal tight
set(gca, 'FontSize', 11)

%% FUNCTIONS

function R2 = calcR2(input1, input2)
    R2 = 1-sum((input1-input2).^2 / sum((input1-mean(input1)).^2));
end
