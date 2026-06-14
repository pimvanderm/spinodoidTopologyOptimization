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
csvPath = fullfile(defaultFolder, 'data', 'temp', 'Homogenization_Stiffness_Ledger.csv');

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

% Create the true physical weight vector for the X-axis coordinate mapping
% Loop 0 (Index 1) maps to weight1 = 0.0; Loop 10 (Index 11) maps to weight1 = 1.0
weightVector = 0.0:0.1:1.0; 

% Define the custom string tick labels for the transition breakdown
% Format: [Weight Section 1 %] / [Weight Section 2 %]
tickLabels = {'0/100', '10/90', '20/80', '30/70', '40/60', '50/50', ...
              '60/40', '70/30', '80/20', '90/10', '100/0'};

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
% Plot 1: E1 Modulus Trajectory vs Hypothesis (Mapped to Weight Coordinate)
figure('Name', 'E1 Modulus Tracking', 'Color', 'w', 'Position', [100, 100, 850, 500]);
hold on;
plot(weightVector, E1_seeds(1,:), 'Color', [0.7 0.7 0.7], 'LineStyle', '--', 'DisplayName', 'Seed 1');
plot(weightVector, E1_seeds(2,:), 'Color', [0.7 0.7 0.7], 'LineStyle', ':', 'DisplayName', 'Seed 2');
plot(weightVector, E1_seeds(3,:), 'Color', [0.7 0.7 0.7], 'LineStyle', '-.', 'DisplayName', 'Seed 3');

plot(weightVector, E1_avg, 'k-o', 'LineWidth', 2.5, 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Averaged E_1');
yline(E1_hyp, 'r-', 'LineWidth', 2, 'DisplayName', 'Hypothesis [15,15,15]');

title('Directional Young''s Modulus (E_1) Across Transition Phase');
xlabel('Transition Weight (Columnar 1 % / Lamellar 2 %)');
ylabel('Stiffness E_1 (MPa)');
xlim([0 1]);

% Apply custom x-axis breaks and structural breakdown strings
set(gca, 'XTick', weightVector, 'XTickLabel', tickLabels);
legend('Location', 'best');
grid on; hold off;

% Plot 2: Inverse NN Prediction based on Averaged Stiffness (Mapped to Weight Coordinate)
figure('Name', 'Inverse NN - Averaged', 'Color', 'w', 'Position', [980, 100, 850, 500]);
hold on;
plot(weightVector, inputArrays_avg(:, 2), '-o', 'LineWidth', 2, 'DisplayName', '\theta_1 (X-axis)');
plot(weightVector, inputArrays_avg(:, 3), '-s', 'LineWidth', 2, 'DisplayName', '\theta_2 (Y-axis)');
plot(weightVector, inputArrays_avg(:, 4), '-^', 'LineWidth', 2, 'DisplayName', '\theta_3 (Z-axis)');

title('Averaged Recovered Conical Half-Angles (\theta)');
ylabel('Angle (Degrees)');
xlabel('Transition Weight (Columnar %/ Lamellar %)');
ylim([0 90]);
xlim([0 1]);

% Apply custom x-axis breaks and structural breakdown strings
set(gca, 'XTick', weightVector, 'XTickLabel', tickLabels);
legend('Location', 'best');
grid on; hold off;