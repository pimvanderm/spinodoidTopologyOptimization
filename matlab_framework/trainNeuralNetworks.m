%% Inverse-Designed Spinodoid Metamaterials - Neural Network Implementation
% Based on Kumar et al. (2020) npj Computational Materials
% Paper: https://doi.org/10.1038/s41524-020-0341-6

clear; clc; close all;

%% 1. READ AND VISUALIZE DATA (with normalization)
fprintf('Loading datasets...\n');
trainData = readtable("spinodoid-train.csv");
testData  = readtable("spinodoid-test.csv");

% Convert to arrays
trainData = table2array(trainData);
testData  = table2array(testData);

% Split inputs and outputs
ThetaTrain_raw = trainData(:,1:4);   % Raw design parameters
STrain_raw     = trainData(:,5:end); % Raw stiffness

ThetaTest_raw = testData(:,1:4);
STest_raw     = testData(:,5:end);

% NORMALIZE INPUTS
% Define normalization parameters based on your data ranges
% [rho, theta1, theta2, theta3]
inputMin = [0, 0, 0, 0];
inputMax = [1, 90, 90, 90];

% Normalize to [0, 1]
ThetaTrain = (ThetaTrain_raw - inputMin) ./ (inputMax - inputMin);
ThetaTest = (ThetaTest_raw - inputMin) ./ (inputMax - inputMin);

% Optional: Normalize outputs too (often helps)
outputMin = min(STrain_raw, [], 1);
outputMax = max(STrain_raw, [], 1);

STrain = (STrain_raw - outputMin) ./ (outputMax - outputMin);
STest = (STest_raw - outputMin) ./ (outputMax - outputMin);

% Store normalization parameters for later use
normParams.input.min = inputMin;
normParams.input.max = inputMax;
normParams.output.min = outputMin;
normParams.output.max = outputMax;

save('normalizationParams.mat', 'normParams');

fprintf('Training samples: %d\n', size(ThetaTrain,1));
fprintf('Test samples: %d\n', size(ThetaTest,1));
fprintf('\nNormalization applied:\n');
fprintf('  Inputs normalized to [0, 1]\n');
fprintf('  Outputs normalized to [0, 1]\n\n');

%% 2. SETUP FORWARD NEURAL NETWORK (f-NN)
% Architecture based on Fig. 4 of Kumar et al.
% Input: 4 design parameters -> Output: 9 stiffness components
% 6 hidden layers: [128, 128, 64, 64, 32, 32]

fprintf('\n=== FORWARD NETWORK SETUP ===\n');

layersF = [
    featureInputLayer(4, "Name", "input", "Normalization", "zscore")
    
    fullyConnectedLayer(128, "Name", "fc1")
    reluLayer("Name", "relu1")
    
    fullyConnectedLayer(128, "Name", "fc2")
    reluLayer("Name", "relu2")
    
    fullyConnectedLayer(64, "Name", "fc3")
    reluLayer("Name", "relu3")
    
    fullyConnectedLayer(64, "Name", "fc4")
    reluLayer("Name", "relu4")
    
    fullyConnectedLayer(32, "Name", "fc5")
    reluLayer("Name", "relu5")
    
    fullyConnectedLayer(32, "Name", "fc6")
    reluLayer("Name", "relu6")
    
    fullyConnectedLayer(size(STrain,2), "Name", "output")
    regressionLayer("Name", "regression")
];

fprintf('Forward network architecture:\n');
fprintf('  Input: %d parameters\n', 4);
fprintf('  Hidden layers: 128 -> 128 -> 64 -> 64 -> 32 -> 32\n');
fprintf('  Output: %d stiffness components\n\n', size(STrain,2));

%% 3. TRAIN FORWARD MODEL

fprintf('=== TRAINING FORWARD MODEL ===\n');

% Create a function to store losses
lossHistory = [];

function stop = recordLoss(info)
    persistent history
    if info.State == "start"
        history = [];
    elseif info.State == "iteration"
        history(end+1) = info.TrainingLoss;
        assignin('base', 'lossHistory', history);
    end
    stop = false;
end

optionsF = trainingOptions("adam", ...
    MaxEpochs=200, ...
    MiniBatchSize=64, ...
    InitialLearnRate=1e-3, ...
    LearnRateSchedule='piecewise', ...
    LearnRateDropPeriod=30, ...
    LearnRateDropFactor=0.5, ...
    Shuffle="every-epoch", ...
    OutputFcn=@recordLoss, ...
    ValidationData={ThetaTest, STest}, ...
    ValidationFrequency=50, ...
    Plots="training-progress", ...
    Verbose=true, ...
    VerboseFrequency=50);

fNet = trainNetwork(ThetaTrain, STrain, layersF, optionsF);
save(fNet)

%% 4. EVALUATE FORWARD MODEL

fprintf('\n=== FORWARD MODEL EVALUATION ===\n');

% Predictions on test set (normalized)
SPred_test_norm = predict(fNet, ThetaTest);

% Denormalize for evaluation
SPred_test = SPred_test_norm .* (outputMax - outputMin) + outputMin;
STest_denorm = STest_raw;  % Use original raw test data

% Calculate R² for each stiffness component
R2_forward = zeros(1, size(STest_denorm,2));
figure('Position', [100 100 1400 500]);

for i = 1:size(STest_denorm,2)
    SS_res = sum((STest_denorm(:,i) - SPred_test(:,i)).^2);
    SS_tot = sum((STest_denorm(:,i) - mean(STest_denorm(:,i))).^2);
    R2_forward(i) = 1 - SS_res/SS_tot;
    
    subplot(3, 3, i)
    scatter(STest_denorm(:,i), SPred_test(:,i), 20, 'filled', 'MarkerFaceAlpha', 0.6)
    hold on
    plot([min(STest_denorm(:,i)), max(STest_denorm(:,i))], ...
         [min(STest_denorm(:,i)), max(STest_denorm(:,i))], 'r--', 'LineWidth', 2)
    xlabel('True Stiffness', 'FontSize', 10)
    ylabel('Predicted Stiffness', 'FontSize', 10)
    title(sprintf('Component %d (R² = %.4f)', i, R2_forward(i)), 'FontSize', 11)
    grid on
    axis equal tight
end

sgtitle('Forward Model: Predicted vs True Stiffness (Denormalized)', 'FontSize', 16, 'FontWeight', 'bold')

%% 5. SETUP INVERSE NEURAL NETWORK (i-NN)
% Architecture: mirrors forward network but reversed
% Input: 9 stiffness components -> Output: 4 design parameters

fprintf('=== INVERSE NETWORK SETUP ===\n');

layersI = [
    featureInputLayer(size(STrain,2), "Name", "input")
    
    fullyConnectedLayer(100, "Name", "fc1")
    reluLayer("Name", "relu1")
    
    fullyConnectedLayer(100, "Name", "fc2")
    reluLayer("Name", "relu2")
    
    fullyConnectedLayer(100, "Name", "fc3")
    reluLayer("Name", "relu3")
    
    fullyConnectedLayer(100, "Name", "fc4")
    reluLayer("Name", "relu4")
    
    fullyConnectedLayer(100, "Name", "fc5")
    reluLayer("Name", "relu5")
    
    fullyConnectedLayer(100, "Name", "fc6")
    reluLayer("Name", "relu6")
    
    fullyConnectedLayer(size(ThetaTrain,2), "Name", "output")
];

% Create dlnetwork for custom training
% Remove regression layer from forward network for reconstruction
layersF_noReg = fNet.Layers(1:end-1);
dlnetF = dlnetwork(layersF_noReg);

% Freeze forward network by accessing the Learnables table
% This works in MATLAB R2020b and later
if ~isempty(dlnetF.Learnables)
    dlnetF.Learnables.Value(:) = dlnetF.Learnables.Value;
end

dlnetI = dlnetwork(layersI);

fprintf('Inverse network architecture:\n');
fprintf('  Input: %d stiffness components\n', size(STrain,2));
fprintf('  Hidden layers: 6 x 100\n');
fprintf('  Output: %d parameters\n\n', 4);

%% 6. TRAIN INVERSE MODEL (with normalized data)

fprintf('=== TRAINING INVERSE MODEL ===\n');

% Training hyperparameters
numEpochs = 200;
miniBatchSize = 64;
lambda = 0.5;
learnRate = 1e-4;

% Prepare normalized data
numObs = size(STrain,1);
numBatches = floor(numObs/miniBatchSize);

% Adam optimizer parameters
iteration = 0;
averageGrad = [];
averageSqGrad = [];

% Loss tracking
lossHistory = struct();
lossHistory.total = [];
lossHistory.reconstruction = [];
lossHistory.parameter = [];

fprintf('Training with NORMALIZED data\n');
fprintf('Training configuration:\n');
fprintf('  Epochs: %d\n', numEpochs);
fprintf('  Batch size: %d\n', miniBatchSize);
fprintf('  Initial learning rate: %.1e\n', learnRate);
fprintf('  Lambda (parameter loss weight): %.2f\n\n', lambda);

% Training loop
fprintf('Progress:\n');
fig = figure('Position', [100 100 1200 400]);

for epoch = 1:numEpochs
    % Shuffle data
    idx = randperm(numObs);
    
    % Reduce lambda after initial epochs
    if epoch > 40
        lambda = 0;
    end
    
    epochLoss = 0;
    epochReconLoss = 0;
    epochParamLoss = 0;
    
    for i = 1:numBatches
        iteration = iteration + 1;
        
        % Get mini-batch (using NORMALIZED data)
        batchIdx = idx((i-1)*miniBatchSize+1 : i*miniBatchSize);
        SBatch = STrain(batchIdx,:)';  % Normalized stiffness
        ThetaBatch = ThetaTrain(batchIdx,:)';  % Normalized parameters
        
        % Convert to dlarray
        dlS = dlarray(single(SBatch), 'CB');
        dlThetaTrue = dlarray(single(ThetaBatch), 'CB');
        
        % Compute loss and gradients
        [loss, gradients, reconLoss, paramLoss] = dlfeval( ...
            @inverseLossFunction, dlnetI, dlnetF, dlS, dlThetaTrue, lambda);
        
        % Update network using Adam
        [dlnetI, averageGrad, averageSqGrad] = adamupdate( ...
            dlnetI, gradients, averageGrad, averageSqGrad, iteration, learnRate);
        
        % Accumulate losses
        epochLoss = epochLoss + extractdata(loss);
        epochReconLoss = epochReconLoss + extractdata(reconLoss);
        epochParamLoss = epochParamLoss + extractdata(paramLoss);
    end
    
    % Average losses
    lossHistory.total(epoch) = epochLoss / numBatches;
    lossHistory.reconstruction(epoch) = epochReconLoss / numBatches;
    lossHistory.parameter(epoch) = epochParamLoss / numBatches;
    
    % Display progress
    if mod(epoch, 10) == 0 || epoch == 1
        fprintf('Epoch %3d/%d | Loss: %.4e | Recon: %.4e | Param: %.4e | λ: %.3f\n', ...
            epoch, numEpochs, lossHistory.total(epoch), ...
            lossHistory.reconstruction(epoch), lossHistory.parameter(epoch), lambda);
    end
    
    % Update plot every 5 epochs
    if mod(epoch, 5) == 0
        figure(fig);
        subplot(1,2,1)
        semilogy(1:epoch, lossHistory.total(1:epoch), 'b-', 'LineWidth', 2)
        hold on
        semilogy(1:epoch, lossHistory.reconstruction(1:epoch), 'r--', 'LineWidth', 1.5)
        semilogy(1:epoch, lossHistory.parameter(1:epoch), 'g:', 'LineWidth', 1.5)
        hold off
        xlabel('Epoch', 'FontSize', 12)
        ylabel('Loss (log scale)', 'FontSize', 12)
        title('Training Loss', 'FontSize', 14)
        legend('Total', 'Reconstruction', 'Parameter', 'Location', 'best')
        grid on
        
        subplot(1,2,2)
        plot(1:epoch, lossHistory.total(1:epoch), 'b-', 'LineWidth', 2)
        xlabel('Epoch', 'FontSize', 12)
        ylabel('Total Loss', 'FontSize', 12)
        title('Total Loss (Linear Scale)', 'FontSize', 14)
        grid on
        
        drawnow
    end
end

fprintf('\nTraining complete!\n\n');
save dlnetF
save dlnetI

%% 7. EVALUATE INVERSE MODEL

fprintf('=== INVERSE MODEL EVALUATION ===\n');

% Predict design parameters for test set (normalized inputs and outputs)
dlSTest = dlarray(single(STest'), 'CB');
ThetaPred_test_norm = extractdata(forward(dlnetI, dlSTest))';

% Denormalize predicted parameters
ThetaPred_test = ThetaPred_test_norm .* (inputMax - inputMin) + inputMin;

% Reconstruct stiffness using forward network (still normalized)
SRecon_test_norm = predict(fNet, ThetaPred_test_norm);

% Denormalize reconstructed stiffness
SRecon_test = SRecon_test_norm .* (outputMax - outputMin) + outputMin;

% Calculate R² for reconstruction (using denormalized values)
R2_inverse_recon = zeros(1, size(STest_raw,2));
for i = 1:size(STest_raw,2)
    SS_res = sum((STest_raw(:,i) - SRecon_test(:,i)).^2);
    SS_tot = sum((STest_raw(:,i) - mean(STest_raw(:,i))).^2);
    R2_inverse_recon(i) = 1 - SS_res/SS_tot;
end

% Visualization
figure('Position', [100 100 1400 900]);

% Stiffness reconstruction (denormalized)
for i = 1:9
    subplot(4, 3, i)
    scatter(STest_raw(:,i), SRecon_test(:,i), 20, 'filled', 'MarkerFaceAlpha', 0.6)
    hold on
    plot([min(STest_raw(:,i)), max(STest_raw(:,i))], ...
         [min(STest_raw(:,i)), max(STest_raw(:,i))], 'r--', 'LineWidth', 2)
    xlabel('Queried Stiffness', 'FontSize', 10)
    ylabel('Reconstructed Stiffness', 'FontSize', 10)
    title(sprintf('C_{%d} (R² = %.4f)', i, R2_inverse_recon(i)), 'FontSize', 11)
    grid on
    axis equal tight
end

% Parameter predictions (denormalized)
paramNames = {'\rho', '\theta_1 (°)', '\theta_2 (°)', '\theta_3 (°)'};
for i = 1:4
    subplot(4, 3, 9+i)
    scatter(ThetaTest_raw(:,i), ThetaPred_test(:,i), 20, 'filled', 'MarkerFaceAlpha', 0.6)
    hold on
    plot([min(ThetaTest_raw(:,i)), max(ThetaTest_raw(:,i))], ...
         [min(ThetaTest_raw(:,i)), max(ThetaTest_raw(:,i))], 'r--', 'LineWidth', 2)
    xlabel(['True ' paramNames{i}], 'FontSize', 10)
    ylabel(['Predicted ' paramNames{i}], 'FontSize', 10)
    
    % R² for parameters
    SS_res = sum((ThetaTest_raw(:,i) - ThetaPred_test(:,i)).^2);
    SS_tot = sum((ThetaTest_raw(:,i) - mean(ThetaTest_raw(:,i))).^2);
    R2_param = 1 - SS_res/SS_tot;
    
    title(sprintf('%s (R² = %.3f*)', paramNames{i}, R2_param), 'FontSize', 11)
    grid on
    axis equal tight
end

sgtitle('Inverse Model: Reconstruction and Parameter Prediction (Denormalized)', ...
    'FontSize', 16, 'FontWeight', 'bold')

%% 8. DEMONSTRATION: INVERSE DESIGN EXAMPLES

fprintf('=== INVERSE DESIGN DEMONSTRATIONS ===\n\n');

% Example 1: Design for specific test case
idx_demo = round(size(ThetaTest,1)*rand);
S_target = STest_raw(idx_demo,:)';
S_target_norm = normalizeOutputs(S_target,normParams);
Theta_true = ThetaTest_raw(idx_demo,:);
Theta_true_norm = normalizeInputs(Theta_true,normParams);

dlS_target = dlarray(single(S_target), 'CB');
Theta_predicted_norm = extractdata(forward(dlnetI, dlS_target))';
S_reconstructed_norm = predict(fNet, Theta_predicted_norm);
S_reconstructed = denormalizeOutputs(S_reconstructed_norm,normParams);
Theta_predicted = denormalizeInputs(Theta_predicted_norm,normParams);

fprintf('Example 1: Test sample #%d\n', idx_demo);
fprintf('True parameters:      [%.3f, %.1f°, %.1f°, %.1f°]\n', ...
    Theta_true(1), Theta_true(2), Theta_true(3), Theta_true(4));
fprintf('Predicted parameters: [%.3f, %.1f°, %.1f°, %.1f°]\n', ...
    Theta_predicted(1), Theta_predicted(2), Theta_predicted(3), Theta_predicted(4));
fprintf('True properties:      [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n', S_target);
fprintf('Predicted properties:      [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n', S_reconstructed);
fprintf('Stiffness reconstruction error: %.2e\n\n', norm(S_target - S_reconstructed'));

%% Example 2: Multiple random test cases
num_examples = 5;
fprintf('Example 2: Reconstruction quality for %d random samples\n', num_examples);
fprintf('%-8s %-25s %-15s\n', 'Sample', 'Relative Error', 'Max Component Error');
fprintf('%s\n', repmat('-', 1, 50));

for i = 1:num_examples
    idx = randi(size(STest,1));
    S_query = STest(idx,:)';
    
    dlS_query = dlarray(single(S_query), 'CB');
    Theta_pred = extractdata(forward(dlnetI, dlS_query))';
    S_recon = predict(fNet, Theta_pred);
    
    rel_error = norm(S_query - S_recon') / norm(S_query);
    max_comp_error = max(abs(S_query - S_recon'));
    
    fprintf('%-8d %-25.4e %-15.4e\n', idx, rel_error, max_comp_error);
end

fprintf('\n=== COMPLETE ===\n');

%% DEMONSTRATION: Design for Target Stiffness

% Load normalization parameters
load('normalizationParams.mat', 'normParams');

% Target stiffness (raw physical values)
S_target_raw = STest_raw(42,:)';

% Normalize
S_target_norm = normalizeOutputs(S_target_raw', normParams)';

% Predict (inverse model works with normalized data)
dlS_target = dlarray(single(S_target_norm), 'CB');
Theta_pred_norm = extractdata(forward(dlnetI, dlS_target))';

% Denormalize prediction
Theta_pred_raw = denormalizeInputs(Theta_pred_norm, normParams);

fprintf('Target stiffness query (sample #42):\n');
fprintf('Predicted parameters:\n');
fprintf('  ρ = %.3f\n', Theta_pred_raw(1));
fprintf('  θ₁ = %.1f°\n', Theta_pred_raw(2));
fprintf('  θ₂ = %.1f°\n', Theta_pred_raw(3));
fprintf('  θ₃ = %.1f°\n', Theta_pred_raw(4));

%% HELPER FUNCTION: Inverse Loss Function

%% HELPER FUNCTIONS

function normalized = normalizeInputs(raw, normParams)
    % Normalize design parameters to [0, 1]
    normalized = (raw - normParams.input.min) ./ (normParams.input.max - normParams.input.min);
end

function denormalized = denormalizeInputs(normalized, normParams)
    % Denormalize design parameters back to original scale
    denormalized = normalized .* (normParams.input.max - normParams.input.min) + normParams.input.min;
end

function normalized = normalizeOutputs(raw, normParams)
    % Normalize stiffness to [0, 1]
    normalized = (raw - normParams.output.min) ./ (normParams.output.max - normParams.output.min);
end

function denormalized = denormalizeOutputs(normalized, normParams)
    % Denormalize stiffness back to original scale
    denormalized = normalized .* (normParams.output.max - normParams.output.min) + normParams.output.min;
end

function [loss, gradients, reconLoss, paramLoss] = inverseLossFunction(dlnetI, dlnetF, dlS, dlThetaTrue, lambda)
    % Forward pass through inverse network
    dlThetaPred = forward(dlnetI, dlS);
    
    % Reconstruct stiffness through frozen forward network
    dlSRecon = forward(dlnetF, dlThetaPred);
    
    % Reconstruction loss (primary objective)
    reconLoss = mean((dlSRecon - dlS).^2, 'all');
    
    % Parameter prediction loss (regularization)
    paramLoss = mean((dlThetaPred - dlThetaTrue).^2, 'all');
    
    % Total loss
    loss = reconLoss + lambda * paramLoss;
    
    % Compute gradients with respect to inverse network parameters
    gradients = dlgradient(loss, dlnetI.Learnables);
end

%% Example: Young's Modulus Surface

% --- Input Data ---
% Assuming 'output' comes from your dlnetF
% Indices: 1:C11, 2:C12, 3:C13, 4:C22, 5:C23, 6:C33, 7:C44, 8:C55, 9:C66
input = [0.293, 34.9, 51.9, 21.2]
prop = forward(dlnetF,input);

% 1. Construct the Stiffness Matrix (C)
C = zeros(6,6);
C(1,1) = prop(1); C(1,2) = prop(2); C(1,3) = prop(3);
C(2,1) = prop(2); C(2,2) = prop(4); C(2,3) = prop(5);
C(3,1) = prop(3); C(3,2) = prop(5); C(3,3) = prop(6);
C(4,4) = prop(7); C(5,5) = prop(8); C(6,6) = prop(9);

% 2. Calculate the Compliance Matrix (S)
% Directional Young's Modulus is defined via the compliance components
S = inv(C);

% 3. Generate Spherical Grid
[phi, theta] = meshgrid(linspace(0, 2*pi, 100), linspace(0, pi, 100));

% Directional unit vectors
n1 = sin(theta) .* cos(phi);
n2 = sin(theta) .* sin(phi);
n3 = cos(theta);

% 4. Calculate Directional Young's Modulus E(n)
% Using the orthotropic expansion of the 4th order compliance tensor
invE = S(1,1)*n1.^4 + S(2,2)*n2.^4 + S(3,3)*n3.^4 ...
     + (2*S(1,2) + S(6,6))*n1.^2 .* n2.^2 ...
     + (2*S(1,3) + S(5,5))*n1.^2 .* n3.^2 ...
     + (2*S(2,3) + S(4,4))*n2.^2 .* n3.^2;

E = 1 ./ invE;

% 5. Convert to Cartesian Coordinates for Plotting
X = E .* n1;
Y = E .* n2;
Z = E .* n3;

% 6. Visualization
figure('Color', 'w');
s = surf(X, Y, Z, E);
shading interp; 
colormap(jet);
colorbar;
axis equal;
grid on;
xlabel('E_1'); ylabel('E_2'); zlabel('E_3');
title('Directional Young''s Modulus Surface');
view(45, 30);
light('Position',[1 1 1],'Style','infinite');
lighting gouraud;

%% Example: Comparing two Surfaces

% --- Input Data ---
% Prop1: Target Property Set (e.g. your Ground Truth)
% Prop2: Reconstructed Property Set (from forward(dlnetF, input))
% Sequence: [C11, C12, C13, C22, C23, C33, C44, C55, C66]
input1 = [0.305, 30.3, 33.4, 24.4];
prop_target = forward(dlnetF,input1);
input2 = [0.293, 34.9, 51.9, 21.2];
prop_recon = forward(dlnetF,input2);
% 1. Create a function-like setup to process both
props = {prop_target, prop_recon};
colors = {'#0072BD', '#D95319'}; % Blue for Target, Orange for Recon
names = {'Target Property', 'Reconstructed'};

figure('Color', 'w', 'Name', 'Reconstruction Loss Visualization');
hold on;

for i = 1:2
    p = props{i};
    
    % Build Stiffness Matrix C
    C = zeros(6,6);
    C(1,1)=p(1); C(1,2)=p(2); C(1,3)=p(3);
    C(2,1)=p(2); C(2,2)=p(4); C(2,3)=p(5);
    C(3,1)=p(3); C(3,2)=p(5); C(3,3)=p(6);
    C(4,4)=p(7); C(5,5)=p(8); C(6,6)=p(9);
    
    % Convert to Compliance S
    S = inv(C);
    
    % Generate Mesh
    [phi, theta] = meshgrid(linspace(0, 2*pi, 80), linspace(0, pi, 80));
    n1 = sin(theta) .* cos(phi);
    n2 = sin(theta) .* sin(phi);
    n3 = cos(theta);
    
    % Directional Young's Modulus E(n)
    invE = S(1,1)*n1.^4 + S(2,2)*n2.^4 + S(3,3)*n3.^4 ...
         + (2*S(1,2) + S(6,6))*n1.^2 .* n2.^2 ...
         + (2*S(1,3) + S(5,5))*n1.^2 .* n3.^2 ...
         + (2*S(2,3) + S(4,4))*n2.^2 .* n3.^2;
    E = 1 ./ invE;
    
    % Cartesian Coordinates
    X = E .* n1; Y = E .* n2; Z = E .* n3;
    
    % Plot with Transparency
    surf(X, Y, Z, 'FaceColor', colors{i}, 'EdgeColor', 'none', ...
         'FaceAlpha', 0.4, 'DisplayName', names{i});
end

% 2. Styling the Comparison
view(3); axis equal; grid on;
light('Position', [1 1 1]); lighting gouraud;
legend('Location', 'northeastoutside');
title('Structural Property Comparison: Target vs. Recon');
xlabel('E_x'); ylabel('E_y'); zlabel('E_z');

% Adding a simple "Loss" indicator in the corner
recon_loss = norm(prop_target - prop_recon);
annotation('textbox', [0.15, 0.8, 0.1, 0.1], 'String', ...
    sprintf('Recon Loss (L2): %.4f', recon_loss), 'FitBoxToText', 'on');