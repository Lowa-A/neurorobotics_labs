%% Lab 08 - Feature Selection and Classification
% This script:
% 1. Loads and concatenates processed .mat files
% 2. Extracts windows from Cue to end of Continuous Feedback
% 3. Computes Fisher Score for feature selection
% 4. Trains LDA classifier and evaluates accuracy

clearvars; clc; close all;

%% Step 1: Load and concatenate processed .mat files
fprintf('========================================\n');
fprintf('Step 1: Loading and concatenating data\n');
fprintf('========================================\n');

mat_files = {
    'ah7.20170613.161402.offline.mi.mi_bhbf.mat';
    'ah7.20170613.162331.offline.mi.mi_bhbf.mat';
    'ah7.20170613.162934.offline.mi.mi_bhbf.mat'
};

% Initialize variables
PSD_all = cell(length(mat_files), 1);
events_all = cell(length(mat_files), 1);

% Load each file separately (for individual Fisher score visualization)
for i = 1:length(mat_files)
    fprintf('Loading file %d/%d: %s\n', i, length(mat_files), mat_files{i});
    data = load(mat_files{i});
    PSD_all{i} = data.PSD;
    events_all{i} = data.events;
    
    if i == 1
        f = data.f;
        samplerate = data.samplerate;
        channels = data.channels;
        wlength = data.wlength;
        wshift = data.wshift;
    end
end

% Concatenate all PSDs
PSD_concat = [];
events_concat = struct('TYP', [], 'POS', [], 'DUR', []);
cumulative_windows = 0;

for i = 1:length(mat_files)
    PSD_concat = cat(1, PSD_concat, PSD_all{i});
    events_concat.TYP = [events_concat.TYP; events_all{i}.TYP];
    events_concat.POS = [events_concat.POS; events_all{i}.POS + cumulative_windows];
    events_concat.DUR = [events_concat.DUR; events_all{i}.DUR];
    cumulative_windows = cumulative_windows + size(PSD_all{i}, 1);
end

fprintf('\nConcatenated PSD size: [%d windows x %d frequencies x %d channels]\n', ...
        size(PSD_concat, 1), size(PSD_concat, 2), size(PSD_concat, 3));
fprintf('Total events: %d\n', length(events_concat.TYP));

%% Step 2: Create label vectors for each file separately
fprintf('\n========================================\n');
fprintf('Step 2: Creating label vectors\n');
fprintf('========================================\n');

% Event type definitions
EVENT_FIXATION = 786;
EVENT_CUE_HAND = 773;
EVENT_CUE_FEET = 771;
EVENT_FEEDBACK = 781;

% Create label vectors for each file
Ck_all = cell(length(mat_files), 1);
CFbK_all = cell(length(mat_files), 1);

for file_idx = 1:length(mat_files)
    num_windows = size(PSD_all{file_idx}, 1);
    Ck = zeros(num_windows, 1);     % Cue type at each window
    CFbK = zeros(num_windows, 1);   % Feedback windows indicator
    
    % First pass: identify cue and feedback periods
    cue_positions = [];
    cue_types = [];
    
    for i = 1:length(events_all{file_idx}.TYP)
        pos = events_all{file_idx}.POS(i);
        dur = events_all{file_idx}.DUR(i);
        typ = events_all{file_idx}.TYP(i);
        
        if typ == EVENT_CUE_HAND || typ == EVENT_CUE_FEET
            cue_positions = [cue_positions; pos];
            cue_types = [cue_types; typ];
        elseif typ == EVENT_FEEDBACK
            % Mark feedback period
            win_end = min(pos + dur - 1, num_windows);
            CFbK(pos:win_end) = typ;
            
            % Find the cue that belongs to this feedback
            % (the most recent cue before this feedback)
            cue_idx = find(cue_positions < pos, 1, 'last');
            if ~isempty(cue_idx)
                % Propagate cue label through feedback period
                Ck(pos:win_end) = cue_types(cue_idx);
            end
        end
    end
    
    Ck_all{file_idx} = Ck;
    CFbK_all{file_idx} = CFbK;
    
    fprintf('File %d - Cue windows (hand): %d, (feet): %d, Feedback: %d\n', ...
            file_idx, sum(Ck == EVENT_CUE_HAND), sum(Ck == EVENT_CUE_FEET), sum(CFbK == EVENT_FEEDBACK));
end

% Also create concatenated label vectors
num_windows_concat = size(PSD_concat, 1);
Ck_concat = zeros(num_windows_concat, 1);
CFbK_concat = zeros(num_windows_concat, 1);

% First pass: identify cue and feedback periods
cue_positions = [];
cue_types = [];

for i = 1:length(events_concat.TYP)
    pos = events_concat.POS(i);
    dur = events_concat.DUR(i);
    typ = events_concat.TYP(i);
    
    if typ == EVENT_CUE_HAND || typ == EVENT_CUE_FEET
        cue_positions = [cue_positions; pos];
        cue_types = [cue_types; typ];
    elseif typ == EVENT_FEEDBACK
        win_end = min(pos + dur - 1, num_windows_concat);
        CFbK_concat(pos:win_end) = typ;
        
        % Find the cue that belongs to this feedback
        cue_idx = find(cue_positions < pos, 1, 'last');
        if ~isempty(cue_idx)
            % Propagate cue label through feedback period
            Ck_concat(pos:win_end) = cue_types(cue_idx);
        end
    end
end

%% Step 3: Extract windows from Cue to end of Continuous Feedback
fprintf('\n========================================\n');
fprintf('Step 3: Extracting task windows\n');
fprintf('========================================\n');

% Extract windows where feedback is active (uses concatenated data)
feedback_idx = find(CFbK_concat == EVENT_FEEDBACK);
P = PSD_concat(feedback_idx, :, :);  % [windows x freqs x channels]
Pk = Ck_concat(feedback_idx);        % Cue labels for these windows

% Remove any unlabeled windows
valid_idx = (Pk == EVENT_CUE_HAND | Pk == EVENT_CUE_FEET);
P = P(valid_idx, :, :);
Pk = Pk(valid_idx);

% Convert numeric labels to categorical for classification
% 773 (hand) -> 1, 771 (feet) -> 2
Pk_labels = zeros(size(Pk));
Pk_labels(Pk == EVENT_CUE_HAND) = 1;
Pk_labels(Pk == EVENT_CUE_FEET) = 2;

fprintf('Extracted %d windows for classification\n', size(P, 1));
fprintf('  Hand MI windows: %d\n', sum(Pk == EVENT_CUE_HAND));
fprintf('  Feet MI windows: %d\n', sum(Pk == EVENT_CUE_FEET));

%% Step 4: Reshape data for Fisher Score computation
fprintf('\n========================================\n');
fprintf('Step 4: Reshaping data\n');
fprintf('========================================\n');

% Reshape from [windows x freqs x channels] to [windows x features]
% where features = freqs * channels
num_windows_task = size(P, 1);
num_freqs = size(P, 2);
num_channels = size(P, 3);
num_features = num_freqs * num_channels;

P_reshaped = reshape(P, num_windows_task, num_features);

fprintf('Reshaped data: [%d windows x %d features]\n', size(P_reshaped));
fprintf('Feature dimensions: %d frequencies x %d channels = %d features\n', ...
        num_freqs, num_channels, num_features);

%% Step 5: Compute Fisher Score for each feature
fprintf('\n========================================\n');
fprintf('Step 5: Computing Fisher Scores\n');
fprintf('========================================\n');

% Apply log transformation before computing Fisher scores
fprintf('Applying log transformation to PSD...\n');
P_reshaped_log = P_reshaped;%log(P_reshaped + eps);

% Separate data by class
class1_idx = (Pk == EVENT_CUE_HAND);
class2_idx = (Pk == EVENT_CUE_FEET);

P_class1 = P_reshaped_log(class1_idx, :);
P_class2 = P_reshaped_log(class2_idx, :);

% Compute means and variances for each class
mu1 = mean(P_class1, 1);  % [1 x features]
mu2 = mean(P_class2, 1);  % [1 x features]
var1 = var(P_class1, 0, 1);  % [1 x features]
var2 = var(P_class2, 0, 1);  % [1 x features]

% Fisher Score: (mu1 - mu2)^2 / (var1 + var2)
fisher_scores = (mu1 - mu2).^2 ./ (var1 + var2 + eps);

fprintf('Fisher scores computed for %d features\n', length(fisher_scores));
fprintf('Min Fisher score: %.4f\n', min(fisher_scores));
fprintf('Max Fisher score: %.4f\n', max(fisher_scores));

%% Step 6: Compute and Visualize Fisher Score maps for each file
fprintf('\n========================================\n');
fprintf('Step 6: Computing Fisher Scores per file\n');
fprintf('========================================\n');

% Compute Fisher scores for each file separately
fisher_maps = zeros(num_freqs, num_channels, length(mat_files));

for file_idx = 1:length(mat_files)
    % Extract task windows for this file
    feedback_idx_file = find(CFbK_all{file_idx} == EVENT_FEEDBACK);
    P_file = PSD_all{file_idx}(feedback_idx_file, :, :);
    Pk_file = Ck_all{file_idx}(feedback_idx_file);
    
    % Remove unlabeled windows
    valid_idx_file = (Pk_file == EVENT_CUE_HAND | Pk_file == EVENT_CUE_FEET);
    P_file = P_file(valid_idx_file, :, :);
    Pk_file = Pk_file(valid_idx_file);
    
    % Reshape to [windows x features]
    P_file_reshaped = reshape(P_file, size(P_file, 1), num_features);
    
    % Apply log transformation
    P_file_reshaped_log = log(P_file_reshaped + eps);
    
    % Separate by class
    class1_idx_file = (Pk_file == EVENT_CUE_HAND);
    class2_idx_file = (Pk_file == EVENT_CUE_FEET);
    
    P_class1_file = P_file_reshaped_log(class1_idx_file, :);
    P_class2_file = P_file_reshaped_log(class2_idx_file, :);
    
    % Compute Fisher scores
    mu1_file = mean(P_class1_file, 1);
    mu2_file = mean(P_class2_file, 1);
    var1_file = var(P_class1_file, 0, 1);
    var2_file = var(P_class2_file, 0, 1);
    
    fisher_scores_file = (mu1_file - mu2_file).^2 ./ (var1_file + var2_file + eps);
    fisher_maps(:, :, file_idx) = reshape(fisher_scores_file, num_freqs, num_channels);
    
    fprintf('File %d: Fisher scores computed\n', file_idx);
end

% Visualize Fisher Score maps for each file (Figure 1)
figure('Position', [100, 100, 1500, 400]);

% Define channel labels matching the image
channel_labels = {'Fz', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'C3', 'C1', 'Cz', 'C2', 'C4', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4'};

for file_idx = 1:length(mat_files)
    subplot(1, 3, file_idx);
    imagesc(f, 1:num_channels, fisher_maps(:, :, file_idx)');
    colorbar;
    xlabel('Hz');
    ylabel('channel');
    title(sprintf('Calibration run %d', file_idx));
    colormap('jet');
    
    % Set y-axis to show channel labels
    set(gca, 'YTick', 1:num_channels);
    set(gca, 'YTickLabel', channel_labels);
    set(gca, 'YDir', 'reverse');  % Invert y-axis
end

sgtitle('Fisher score');


%% Step 7: Select most discriminative features
fprintf('\n========================================\n');
fprintf('Step 7: Selecting top features\n');
fprintf('========================================\n');

% Select top N features based on Fisher score
num_selected_features = 3;  % Adjust this number as needed
[~, sorted_idx] = sort(fisher_scores, 'descend');
selected_features_idx = sorted_idx(1:num_selected_features);

% Extract selected features
P_selected = P_reshaped(:, selected_features_idx);

fprintf('Selected %d most discriminative features\n', num_selected_features);

% Convert feature indices back to channel/frequency pairs
[selected_freqs, selected_channels] = ind2sub([num_freqs, num_channels], selected_features_idx);
fprintf('\nTop 10 selected features:\n');
for i = 1:min(10, num_selected_features)
    fprintf('  Feature %d: Channel %d, Frequency %.1f Hz (Fisher=%.4f)\n', ...
            i, selected_channels(i), f(selected_freqs(i)), fisher_scores(selected_features_idx(i)));
end

%% Step 8: Train QDA classifier
fprintf('\n========================================\n');
fprintf('Step 8: Training QDA classifier\n');
fprintf('========================================\n');

% Create label index for training (only continuous feedback data)
LabelIdx = CFbK_concat == EVENT_FEEDBACK;

% Extract features for training (only feedback period)
% Apply log transformation to PSD values
F = log(reshape(PSD_concat, size(PSD_concat, 1), num_features) + eps);  % Log transform
F_selected = F(:, selected_features_idx);

% Get labels for feedback period
Ck_train = Ck_concat(LabelIdx);

% Remove unlabeled windows
valid_train_idx = (Ck_train == EVENT_CUE_HAND | Ck_train == EVENT_CUE_FEET);
F_train = F_selected(LabelIdx, :);
F_train = F_train(valid_train_idx, :);
Ck_train = Ck_train(valid_train_idx);

% Convert to class labels (1 for hand, 2 for feet)
Ck_train_labels = zeros(size(Ck_train));
Ck_train_labels(Ck_train == EVENT_CUE_HAND) = 1;
Ck_train_labels(Ck_train == EVENT_CUE_FEET) = 2;

% Train Quadratic Discriminant Analysis model
qda_model = fitcdiscr(F_train, Ck_train_labels, 'DiscrimType', 'quadratic');

fprintf('QDA model trained successfully\n');
fprintf('Training data size: [%d windows x %d features]\n', size(F_train));

%% Step 9: Evaluate model (training accuracy)
fprintf('\n========================================\n');
fprintf('Step 9: Evaluating classifier\n');
fprintf('========================================\n');

% Predict on training data (continuous feedback period)
% Returns both predicted class (Gk) and posterior probabilities (pp)
[Gk, pp] = predict(qda_model, F_train);

fprintf('Prediction outputs:\n');
fprintf('  Gk (predicted classes): [%d x 1]\n', length(Gk));
fprintf('  pp (posterior probabilities): [%d x %d]\n', size(pp, 1), size(pp, 2));

% Compute single sample accuracy (percentage of correct decisions during feedback)
overall_accuracy = sum(Gk == Ck_train_labels) / length(Ck_train_labels) * 100;

% Compute per-class accuracy
hand_idx = (Ck_train_labels == 1);
feet_idx = (Ck_train_labels == 2);

hand_accuracy = sum(Gk(hand_idx) == Ck_train_labels(hand_idx)) / sum(hand_idx) * 100;
feet_accuracy = sum(Gk(feet_idx) == Ck_train_labels(feet_idx)) / sum(feet_idx) * 100;

fprintf('\nClassification Results (Single Sample Accuracy):\n');
fprintf('  Overall Accuracy: %.2f%%\n', overall_accuracy);
fprintf('  Hand MI Accuracy: %.2f%% (%d/%d correct)\n', hand_accuracy, sum(Gk(hand_idx) == Ck_train_labels(hand_idx)), sum(hand_idx));
fprintf('  Feet MI Accuracy: %.2f%% (%d/%d correct)\n', feet_accuracy, sum(Gk(feet_idx) == Ck_train_labels(feet_idx)), sum(feet_idx));
fprintf('  Total windows evaluated: %d\n', length(Gk));

%% Step 10: Visualize classifier space (Figure 2) and accuracy (Figure 3)
fprintf('\n========================================\n');
fprintf('Step 10: Visualizing results\n');
fprintf('========================================\n');

% Figure 2: Classifier space of two features (if we have at least 2 features)
if num_selected_features >= 2
    figure('Position', [100, 100, 800, 600]);
    
    % Use first two selected features for visualization
    feature1 = F_train(:, 2);
    feature2 = F_train(:, 1);
    
    % Plot both classes
    hand_mask = (Ck_train_labels == 1);
    feet_mask = (Ck_train_labels == 2);
    
    hold on;
    scatter(feature1(feet_mask), feature2(feet_mask), 30, 'ko', 'DisplayName', 'both feet');
    scatter(feature1(hand_mask), feature2(hand_mask), 30, 'bv', 'DisplayName', 'both hands');
    
    % Plot decision boundary
    x_range = linspace(min(feature1), max(feature1), 100);
    y_range = linspace(min(feature2), max(feature2), 100);
    [X_grid, Y_grid] = meshgrid(x_range, y_range);
    
    % Create grid of points for boundary
    grid_features = zeros(length(x_range) * length(y_range), num_selected_features);
    grid_features(:, 1) = X_grid(:);
    grid_features(:, 2) = Y_grid(:);
    % Fill remaining features with median values
    for i = 3:num_selected_features
        grid_features(:, i) = median(F_train(:, i));
    end
    
    % Predict on grid
    [grid_pred, ~] = predict(qda_model, grid_features);
    Z = reshape(grid_pred, length(y_range), length(x_range));
    
    % Plot boundary
    contour(X_grid, Y_grid, Z, [1.5 1.5], 'r-', 'LineWidth', 2, 'DisplayName', 'Boundary');
    
    xlabel(sprintf('C%d@%dHz', selected_channels(2), round(f(selected_freqs(2)))));
    ylabel(sprintf('C%d@%dHz', selected_channels(1), round(f(selected_freqs(1)))));
    title('Classifier space of two features');
    legend('show');
    grid on;
    hold off;
end

% Figure 3: Single sample accuracy bar plot
figure('Position', [100, 100, 800, 600]);

% Bar plot of accuracies
accuracies = [overall_accuracy, hand_accuracy, feet_accuracy];
bar_labels = {'overall', 'both hands', 'both feet'};

bar(accuracies);
set(gca, 'XTickLabel', bar_labels);
ylabel('accuracy [%]');
title('Single sample accuracy on trainset');
ylim([0 100]);
grid on;

% Add value labels on bars
for i = 1:length(accuracies)
    if accuracies(i) > 5
        text(i, accuracies(i) - 5, sprintf('%.4f', accuracies(i)), ...
             'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10);
    end
end

fprintf('\nClassification complete!\n');
