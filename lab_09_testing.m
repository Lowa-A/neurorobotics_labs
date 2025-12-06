%% Lab 09 - Testing Script: Model Evaluation with Control Framework
% This script implements the TESTING phase of the BCI pipeline:
% 1. Loads trained model from lab_09_training.m
% 2. Loads evaluation data (online files)
% 3. Extracts selected features from test data
% 4. Tests decoder (single sample accuracy)
% 5. Applies exponential smoothing control framework (trial-based accuracy)

clearvars; clc; close all;

%% Configuration
fprintf('========================================\n');
fprintf('Lab 09 - Testing Phase\n');
fprintf('========================================\n');

% Load trained model
model_filename = 'ah7.20241209.bhbf.mat';

% Evaluation files (online test data) - GDF files
evaluation_gdf_files = {
    'ah7.20170613.170929.online.mi.mi_bhbf.ema.gdf';
    'ah7.20170613.171649.online.mi.mi_bhbf.dynamic.gdf';
    'ah7.20170613.172356.online.mi.mi_bhbf.dynamic.gdf';
    'ah7.20170613.173100.online.mi.mi_bhbf.ema.gdf'
};

% Corresponding MAT files after processing
evaluation_files = {
    'ah7.20170613.170929.online.mi.mi_bhbf.ema.mat';
    'ah7.20170613.171649.online.mi.mi_bhbf.dynamic.mat';
    'ah7.20170613.172356.online.mi.mi_bhbf.dynamic.mat';
    'ah7.20170613.173100.online.mi.mi_bhbf.ema.mat'
};

% Control framework parameters
alpha = 0.95;                   % Integration parameter for exponential smoothing
threshold_class1 = 0.3;         % Threshold for class 1 (hand) - D(t) <= threshold
threshold_class2 = 0.7;         % Threshold for class 2 (feet) - D(t) >= threshold
use_rejection = true;           % Enable rejection (don't force decision if no threshold crossed)
decision_strategy = 'threshold';  % 'threshold' or 'end-of-trial'

fprintf('\nControl Framework Parameters:\n');
fprintf('  Alpha (integration): %.2f\n', alpha);
fprintf('  Threshold Class 1 (Hand): D(t) <= %.2f\n', threshold_class1);
fprintf('  Threshold Class 2 (Feet): D(t) >= %.2f\n', threshold_class2);
fprintf('  Rejection: %s\n', mat2str(use_rejection));
fprintf('  Decision Strategy: %s\n', decision_strategy);

%% Step 0: Pre-process online GDF files (if needed)
fprintf('\n========================================\n');
fprintf('Step 0: Pre-processing evaluation data\n');
fprintf('========================================\n');

% Check if MAT files exist, if not, process GDF files
files_to_process = [];
for i = 1:length(evaluation_files)
    if ~exist(evaluation_files{i}, 'file')
        files_to_process = [files_to_process; i];
        fprintf('MAT file not found: %s\n', evaluation_files{i});
    end
end

if ~isempty(files_to_process)
    fprintf('\nProcessing %d GDF files using PSD_files.m...\n', length(files_to_process));
    for i = files_to_process'
        fprintf('\nProcessing file %d/%d: %s\n', i, length(evaluation_gdf_files), evaluation_gdf_files{i});
        PSD_files(evaluation_gdf_files{i});
    end
    fprintf('\nAll evaluation files processed successfully!\n');
else
    fprintf('All MAT files already exist. Skipping pre-processing.\n');
end

%% Step 1: Load trained model
fprintf('\n========================================\n');
fprintf('Step 1: Loading trained model\n');
fprintf('========================================\n');

if ~exist(model_filename, 'file')
    error('Trained model not found: %s\nPlease run lab_09_training.m first!', model_filename);
end

fprintf('Loading model from: %s\n', model_filename);
model_data = load(model_filename);

% Extract model components
Model = model_data.Model;
selected_features_idx = model_data.selected_features_idx;
selected_channels = model_data.selected_channels;
selected_freqs = model_data.selected_freqs;
f = model_data.f;
num_features = model_data.num_features;
num_freqs = model_data.num_freqs;
num_channels = model_data.num_channels;
EVENT_CUE_HAND = model_data.EVENT_CUE_HAND;
EVENT_CUE_FEET = model_data.EVENT_CUE_FEET;
EVENT_FEEDBACK = model_data.EVENT_FEEDBACK;

fprintf('\nModel loaded successfully!\n');
fprintf('Selected features: %d\n', length(selected_features_idx));
for i = 1:length(selected_features_idx)
    fprintf('  Feature %d: Channel %d, Frequency %.1f Hz\n', ...
            i, selected_channels(i), f(selected_freqs(i)));
end

%% Step 2: Load evaluation data
fprintf('\n========================================\n');
fprintf('Step 2: Loading evaluation data\n');
fprintf('========================================\n');

PSD_eval_all = cell(length(evaluation_files), 1);
events_eval_all = cell(length(evaluation_files), 1);

for i = 1:length(evaluation_files)
    fprintf('Loading file %d/%d: %s\n', i, length(evaluation_files), evaluation_files{i});
    data = load(evaluation_files{i});
    PSD_eval_all{i} = data.PSD;
    events_eval_all{i} = data.events;
end

% Concatenate all evaluation data
PSD_eval_concat = [];
events_eval_concat = struct('TYP', [], 'POS', [], 'DUR', []);
cumulative_windows = 0;

for i = 1:length(evaluation_files)
    PSD_eval_concat = cat(1, PSD_eval_concat, PSD_eval_all{i});
    events_eval_concat.TYP = [events_eval_concat.TYP; events_eval_all{i}.TYP];
    events_eval_concat.POS = [events_eval_concat.POS; events_eval_all{i}.POS + cumulative_windows];
    events_eval_concat.DUR = [events_eval_concat.DUR; events_eval_all{i}.DUR];
    cumulative_windows = cumulative_windows + size(PSD_eval_all{i}, 1);
end

fprintf('\nConcatenated evaluation PSD size: [%d windows x %d frequencies x %d channels]\n', ...
        size(PSD_eval_concat, 1), size(PSD_eval_concat, 2), size(PSD_eval_concat, 3));

%% Step 3: Create label vectors for evaluation data
fprintf('\n========================================\n');
fprintf('Step 3: Creating label vectors\n');
fprintf('========================================\n');

num_windows_eval = size(PSD_eval_concat, 1);
Ck_eval = zeros(num_windows_eval, 1);
CFbK_eval = zeros(num_windows_eval, 1);

% Track cue positions to propagate labels through feedback
% We track ALL cues (hand, feet, and rest) to properly label feedback periods
EVENT_REST = 783;
cue_positions = [];
cue_types = [];

for i = 1:length(events_eval_concat.TYP)
    pos = events_eval_concat.POS(i);
    dur = events_eval_concat.DUR(i);
    typ = events_eval_concat.TYP(i);
    
    if typ == EVENT_CUE_HAND || typ == EVENT_CUE_FEET || typ == EVENT_REST
        cue_positions = [cue_positions; pos];
        cue_types = [cue_types; typ];
    elseif typ == EVENT_FEEDBACK
        win_end = min(pos + dur - 1, num_windows_eval);
        CFbK_eval(pos:win_end) = typ;
        
        % Propagate cue label through feedback period
        cue_idx = find(cue_positions < pos, 1, 'last');
        if ~isempty(cue_idx)
            Ck_eval(pos:win_end) = cue_types(cue_idx);
        end
    end
end

fprintf('Total feedback windows: %d\n', sum(CFbK_eval == EVENT_FEEDBACK));
fprintf('Hand MI windows: %d\n', sum(Ck_eval == EVENT_CUE_HAND));
fprintf('Feet MI windows: %d\n', sum(Ck_eval == EVENT_CUE_FEET));

%% Step 4: Extract selected features from test data
fprintf('\n========================================\n');
fprintf('Step 4: Extracting test features\n');
fprintf('========================================\n');

% Reshape RAW PSD first
fprintf('Reshaping test PSD...\n');
P_eval_all = reshape(PSD_eval_concat, size(PSD_eval_concat, 1), num_features);

% Apply log transformation
fprintf('Applying log transformation to test PSD...\n');
F_eval_all = log(P_eval_all + eps);

% Extract feedback period windows
feedback_idx = find(CFbK_eval == EVENT_FEEDBACK);
Ck_test = Ck_eval(feedback_idx);

% Remove unlabeled windows
valid_idx = (Ck_test == EVENT_CUE_HAND | Ck_test == EVENT_CUE_FEET);

F_test_all = F_eval_all(feedback_idx, :);
F_test_all = F_test_all(valid_idx, :);
Ck_test = Ck_test(valid_idx);

% Extract ONLY the selected features (same as training)
F_test = F_test_all(:, selected_features_idx);

% Convert to class labels
Ck_test_labels = zeros(size(Ck_test));
Ck_test_labels(Ck_test == EVENT_CUE_HAND) = 1;
Ck_test_labels(Ck_test == EVENT_CUE_FEET) = 2;

fprintf('Test data: [%d windows x %d features]\n', size(F_test));

%% Step 5: Test decoder (Single Sample Accuracy)
fprintf('\n========================================\n');
fprintf('Step 5: Testing decoder (Single Sample)\n');
fprintf('========================================\n');

% Predict on test data
[Gk, pp] = predict(Model, F_test);

% Compute single sample accuracy
overall_acc = sum(Gk == Ck_test_labels) / length(Ck_test_labels) * 100;
hand_acc = sum(Gk(Ck_test_labels == 1) == 1) / sum(Ck_test_labels == 1) * 100;
feet_acc = sum(Gk(Ck_test_labels == 2) == 2) / sum(Ck_test_labels == 2) * 100;

fprintf('\nSingle Sample Accuracy (Test Data):\n');
fprintf('  Overall: %.2f%%\n', overall_acc);
fprintf('  Hand MI: %.2f%%\n', hand_acc);
fprintf('  Feet MI: %.2f%%\n', feet_acc);

%% Step 6: Apply Control Framework (Exponential Smoothing)
fprintf('\n========================================\n');
fprintf('Step 6: Applying control framework\n');
fprintf('========================================\n');

% Find all feedback trials
feedback_events = find(events_eval_concat.TYP == EVENT_FEEDBACK);
num_trials = length(feedback_events);

fprintf('Total feedback events: %d (includes rest trials)\n', num_trials);

% Initialize storage for trial-based results
trial_decisions = zeros(num_trials, 1);      % Final decision per trial
trial_true_labels = zeros(num_trials, 1);    % True labels per trial
trial_decision_times = zeros(num_trials, 1); % Time to decision (in windows)
D_all = cell(num_trials, 1);                 % Store D(t) trajectory for each trial

% Process each trial
for trial_idx = 1:num_trials
    % Get trial window range
    trial_start = events_eval_concat.POS(feedback_events(trial_idx));
    trial_dur = events_eval_concat.DUR(feedback_events(trial_idx));
    trial_end = min(trial_start + trial_dur - 1, num_windows_eval);
    
    % Get true label for this trial
    trial_label_idx = find(cue_positions < trial_start, 1, 'last');
    if ~isempty(trial_label_idx)
        trial_true_label = cue_types(trial_label_idx);
        if trial_true_label == EVENT_CUE_HAND
            trial_true_labels(trial_idx) = 1;
        elseif trial_true_label == EVENT_CUE_FEET
            trial_true_labels(trial_idx) = 2;
        % else: rest trial, leave as 0 (will be filtered out)
        end
    end
    
    % Get features for this trial
    trial_windows = trial_start:trial_end;
    F_trial = F_eval_all(trial_windows, selected_features_idx);
    
    % Predict posterior probabilities for trial
    [~, pp_trial] = predict(Model, F_trial);
    
    % Apply exponential smoothing - D is [nwindows x 2] for both classes
    D = zeros(length(trial_windows), 2);
    D(1, :) = [0.5, 0.5];  % Reset at beginning of trial
    
    decision_made = false;
    decision_time = length(trial_windows);
    
    for t = 2:length(trial_windows)
        % Exponential smoothing: D(t) = D(t-1) * alpha + pp(t) * (1-alpha)
        D(t, :) = D(t-1, :) * alpha + pp_trial(t, :) * (1 - alpha);
        
        % Use probability/control of class 2 (feet) as decision variable
        prob_feet = D(t, 2);
        
        % THRESHOLD STRATEGY: Check if threshold reached during trial
        if strcmp(decision_strategy, 'threshold') && ~decision_made
            if prob_feet <= threshold_class1  % Low probability of feet
                trial_decisions(trial_idx) = 1;  % Hand
                decision_made = true;
                decision_time = t;
            elseif prob_feet >= threshold_class2  % High probability of feet
                trial_decisions(trial_idx) = 2;  % Feet
                decision_made = true;
                decision_time = t;
            end
        end
    end
    
    % Make decision based on strategy
    if strcmp(decision_strategy, 'end-of-trial')
        % END-OF-TRIAL STRATEGY: Use final integrated probability
        final_prob = D(end, 2);  % Use feet probability
        
        if final_prob <= threshold_class1
            trial_decisions(trial_idx) = 1;  % Hand
            decision_made = true;
        elseif final_prob >= threshold_class2
            trial_decisions(trial_idx) = 2;  % Feet
            decision_made = true;
        else
            % Final probability in uncertain region
            if use_rejection
                trial_decisions(trial_idx) = 0;  % Rejected
            else
                % Force decision based on which class is closer
                if final_prob < 0.5
                    trial_decisions(trial_idx) = 1;  % Hand
                else
                    trial_decisions(trial_idx) = 2;  % Feet
                end
                decision_made = true;
            end
        end
        decision_time = length(trial_windows);
    else
        % THRESHOLD STRATEGY: Handle case where no threshold crossed
        if ~decision_made
            if use_rejection
                trial_decisions(trial_idx) = 0;  % Rejected (no decision)
            else
                % Force decision based on final probability
                if D(end, 1) >= 0.5
                    trial_decisions(trial_idx) = 1;  % Hand
                else
                    trial_decisions(trial_idx) = 2;  % Feet
                end
            end
            decision_time = length(trial_windows);
        end
    end
    
    trial_decision_times(trial_idx) = decision_time;
    D_all{trial_idx} = D(:, 2);  % Store class 2 (feet) probability for plotting
end

% Compute trial-based accuracy
valid_trials = (trial_true_labels > 0);  % Exclude rest trials
trial_decisions_valid = trial_decisions(valid_trials);
trial_true_labels_valid = trial_true_labels(valid_trials);
num_valid_trials = length(trial_true_labels_valid);

fprintf('Motor imagery trials evaluated: %d (excluding %d rest trials)\n', num_valid_trials, num_trials - num_valid_trials);

% Identify rejected trials (decision = 0)
rejected_trials = (trial_decisions_valid == 0);
num_rejected = sum(rejected_trials);
rejection_rate = num_rejected / num_valid_trials * 100;

% Compute accuracy WITHOUT rejection (force all trials to have decisions)
trial_decisions_no_reject = trial_decisions_valid;
for i = 1:length(trial_decisions_no_reject)
    if trial_decisions_no_reject(i) == 0
        % Use final probability for rejected trials
        trial_idx = find(valid_trials, i, 'first');
        trial_idx = trial_idx(end);
        if D_all{trial_idx}(end) >= 0.5
            trial_decisions_no_reject(i) = 1;
        else
            trial_decisions_no_reject(i) = 2;
        end
    end
end

overall_acc_no_reject = sum(trial_decisions_no_reject == trial_true_labels_valid) / num_valid_trials * 100;
hand_acc_no_reject = sum(trial_decisions_no_reject(trial_true_labels_valid == 1) == 1) / sum(trial_true_labels_valid == 1) * 100;
feet_acc_no_reject = sum(trial_decisions_no_reject(trial_true_labels_valid == 2) == 2) / sum(trial_true_labels_valid == 2) * 100;

fprintf('\nTrial-Based Accuracy (WITHOUT Rejection):\n');
fprintf('  Overall: %.2f%% (%d/%d correct)\n', overall_acc_no_reject, sum(trial_decisions_no_reject == trial_true_labels_valid), num_valid_trials);
fprintf('  Hand MI: %.2f%% (%d/%d correct)\n', hand_acc_no_reject, sum(trial_decisions_no_reject(trial_true_labels_valid == 1) == 1), sum(trial_true_labels_valid == 1));
fprintf('  Feet MI: %.2f%% (%d/%d correct)\n', feet_acc_no_reject, sum(trial_decisions_no_reject(trial_true_labels_valid == 2) == 2), sum(trial_true_labels_valid == 2));

% Compute accuracy WITH rejection (only on decided trials)
if num_rejected > 0
    decided_trials = ~rejected_trials;
    trial_decisions_decided = trial_decisions_valid(decided_trials);
    trial_true_labels_decided = trial_true_labels_valid(decided_trials);
    
    overall_acc_with_reject = sum(trial_decisions_decided == trial_true_labels_decided) / length(trial_decisions_decided) * 100;
    hand_decided = (trial_true_labels_decided == 1);
    feet_decided = (trial_true_labels_decided == 2);
    hand_acc_with_reject = sum(trial_decisions_decided(hand_decided) == 1) / sum(hand_decided) * 100;
    feet_acc_with_reject = sum(trial_decisions_decided(feet_decided) == 2) / sum(feet_decided) * 100;
    
    fprintf('\nTrial-Based Accuracy (WITH Rejection):\n');
    fprintf('  Rejected trials: %d/%d (%.2f%%)\n', num_rejected, num_valid_trials, rejection_rate);
    fprintf('  Overall: %.2f%% (%d/%d correct)\n', overall_acc_with_reject, sum(trial_decisions_decided == trial_true_labels_decided), length(trial_decisions_decided));
    fprintf('  Hand MI: %.2f%% (%d/%d correct)\n', hand_acc_with_reject, sum(trial_decisions_decided(hand_decided) == 1), sum(hand_decided));
    fprintf('  Feet MI: %.2f%% (%d/%d correct)\n', feet_acc_with_reject, sum(trial_decisions_decided(feet_decided) == 2), sum(feet_decided));
else
    fprintf('\nNo trials rejected (all trials reached decision thresholds)\n');
    overall_acc_with_reject = overall_acc_no_reject;
end

% Average decision time (only for decided trials)
decided_idx = find(valid_trials);
decided_idx = decided_idx(~rejected_trials);
if ~isempty(decided_idx)
    avg_decision_time = mean(trial_decision_times(decided_idx));
    fprintf('\nAverage decision time (decided trials): %.1f windows (%.2f seconds)\n', avg_decision_time, avg_decision_time * 0.0625);
end

%% Step 7: Visualize Results
fprintf('\n========================================\n');
fprintf('Step 7: Visualizing results\n');
fprintf('========================================\n');

% Figure 1: Single sample accuracy on test set
figure('Position', [100, 100, 600, 500]);
accuracies_single = [overall_acc, hand_acc, feet_acc];
bar(accuracies_single, 'FaceColor', [0.2 0.5 0.8]);
set(gca, 'XTickLabel', {'overall', 'both hands', 'both feet'});
ylabel('accuracy [%]');
title('Single sample accuracy on test set');
ylim([0 100]);
grid on;

% Figure 2+: Evidence accumulation for last 30 trials (separate figures)
% Determine which trials to plot (last 30 valid trials)
% valid_trial_indices = find(trial_true_labels > 0);
% num_valid_trials = length(valid_trial_indices);
% trials_to_plot = valid_trial_indices(max(1, num_valid_trials-29):num_valid_trials);
% num_plots = length(trials_to_plot);
% 
% for plot_idx = 1:num_plots
%     trial_idx_plot = trials_to_plot(plot_idx);
% 
%     figure('Position', [100 + mod(plot_idx-1, 10)*50, 100 + floor((plot_idx-1)/10)*50, 800, 600]);
% 
%     % Get raw probabilities for this trial
%     trial_start = events_eval_concat.POS(feedback_events(trial_idx_plot));
%     trial_dur = events_eval_concat.DUR(feedback_events(trial_idx_plot));
%     trial_end = min(trial_start + trial_dur - 1, num_windows_eval);
%     trial_windows = trial_start:trial_end;
%     F_trial = F_eval_all(trial_windows, selected_features_idx);
%     [~, pp_trial] = predict(Model, F_trial);
% 
%     % Plot raw probabilities as circles
%     plot(1:length(pp_trial), pp_trial(:, 1), 'o', 'Color', [0.7 0.7 0.7], 'MarkerSize', 5);
%     hold on;
% 
%     % Plot integrated probabilities
%     plot(1:length(D_all{trial_idx_plot}), D_all{trial_idx_plot}, 'b-', 'LineWidth', 2);
% 
%     % Thresholds
%     yline(threshold_class1, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
%     yline(threshold_class2, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
%     yline(0.5, ':', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
% 
%     % Mark decision point if made before end
%     if trial_decision_times(trial_idx_plot) < length(D_all{trial_idx_plot})
%         xline(trial_decision_times(trial_idx_plot), 'r--', 'LineWidth', 1.5);
%     end
% 
%     % Add threshold labels
%     text(length(pp_trial)*1.02, threshold_class1, 'Th_1', 'FontSize', 10);
%     text(length(pp_trial)*1.02, threshold_class2, 'Th_2', 'FontSize', 10);
% 
%     % Title with trial info
%     true_class = trial_true_labels(trial_idx_plot);
%     pred_class = trial_decisions(trial_idx_plot);
% 
%     if pred_class == 0
%         marker = 'REJECTED';
%         title_color = [0.8 0.8 0.0];  % Yellow for rejected
%     elseif true_class == pred_class
%         marker = '✓';
%         title_color = 'green';
%     else
%         marker = '✗';
%         title_color = 'red';
%     end
% 
%     class_name = 'both hands';
%     if true_class == 2
%         class_name = 'both feet';
%     end
% 
%     title(sprintf('Trial %d - Class %s %s', trial_idx_plot, class_name, marker), ...
%           'FontSize', 12, 'Color', title_color);
%     xlabel('Window');
%     ylabel('P(class)');
%     ylim([0 1]);
%     xlim([0 length(pp_trial)+5]);
%     grid on;
% end

% Figure 3: Trial accuracy comparison
figure('Position', [100, 100, 600, 500]);
if num_rejected > 0
    accuracies_trial_plot = [overall_acc_no_reject, overall_acc_with_reject];
else
    accuracies_trial_plot = [overall_acc_no_reject, overall_acc_no_reject];
end
bar(accuracies_trial_plot, 'FaceColor', [0.2 0.5 0.8]);
set(gca, 'XTickLabel', {'without rejection', 'with rejection'});
ylabel('accuracy [%]');
title('Trial accuracy on test');
ylim([0 100]);
grid on;

fprintf('\nTesting complete!\n');
fprintf('========================================\n');
fprintf('Summary:\n');
fprintf('  Single Sample Accuracy: %.2f%%\n', overall_acc);
fprintf('  Trial-Based Accuracy (no rejection): %.2f%%\n', overall_acc_no_reject);
if exist('overall_acc_with_reject', 'var') && num_rejected > 0
    fprintf('  Trial-Based Accuracy (with rejection): %.2f%%\n', overall_acc_with_reject);
    fprintf('  Rejection Rate: %.2f%%\n', rejection_rate);
    fprintf('  Improvement: %.2f%%\n', overall_acc_with_reject - overall_acc);
else
    fprintf('  Improvement: %.2f%%\n', overall_acc_no_reject - overall_acc);
end
fprintf('========================================\n');
