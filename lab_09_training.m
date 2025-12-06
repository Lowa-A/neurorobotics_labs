%% Lab 09 - Training Script: Feature Selection and Classifier Training
% This script implements the TRAINING phase of the BCI pipeline:
% 1. Loads calibration data (offline files 1-3)
% 2. Performs feature selection using Fisher Score
% 3. Trains QDA classifier
% 4. Evaluates training accuracy
% 5. Saves the trained model for testing

clearvars; clc; close all;

%% Configuration
fprintf('========================================\n');
fprintf('Lab 09 - Training Phase\n');
fprintf('========================================\n');

% Calibration files (offline training data)
calibration_files = {
    'ah7.20170613.161402.offline.mi.mi_bhbf.mat';
    'ah7.20170613.162331.offline.mi.mi_bhbf.mat';
    'ah7.20170613.162934.offline.mi.mi_bhbf.mat'
};

% Number of features to select
num_selected_features = 5;

% Output model filename
model_filename = 'ah7.20241209.bhbf.mat';

%% Step 1: Load calibration data
fprintf('\n========================================\n');
fprintf('Step 1: Loading calibration data\n');
fprintf('========================================\n');

PSD_all = cell(length(calibration_files), 1);
events_all = cell(length(calibration_files), 1);

for i = 1:length(calibration_files)
    fprintf('Loading file %d/%d: %s\n', i, length(calibration_files), calibration_files{i});
    data = load(calibration_files{i});
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

% Concatenate all calibration data
PSD_concat = [];
events_concat = struct('TYP', [], 'POS', [], 'DUR', []);
cumulative_windows = 0;

for i = 1:length(calibration_files)
    PSD_concat = cat(1, PSD_concat, PSD_all{i});
    events_concat.TYP = [events_concat.TYP; events_all{i}.TYP];
    events_concat.POS = [events_concat.POS; events_all{i}.POS + cumulative_windows];
    events_concat.DUR = [events_concat.DUR; events_all{i}.DUR];
    cumulative_windows = cumulative_windows + size(PSD_all{i}, 1);
end

fprintf('\nConcatenated PSD size: [%d windows x %d frequencies x %d channels]\n', ...
        size(PSD_concat, 1), size(PSD_concat, 2), size(PSD_concat, 3));

%% Step 2: Create label vectors
fprintf('\n========================================\n');
fprintf('Step 2: Creating label vectors\n');
fprintf('========================================\n');

% Event type definitions
EVENT_FIXATION = 786;
EVENT_CUE_HAND = 773;
EVENT_CUE_FEET = 771;
EVENT_REST = 783;      % Rest trials (not present in offline calibration data)
EVENT_FEEDBACK = 781;

num_windows_concat = size(PSD_concat, 1);
Ck_concat = zeros(num_windows_concat, 1);
CFbK_concat = zeros(num_windows_concat, 1);

% Track cue positions to propagate labels through feedback
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
        
        % Propagate cue label through feedback period
        cue_idx = find(cue_positions < pos, 1, 'last');
        if ~isempty(cue_idx)
            Ck_concat(pos:win_end) = cue_types(cue_idx);
        end
    end
end

fprintf('Total feedback windows: %d\n', sum(CFbK_concat == EVENT_FEEDBACK));
fprintf('Hand MI windows: %d\n', sum(Ck_concat == EVENT_CUE_HAND));
fprintf('Feet MI windows: %d\n', sum(Ck_concat == EVENT_CUE_FEET));

%% Step 3: Extract training windows
fprintf('\n========================================\n');
fprintf('Step 3: Extracting training data\n');
fprintf('========================================\n');

% Extract windows during feedback period
feedback_idx = find(CFbK_concat == EVENT_FEEDBACK);
Ck_feedback = Ck_concat(feedback_idx);

% Remove unlabeled windows
valid_idx = (Ck_feedback == EVENT_CUE_HAND | Ck_feedback == EVENT_CUE_FEET);

% Dimensions
num_freqs = size(PSD_concat, 2);
num_channels = size(PSD_concat, 3);
num_features = num_freqs * num_channels;

% Reshape PSD without log transformation (for Fisher score computation)
fprintf('Preparing raw PSD features...\n');
F_all_raw = reshape(PSD_concat, size(PSD_concat, 1), num_features);
F_train_all_raw = F_all_raw(feedback_idx, :);
F_train_all_raw = F_train_all_raw(valid_idx, :);
Ck_train = Ck_feedback(valid_idx);

fprintf('Training data: [%d windows x %d features]\n', size(F_train_all_raw));

%% Step 4: Feature Selection using Fisher Score (on RAW PSD)
fprintf('\n========================================\n');
fprintf('Step 4: Feature selection (Fisher Score)\n');
fprintf('========================================\n');

% Separate by class
class_hand = (Ck_train == EVENT_CUE_HAND);
class_feet = (Ck_train == EVENT_CUE_FEET);

F_hand = F_train_all_raw(class_hand, :);
F_feet = F_train_all_raw(class_feet, :);

% Compute Fisher Score for each feature on raw PSD
mu_hand = mean(F_hand, 1);
mu_feet = mean(F_feet, 1);
var_hand = var(F_hand, 0, 1);
var_feet = var(F_feet, 0, 1);

fisher_scores = (mu_hand - mu_feet).^2 ./ (var_hand + var_feet + eps);

% Select top features by Fisher score
[~, sorted_idx] = sort(fisher_scores, 'descend');
selected_features_idx = sorted_idx(1:num_selected_features);

% Convert to channel/frequency pairs
[selected_freqs, selected_channels] = ind2sub([num_freqs, num_channels], selected_features_idx);

fprintf('\nSelected %d features:\n', num_selected_features);
for i = 1:num_selected_features
    fprintf('  Feature %d: Channel %d, Frequency %.1f Hz (Fisher=%.4f)\n', ...
            i, selected_channels(i), f(selected_freqs(i)), fisher_scores(selected_features_idx(i)));
end

% Now apply log transformation only to selected features for training
fprintf('Applying log transformation to selected features...\n');
F_train_all_log = log(F_train_all_raw + eps);
F_train = F_train_all_log(:, selected_features_idx);

%% Step 5: Train QDA Classifier
fprintf('\n========================================\n');
fprintf('Step 5: Training QDA classifier\n');
fprintf('========================================\n');

% Convert to class labels (1=hand, 2=feet)
Ck_train_labels = zeros(size(Ck_train));
Ck_train_labels(Ck_train == EVENT_CUE_HAND) = 1;
Ck_train_labels(Ck_train == EVENT_CUE_FEET) = 2;

% Train classifier
Model = fitcdiscr(F_train, Ck_train_labels, 'DiscrimType', 'quadratic');

fprintf('QDA classifier trained successfully\n');

%% Step 6: Evaluate training accuracy
fprintf('\n========================================\n');
fprintf('Step 6: Evaluating training accuracy\n');
fprintf('========================================\n');

% Predict on training data
[Gk, pp] = predict(Model, F_train);

% Compute accuracies
overall_acc = sum(Gk == Ck_train_labels) / length(Ck_train_labels) * 100;
hand_acc = sum(Gk(Ck_train_labels == 1) == 1) / sum(Ck_train_labels == 1) * 100;
feet_acc = sum(Gk(Ck_train_labels == 2) == 2) / sum(Ck_train_labels == 2) * 100;

fprintf('\nTraining Accuracy:\n');
fprintf('  Overall: %.2f%%\n', overall_acc);
fprintf('  Hand MI: %.2f%%\n', hand_acc);
fprintf('  Feet MI: %.2f%%\n', feet_acc);

%% Step 7: Save trained model (visualization disabled)
fprintf('\n========================================\n');
fprintf('Step 7: Saving trained model\n');
fprintf('========================================\n');

% Package all necessary information for testing
fprintf('Saving model to: %s\n', model_filename);

training_accuracy = struct('overall', overall_acc, 'hand', hand_acc, 'feet', feet_acc);

save(model_filename, ...
     'Model', ...
     'selected_features_idx', ...
     'selected_channels', ...
     'selected_freqs', ...
     'f', ...
     'channels', ...
     'num_features', ...
     'num_freqs', ...
     'num_channels', ...
     'fisher_scores', ...
     'EVENT_CUE_HAND', ...
     'EVENT_CUE_FEET', ...
     'EVENT_FEEDBACK', ...
     'training_accuracy');

fprintf('\nModel saved successfully!\n');
fprintf('\nTraining complete!\n');
fprintf('========================================\n');
fprintf('Model ready for testing phase\n');
fprintf('Use this model with evaluation (online) data\n');
fprintf('========================================\n');
