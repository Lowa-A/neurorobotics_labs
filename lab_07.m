%% Lab 07 - ERD/ERS Analysis using PSD data
% This script:
% 1. Processes three offline GDF files using PSD_files function
% 2. Loads and concatenates the processed .mat files
% 3. Extracts trials from fixation to end of continuous feedback
% 4. Computes ERD/ERS using fixation as reference
% 5. Visualizes ERD/ERS for motor imagery classes

clearvars; clc; close all;

%% Step 1: Process GDF files with PSD_files function
fprintf('========================================\n');
fprintf('Step 1: Processing GDF files\n');
fprintf('========================================\n');

gdf_files = {
    'ah7.20170613.161402.offline.mi.mi_bhbf.gdf';
    'ah7.20170613.162331.offline.mi.mi_bhbf.gdf';
    'ah7.20170613.162934.offline.mi.mi_bhbf.gdf'
};

% Process each GDF file
for i = 1:length(gdf_files)
    fprintf('\nProcessing file %d/%d: %s\n', i, length(gdf_files), gdf_files{i});
    PSD_files(gdf_files{i});
end

%% Step 2: Load and concatenate processed .mat files
fprintf('\n========================================\n');
fprintf('Step 2: Loading and concatenating data\n');
fprintf('========================================\n');

mat_files = {
    'ah7.20170613.161402.offline.mi.mi_bhbf.mat';
    'ah7.20170613.162331.offline.mi.mi_bhbf.mat';
    'ah7.20170613.162934.offline.mi.mi_bhbf.mat'
};

% Initialize variables
PSD_concat = [];
events_all = struct('TYP', [], 'POS', [], 'DUR', []);
cumulative_windows = 0;

% Load and concatenate
for i = 1:length(mat_files)
    fprintf('Loading file %d/%d: %s\n', i, length(mat_files), mat_files{i});
    data = load(mat_files{i});
    
    % Concatenate PSD along time (windows) dimension
    if isempty(PSD_concat)
        PSD_concat = data.PSD;
        f = data.f;
        samplerate = data.samplerate;
        channels = data.channels;
        wlength = data.wlength;
        wshift = data.wshift;
    else
        PSD_concat = cat(1, PSD_concat, data.PSD);
    end
    
    % Concatenate events with adjusted positions
    events_all.TYP = [events_all.TYP; data.events.TYP];
    events_all.POS = [events_all.POS; data.events.POS + cumulative_windows];
    events_all.DUR = [events_all.DUR; data.events.DUR];
    
    % Update cumulative window count
    cumulative_windows = cumulative_windows + size(data.PSD, 1);
end

fprintf('\nConcatenated PSD size: [%d windows x %d frequencies x %d channels]\n', ...
        size(PSD_concat, 1), size(PSD_concat, 2), size(PSD_concat, 3));
fprintf('Total events: %d\n', length(events_all.TYP));

%% Step 3: Extract trials (fixation to end of continuous feedback)
fprintf('\n========================================\n');
fprintf('Step 3: Extracting trials\n');
fprintf('========================================\n');

% Event type definitions
EVENT_FIXATION = 786;    % Fixation cross
EVENT_CUE_HAND = 773;    % Both hands MI
EVENT_CUE_FEET = 771;    % Both feet MI
EVENT_FEEDBACK = 781;    % Continuous feedback

% Find all fixation events (start of trial)
fix_idx = find(events_all.TYP == EVENT_FIXATION);
fprintf('Found %d fixation events (trials)\n', length(fix_idx));

% Initialize trial storage
num_trials = length(fix_idx);
Activity = [];   % Will be [windows x frequencies x channels x trials]
Reference = [];  % Will be [windows x frequencies x channels x trials]
trial_labels = zeros(num_trials, 1);  % Store cue type for each trial

% Extract each trial
valid_trials = 0;
for i = 1:num_trials
    % Get fixation event position and duration
    fix_pos = events_all.POS(fix_idx(i));
    fix_dur = events_all.DUR(fix_idx(i));
    
    % Find the continuous feedback event for this trial
    % It should be after the fixation
    fb_candidates = find(events_all.TYP == EVENT_FEEDBACK & events_all.POS > fix_pos);
    
    if isempty(fb_candidates)
        fprintf('Warning: No feedback event found for trial %d\n', i);
        continue;
    end
    
    fb_idx_local = fb_candidates(1);  % Take the first feedback after fixation
    fb_pos = events_all.POS(fb_idx_local);
    fb_dur = events_all.DUR(fb_idx_local);
    
    % Find the cue event (hand or feet) between fixation and feedback
    cue_candidates = find((events_all.TYP == EVENT_CUE_HAND | events_all.TYP == EVENT_CUE_FEET) ...
                          & events_all.POS > fix_pos & events_all.POS < fb_pos);
    
    if isempty(cue_candidates)
        fprintf('Warning: No cue event found for trial %d\n', i);
        continue;
    end
    
    cue_type = events_all.TYP(cue_candidates(1));
    
    % Calculate trial range: from fixation start to end of feedback
    trial_start = fix_pos;
    trial_end = fb_pos + fb_dur - 1;
    
    % Reference period: fixation only
    ref_start = fix_pos;
    ref_end = fix_pos + fix_dur - 1;
    
    % Check if indices are valid
    if trial_start < 1 || trial_end > size(PSD_concat, 1) || ref_end > size(PSD_concat, 1)
        fprintf('Warning: Invalid indices for trial %d\n', i);
        continue;
    end
    
    % Extract trial data
    trial_data = PSD_concat(trial_start:trial_end, :, :);  % [windows x freqs x channels]
    ref_data = PSD_concat(ref_start:ref_end, :, :);        % [windows x freqs x channels]
    
    % Store in 4D matrices
    valid_trials = valid_trials + 1;
    
    if valid_trials == 1
        % Initialize matrices with proper size
        max_trial_length = trial_end - trial_start + 1;
        max_ref_length = ref_end - ref_start + 1;
        Activity = zeros(max_trial_length, size(PSD_concat, 2), size(PSD_concat, 3), num_trials);
        Reference = zeros(max_ref_length, size(PSD_concat, 2), size(PSD_concat, 3), num_trials);
    end
    
    % Store data (pad with NaN if necessary)
    trial_length = size(trial_data, 1);
    ref_length = size(ref_data, 1);
    
    Activity(1:trial_length, :, :, valid_trials) = trial_data;
    Reference(1:ref_length, :, :, valid_trials) = ref_data;
    trial_labels(valid_trials) = cue_type;
end

% Trim to valid trials
Activity = Activity(:, :, :, 1:valid_trials);
Reference = Reference(:, :, :, 1:valid_trials);
trial_labels = trial_labels(1:valid_trials);

fprintf('Successfully extracted %d trials\n', valid_trials);
fprintf('Activity size: [%d windows x %d freqs x %d channels x %d trials]\n', size(Activity));
fprintf('Reference size: [%d windows x %d freqs x %d channels x %d trials]\n', size(Reference));

%% Step 4: Compute ERD/ERS
fprintf('\n========================================\n');
fprintf('Step 4: Computing ERD/ERS\n');
fprintf('========================================\n');

% Compute reference baseline (average across time windows for each trial)
Reference_avg = squeeze(mean(Reference, 1));  % [freqs x channels x trials]

% Compute ERD/ERS for each time window using log-ratio formula
% ERD/ERS = log(Activity / Baseline)
ERD = zeros(size(Activity));

for trial = 1:size(Activity, 4)
    for win = 1:size(Activity, 1)
        % Log-ratio ERD/ERS formula
        ERD(win, :, :, trial) = log(squeeze(Activity(win, :, :, trial)) ./ Reference_avg(:, :, trial));
    end
end

fprintf('ERD/ERS computed: [%d windows x %d freqs x %d channels x %d trials]\n', size(ERD));

%% Step 5: Select meaningful channels for motor imagery
fprintf('\n========================================\n');
fprintf('Step 5: Selecting motor cortex channels\n');
fprintf('========================================\n');

% Motor cortex channels (C3, Cz, C4 region and surrounding)
% Typical indices for 16-channel setup: 7, 9, 11 (C3, Cz, C4 region)
selected_channels = [7, 9, 11];
channel_names = {'Ch7', 'Ch9', 'Ch11'};

fprintf('Selected channels: %s\n', strjoin(channel_names, ', '));

%% Step 6: Visualize ERD/ERS averaged across trials
fprintf('\n========================================\n');
fprintf('Step 6: Visualizing ERD/ERS\n');
fprintf('========================================\n');

% Separate trials by class
hand_trials = find(trial_labels == EVENT_CUE_HAND);
feet_trials = find(trial_labels == EVENT_CUE_FEET);

fprintf('Hand MI trials: %d\n', length(hand_trials));
fprintf('Feet MI trials: %d\n', length(feet_trials));

% Average ERD/ERS across trials for each class
ERD_hand = mean(ERD(:, :, :, hand_trials), 4);  % [windows x freqs x channels]
ERD_feet = mean(ERD(:, :, :, feet_trials), 4);  % [windows x freqs x channels]

% Create time axis for windows
time_axis = (0:size(ERD, 1)-1) * wshift;  % in seconds

% Create single figure with all channels in 2x3 subplots
figure('Position', [100, 100, 1600, 800]);

% Visualize for each selected channel
for ch = 1:length(selected_channels)
    channel_idx = selected_channels(ch);
    
    % Both Hands MI - top row
    subplot(2, 3, ch);
    imagesc(time_axis, f, squeeze(ERD_hand(:, :, channel_idx))');
    axis xy;  % Flip y-axis to have low frequencies at bottom
    colorbar;
    clim([-2.5 0.7]);  % Set color scale
    colormap('hot');  % Hot colormap
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(sprintf('Both Hands MI - %s', channel_names{ch}));
    ylim([0 48]);  % Limit to meaningful frequency range
    
    % Both Feet MI - bottom row
    subplot(2, 3, ch + 3);
    imagesc(time_axis, f, squeeze(ERD_feet(:, :, channel_idx))');
    axis xy;  % Flip y-axis to have low frequencies at bottom
    colorbar;
    clim([-2.5 0.7]);  % Set color scale
    colormap('hot');  % Hot colormap
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(sprintf('Both Feet MI - %s', channel_names{ch}));
    ylim([0 48]);  % Limit to meaningful frequency range
end

sgtitle('ERD/ERS Time-Frequency Maps - All Channels');

fprintf('\nVisualization complete!\n');
