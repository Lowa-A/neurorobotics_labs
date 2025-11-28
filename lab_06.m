%% Load and Concatenate Offline GDF Files
% List of offline GDF files to concatenate
offline_files = {
    'ah7.20170613.161402.offline.mi.mi_bhbf.gdf';
    'ah7.20170613.162331.offline.mi.mi_bhbf.gdf';
    'ah7.20170613.162934.offline.mi.mi_bhbf.gdf'
};

fprintf('\n=== Concatenating Offline GDF Files ===\n');

% Initialize variables for concatenation
s_concat = [];
h_concat = [];
EVENT_concat = struct('POS', [], 'TYP', [], 'DUR', []);
cumulative_samples = 0;
run_boundaries = [0];  % Store sample boundaries for each run

% a. Load each offline GDF file and concatenate
for i = 1:length(offline_files)
    fprintf('Loading file %d/%d: %s\n', i, length(offline_files), offline_files{i});
    
    % Load the GDF file
    [s_temp, h_temp] = sload(offline_files{i});
    
    % b. Concatenate the EEG data
    if isempty(s_concat)
        % First file: initialize
        s_concat = s_temp;
        h_concat = h_temp;
        EVENT_concat = h_temp.EVENT;
    else
        % Subsequent files: concatenate EEG data
        s_concat = [s_concat; s_temp];
        
        % c. Concatenate the events (POS, TYP, DUR)
        % Adjust event positions by adding cumulative sample count
        adjusted_POS = h_temp.EVENT.POS + cumulative_samples;
        
        EVENT_concat.POS = [EVENT_concat.POS; adjusted_POS];
        EVENT_concat.TYP = [EVENT_concat.TYP; h_temp.EVENT.TYP];
        EVENT_concat.DUR = [EVENT_concat.DUR; h_temp.EVENT.DUR];
    end
    
    % Update cumulative sample count for next iteration
    cumulative_samples = cumulative_samples + size(s_temp, 1);
    run_boundaries = [run_boundaries, cumulative_samples];  % Store boundary
    
    fprintf('  Samples: %d, Events: %d, Cumulative samples: %d\n', ...
            size(s_temp, 1), length(h_temp.EVENT.POS), cumulative_samples);
end

% Update the concatenated header with merged events
h_concat.EVENT = EVENT_concat;

% Use concatenated data for the rest of the script
s = s_concat;
h = h_concat;

fprintf('\nConcatenation complete!\n');
fprintf('Total samples: %d\n', size(s, 1));
fprintf('Total events: %d\n', length(h.EVENT.POS));
fprintf('===================================\n\n');

%% Extract event information from the GDF file

EVENT = h.EVENT;  % Event structure from header

% Get total number of samples in the recording
total_samples = size(s, 1);
fs = EVENT.SampleRate;
% Define event codes based on the GDF file specification
TRIAL_START = 1;          % 0x0001 - Trial start
FIXATION_CODE = 786;      % 0x0312 - Fixation cross event
CUE_BOTH_HAND = 773;      % 0x0305 - Both Hand
CUE_BOTH_FEET = 771;      % 0x0303 - Both Feet
REST_CODE = 783;          % 0x030F - Rest
FEEDBACK_CODE = 781;      % 0x030D - Continuous feedback
HIT_CODE = 897;           % 0x0381 - Target hit
MISS_CODE = 898;          % 0x0382 - Target miss
EVENT_OFF = 32768;        % 0x8000 - Event OFF

% Display unique event types in the file for verification
unique_events = unique(EVENT.TYP);
fprintf('\nUnique event codes found in GDF file:\n');
for i = 1:length(unique_events)
    fprintf('  Event code: %d (0x%04X)\n', unique_events(i), unique_events(i));
end

% a. Create Tk [trial vector] (1, 2, 3, ... N)
% Count number of trials (using trial start or cue events as trial markers)
trial_starts = EVENT.POS(EVENT.TYP == TRIAL_START);
if isempty(trial_starts)
    % If no trial start events, use cue events
    cue_positions = EVENT.POS(EVENT.TYP == CUE_BOTH_HAND | EVENT.TYP == CUE_BOTH_FEET | EVENT.TYP == REST_CODE);
else
    cue_positions = trial_starts;
end
num_trials = length(cue_positions);
Tk = zeros(1, total_samples);

% Assign trial numbers to each sample
for trial = 1:num_trials
    if trial < num_trials
        trial_start = cue_positions(trial);
        trial_end = cue_positions(trial + 1) - 1;
    else
        trial_start = cue_positions(trial);
        trial_end = total_samples;
    end
    Tk(trial_start:trial_end) = trial;
end

% Create Rk [run index] (1, 2, 3, ... N runs)
% This vector indicates which GDF file (run) each sample belongs to
Rk = zeros(1, total_samples);
num_runs = length(offline_files);

for run = 1:num_runs
    run_start = run_boundaries(run) + 1;
    run_end = run_boundaries(run + 1);
    Rk(run_start:run_end) = run;
end

% b. Create Fk [fixation periods] (0 or event value)
Fk = zeros(1, total_samples);
fixation_on_events = find(EVENT.TYP == FIXATION_CODE);

for i = 1:length(fixation_on_events)
    fix_start = EVENT.POS(fixation_on_events(i));
    fix_duration = EVENT.DUR(fixation_on_events(i));
    fix_end = fix_start + fix_duration - 1;
    Fk(fix_start:min(fix_end, total_samples)) = FIXATION_CODE;
end

% c. Create Ak [cue periods] (0 or event value)
Ak = zeros(1, total_samples);
cue_events = find(EVENT.TYP == CUE_BOTH_HAND | EVENT.TYP == CUE_BOTH_FEET | EVENT.TYP == REST_CODE);

for i = 1:length(cue_events)
    cue_start = EVENT.POS(cue_events(i));
    cue_type = EVENT.TYP(cue_events(i));
    cue_duration = EVENT.DUR(cue_events(i));
    cue_end = cue_start + cue_duration - 1;
    Ak(cue_start:min(cue_end, total_samples)) = cue_type;
end

% d. Create CFk [continuous feedback periods] (0 or event value)
CFk = zeros(1, total_samples);
feedback_on_events = find(EVENT.TYP == FEEDBACK_CODE);

for i = 1:length(feedback_on_events)
    fb_start = EVENT.POS(feedback_on_events(i));
    fb_duration = EVENT.DUR(feedback_on_events(i));
    fb_end = fb_start + fb_duration - 1;
    CFk(fb_start:min(fb_end, total_samples)) = FEEDBACK_CODE;
end

% e. Create Xk [hit/miss periods] (0 or event value)
Xk = zeros(1, total_samples);
hit_events = find(EVENT.TYP == HIT_CODE);
miss_events = find(EVENT.TYP == MISS_CODE);

% Mark hit samples (if they exist)
for i = 1:length(hit_events)
    hit_pos = EVENT.POS(hit_events(i));
    % Mark a window around the hit event (e.g., 0.5 seconds)
    window = round(0.5 * fs);
    hit_start = max(1, hit_pos - window);
    hit_end = min(total_samples, hit_pos + window);
    Xk(hit_start:hit_end) = HIT_CODE;
end

% Mark miss samples (if they exist)
for i = 1:length(miss_events)
    miss_pos = EVENT.POS(miss_events(i));
    % Mark a window around the miss event (e.g., 0.5 seconds)
    window = round(0.5 * fs);
    miss_start = max(1, miss_pos - window);
    miss_end = min(total_samples, miss_pos + window);
    Xk(miss_start:miss_end) = MISS_CODE;
end

% Note: If no hit/miss events exist, Xk remains all zeros
if isempty(hit_events) && isempty(miss_events)
    fprintf('Note: No hit/miss events found in this recording.\n');
end

% Display summary statistics
fprintf('\n=== Label Vector Summary ===\n');
fprintf('Total samples: %d\n', total_samples);
fprintf('Number of runs: %d\n', num_runs);
fprintf('Number of trials: %d\n', num_trials);
fprintf('Number of fixation events: %d\n', sum(Fk > 0));
fprintf('Number of cue samples: %d\n', sum(Ak > 0));
fprintf('Number of feedback samples: %d\n', sum(CFk > 0));
fprintf('Number of hit samples: %d\n', sum(Xk == HIT_CODE));
fprintf('Number of miss samples: %d\n', sum(Xk == MISS_CODE));
fprintf('===========================\n');

%% Spatial Filtering
fprintf('\n=== Applying Spatial Filters ===\n');

% b. Laplacian filter
fprintf('Loading Laplacian filter mask...\n');
% Load the Laplacian mask
load('laplacian16.mat');
fprintf('Laplacian mask size: [%d x %d]\n', size(lap, 1), size(lap, 2));

fprintf('Applying Laplacian filter...\n');
% Apply Laplacian spatial filter only to first 16 channels (EEG channels)
% s_lap = s_eeg * lap where s_eeg = s(:, 1:16)
s_eeg = s(:, 1:16);  % Extract only EEG channels
s = s_eeg * lap;  % [samples x 16] * [16 x 16] = [samples x 16]
fprintf('Laplacian filter applied to first 16 EEG channels.\n');
fprintf('Laplacian filtered data size: [%d x %d]\n', size(s, 1), size(s, 2));

fprintf('Spatial filtering complete!\n');
fprintf('===========================\n');

%% Signal Processing  

% Define frequency bands
mu_band = [8 12];     % μ (mu) band: 8-12 Hz
beta_band = [18 22];   % β (beta) band: 18-22 Hz

% a. Filter the signal in the μ and β bands
% Design Butterworth filters (5th order)
filter_order = 5;

% Nyquist frequency
nyquist_freq = fs / 2;

% μ band filter
[b_mu, a_mu] = butter(filter_order, mu_band / nyquist_freq, 'bandpass');
fprintf('Designed μ band filter (%.1f-%.1f Hz), order %d\n', mu_band(1), mu_band(2), filter_order);

% β band filter
[b_beta, a_beta] = butter(filter_order, beta_band / nyquist_freq, 'bandpass');
fprintf('Designed β band filter (%.1f-%.1f Hz), order %d\n', beta_band(1), beta_band(2), filter_order);

% Verify filter stability (optional - uncomment to visualize)
% fprintf('Displaying filter responses...\n');
% figure; fvtool(b_mu, a_mu); title('μ Band Filter Response');
% figure; fvtool(b_beta, a_beta); title('β Band Filter Response');

% Get number of channels
num_channels = size(s, 2);

% Initialize matrices for filtered signals
s_mu = zeros(size(s));
s_beta = zeros(size(s));

fprintf('Applying zero-phase filtering to all channels...\n');

% Apply zero-phase filtering using filtfilt for each channel
for ch = 1:num_channels
    % μ band
    s_mu(:, ch) = filtfilt(b_mu, a_mu, s(:, ch));
    
    % β band
    s_beta(:, ch) = filtfilt(b_beta, a_beta, s(:, ch));
end

fprintf('Filtering complete.\n');

% b. Rectify the signal (square it)
fprintf('Rectifying signals (squaring)...\n');
s_mu_squared = s_mu .^ 2;
s_beta_squared = s_beta .^ 2;

% c. Apply a moving average using a 1-second window
window_size = round(fs);  % 1 second window
fprintf('Applying moving average (window size: %d samples = 1 second)...\n', window_size);

% Create moving average filter
%ma_filter = ones(window_size, 1) / window_size;

% Initialize matrices for smoothed signals
s_mu_smooth = zeros(size(s_mu_squared));
s_beta_smooth = zeros(size(s_beta_squared));

% Apply moving average to each channel
for ch = 1:num_channels
    % μ band - movmean automatically handles edges
    s_mu_smooth(:, ch) = movmean(s_mu_squared(:, ch), window_size);
    
    % β band
    s_beta_smooth(:, ch) = movmean(s_beta_squared(:, ch), window_size);
end

%% Trial Extraction
fprintf('\n=== Extracting Trials ===\n');

% Extract trials corresponding to the two motor imagery (MI) tasks

% Trial lasts from fixation cross to end of continuous feedback period
% Get all fixation events (these mark trial starts)
fixation_event_indices = find(EVENT.TYP == FIXATION_CODE);
num_trials_extracted = length(fixation_event_indices);

fprintf('Found %d trials (fixation events)\n', num_trials_extracted);

% Extract fixation periods and continuous feedback periods for each trial
% FixData: fixation period [samples x channels x trials]
% TrialData: whole trial period [samples x channels x trials]

% Find maximum fixation duration and trial duration
max_fix_duration = 0;
max_trial_duration = 0;

for i = 1:num_trials_extracted
    fix_idx = fixation_event_indices(i);
    fix_start = EVENT.POS(fix_idx);
    fix_duration = EVENT.DUR(fix_idx);
    
    % Find feedback event for this trial
    fb_idx = find(EVENT.TYP == FEEDBACK_CODE & EVENT.POS > fix_start, 1, 'first');
    if ~isempty(fb_idx)
        fb_start = EVENT.POS(fb_idx);
        fb_duration = EVENT.DUR(fb_idx);
        trial_duration = (fb_start + fb_duration) - fix_start;
        
        max_fix_duration = max(max_fix_duration, fix_duration);
        max_trial_duration = max(max_trial_duration, trial_duration);
    end
end

fprintf('Maximum fixation duration: %d samples\n', max_fix_duration);
fprintf('Maximum trial duration: %d samples\n', max_trial_duration);

% Initialize arrays for fixation and trial data
FixData_mu = zeros(max_fix_duration, num_channels, num_trials_extracted);
TrialData_mu = zeros(max_trial_duration, num_channels, num_trials_extracted);
FixData_beta = zeros(max_fix_duration, num_channels, num_trials_extracted);
TrialData_beta = zeros(max_trial_duration, num_channels, num_trials_extracted);

% Store cue labels for each trial
Ck = zeros(num_trials_extracted, 1);

% Extract fixation and trial periods for each trial
valid_trials = 0;
for i = 1:num_trials_extracted
    fix_idx = fixation_event_indices(i);
    fix_start = EVENT.POS(fix_idx);
    fix_duration = EVENT.DUR(fix_idx);
    fix_end = fix_start + fix_duration - 1;
    
    % Find cue event for this trial (to identify class)
    cue_idx = find((EVENT.TYP == CUE_BOTH_HAND | EVENT.TYP == CUE_BOTH_FEET) & EVENT.POS > fix_start, 1, 'first');
    
    % Find feedback event for this trial
    fb_idx = find(EVENT.TYP == FEEDBACK_CODE & EVENT.POS > fix_start, 1, 'first');
    
    if ~isempty(fb_idx) && ~isempty(cue_idx) && fix_end <= total_samples
        fb_start = EVENT.POS(fb_idx);
        fb_duration = EVENT.DUR(fb_idx);
        fb_end = fb_start + fb_duration - 1;
        trial_end = min(fb_end, total_samples);
        
        % Extract fixation period (reference period) for μ band
        FixData_mu(1:fix_duration, :, valid_trials + 1) = s_mu_smooth(fix_start:fix_end, :);
        
        % Extract whole trial period (fixation to end of feedback) for μ band
        trial_duration = trial_end - fix_start + 1;
        TrialData_mu(1:trial_duration, :, valid_trials + 1) = s_mu_smooth(fix_start:trial_end, :);
        
        % Extract fixation period (reference period) for β band
        FixData_beta(1:fix_duration, :, valid_trials + 1) = s_beta_smooth(fix_start:fix_end, :);
        
        % Extract whole trial period (fixation to end of feedback) for β band
        TrialData_beta(1:trial_duration, :, valid_trials + 1) = s_beta_smooth(fix_start:trial_end, :);
        
        % Store cue label
        Ck(valid_trials + 1) = EVENT.TYP(cue_idx);
        
        valid_trials = valid_trials + 1;
    end
end

% Trim to valid trials
FixData_mu = FixData_mu(:, :, 1:valid_trials);
TrialData_mu = TrialData_mu(:, :, 1:valid_trials);
FixData_beta = FixData_beta(:, :, 1:valid_trials);
TrialData_beta = TrialData_beta(:, :, 1:valid_trials);
Ck = Ck(1:valid_trials);

fprintf('Successfully extracted %d valid trials\n', valid_trials);
fprintf('  Both Hand trials: %d\n', sum(Ck == CUE_BOTH_HAND));
fprintf('  Both Feet trials: %d\n', sum(Ck == CUE_BOTH_FEET));

%% Compute ERD/ERS
fprintf('\n=== Computing ERD/ERS ===\n');

% μ band ERD/ERS
% Compute reference (mean power during fixation period across all trials)
Reference_mu = repmat(mean(FixData_mu, 1), [size(TrialData_mu, 1) 1 1]);

% Compute ERD/ERS: ERD = 100 * (TrialData - Reference) ./ Reference
ERD_mu = 100 * (TrialData_mu - Reference_mu) ./ Reference_mu;

fprintf('μ band ERD/ERS computed successfully\n');
fprintf('ERD_mu size: [%d samples x %d channels x %d trials]\n', size(ERD_mu, 1), size(ERD_mu, 2), size(ERD_mu, 3));

% β band ERD/ERS
% Compute reference (mean power during fixation period across all trials)
Reference_beta = repmat(mean(FixData_beta, 1), [size(TrialData_beta, 1) 1 1]);

% Compute ERD/ERS: ERD = 100 * (TrialData - Reference) ./ Reference
ERD_beta = 100 * (TrialData_beta - Reference_beta) ./ Reference_beta;

fprintf('β band ERD/ERS computed successfully\n');
fprintf('ERD_beta size: [%d samples x %d channels x %d trials]\n', size(ERD_beta, 1), size(ERD_beta, 2), size(ERD_beta, 3));

%% Temporal Visualization
fprintf('\n=== Temporal Visualization of ERD/ERS ===\n');

% Select channel 7
selected_channel = 7;
fprintf('Selected channel: %d\n', selected_channel);

% Get trial indices for each class
hand_trials = find(Ck == CUE_BOTH_HAND);
feet_trials = find(Ck == CUE_BOTH_FEET);

% Extract ERD/ERS for selected channel and each class
% μ band
ERD_mu_hand = squeeze(ERD_mu(:, selected_channel, hand_trials));  % [samples x trials]
ERD_mu_feet = squeeze(ERD_mu(:, selected_channel, feet_trials));

% β band
ERD_beta_hand = squeeze(ERD_beta(:, selected_channel, hand_trials));
ERD_beta_feet = squeeze(ERD_beta(:, selected_channel, feet_trials));

% Compute average and standard error across trials
% μ band
avg_ERD_mu_hand = mean(ERD_mu_hand, 2);
avg_ERD_mu_feet = mean(ERD_mu_feet, 2);
se_ERD_mu_hand = std(ERD_mu_hand, 0, 2) / sqrt(size(ERD_mu_hand, 2));
se_ERD_mu_feet = std(ERD_mu_feet, 0, 2) / sqrt(size(ERD_mu_feet, 2));

% β band
avg_ERD_beta_hand = mean(ERD_beta_hand, 2);
avg_ERD_beta_feet = mean(ERD_beta_feet, 2);
se_ERD_beta_hand = std(ERD_beta_hand, 0, 2) / sqrt(size(ERD_beta_hand, 2));
se_ERD_beta_feet = std(ERD_beta_feet, 0, 2) / sqrt(size(ERD_beta_feet, 2));

% Create time vector
time_vector = (0:size(ERD_mu, 1)-1) / fs;

% Plot μ band ERD/ERS
figure('Position', [100, 100, 1400, 500]);

subplot(1, 2, 1);
hold on;
plot(time_vector, avg_ERD_mu_hand, 'r', 'LineWidth', 2, 'DisplayName', 'Both Hands');
plot(time_vector, avg_ERD_mu_feet, 'b', 'LineWidth', 2, 'DisplayName', 'Both Feet');
hold off;
xlabel('Time [s]');
ylabel('ERD/ERS [%]');
title(sprintf('\\mu Band ERD/ERS - Channel %d', selected_channel));
legend('Location', 'best');
grid on;
yline(0, 'k--', 'LineWidth', 1);

% Plot β band ERD/ERS
subplot(1, 2, 2);
hold on;
plot(time_vector, avg_ERD_beta_hand, 'r', 'LineWidth', 2, 'DisplayName', 'Both Hands');
plot(time_vector, avg_ERD_beta_feet, 'b', 'LineWidth', 2, 'DisplayName', 'Both Feet');
hold off;
xlabel('Time [s]');
ylabel('ERD/ERS [%]');
title(sprintf('\\beta Band ERD/ERS - Channel %d', selected_channel));
legend('Location', 'best');
grid on;
yline(0, 'k--', 'LineWidth', 1);

hold off;
xlabel('Time [s]');
ylabel('ERD/ERS [%]');
title(sprintf('\\beta Band ERD/ERS - Channel %d', selected_channel));
legend('Location', 'best');
grid on;
yline(0, 'k--', 'LineWidth', 1);

sgtitle('Temporal Evolution of ERD/ERS');

fprintf('Temporal visualization complete!\n');

%% Spatial Visualization - Topographic Maps
fprintf('\n=== Spatial Visualization - Topographic Maps ===\n');

% Load channel locations for topographic plotting
load('chanlocs16.mat');
fprintf('Loaded channel locations\n');

% Define time periods for spatial visualization
% Reference Period: fixation period (indices from start to mean fixation duration)
% Activity Period: continuous feedback period

% Get mean fixation duration and feedback period for indexing
mean_fix_duration = round(mean(max_fix_duration));
% Activity period: approximate continuous feedback period (e.g., samples after fixation)
activity_start_idx = mean_fix_duration + round(2 * fs);  % Start ~2s after fixation
activity_end_idx = min(activity_start_idx + round(3 * fs), size(ERD_mu, 1));  % 3s window

fprintf('Reference period: samples 1-%d\n', mean_fix_duration);
fprintf('Activity period: samples %d-%d\n', activity_start_idx, activity_end_idx);

% Compute average ERD/ERS during reference and activity periods for each class
% Class 1: Both Feet (771)
feet_trials = find(Ck == CUE_BOTH_FEET);

% Reference period (fixation) - average across time and trials
ERD_Ref_771 = mean(mean(ERD_mu(1:mean_fix_duration, :, feet_trials), 1), 3);  % [1 x channels]

% Activity period (continuous feedback) - average across time and trials  
ERD_Act_771 = mean(mean(ERD_mu(activity_start_idx:activity_end_idx, :, feet_trials), 1), 3);  % [1 x channels]

% Class 2: Both Hands (773)
hand_trials = find(Ck == CUE_BOTH_HAND);

% Reference period (fixation) - average across time and trials
ERD_Ref_773 = mean(mean(ERD_mu(1:mean_fix_duration, :, hand_trials), 1), 3);  % [1 x channels]

% Activity period (continuous feedback) - average across time and trials  
ERD_Act_773 = mean(mean(ERD_mu(activity_start_idx:activity_end_idx, :, hand_trials), 1), 3);  % [1 x channels]

% Create topographic plots for Both Feet
figure('Position', [100, 100, 1200, 500]);

% Reference period topographic map
subplot(1, 2, 1);
topoplot(squeeze(ERD_Ref_771), chanlocs16);
colorbar;
title('ERD/ERS - Reference Period (Both Feet)');
clim([-50 100]);  % Set color scale

% Activity period topographic map
subplot(1, 2, 2);
topoplot(squeeze(ERD_Act_771), chanlocs16);
colorbar;
title('ERD/ERS - Activity Period (Both Feet)');
clim([-50 100]); 

sgtitle('\mu Band ERD/ERS - Topographic Maps (Both Feet)');

% Create topographic plots for Both Hands
figure('Position', [100, 100, 1200, 500]);

% Reference period topographic map
subplot(1, 2, 1);
topoplot(squeeze(ERD_Ref_773), chanlocs16);
colorbar;
title('ERD/ERS - Reference Period (Both Hands)');
clim([-50 100]); 

% Activity period topographic map
subplot(1, 2, 2);
topoplot(squeeze(ERD_Act_773), chanlocs16);
colorbar;
title('ERD/ERS - Activity Period (Both Hands)');
clim([-50 100]); 

sgtitle('\mu Band ERD/ERS - Topographic Maps (Both Hands)');

fprintf('Spatial visualization complete!\n');
