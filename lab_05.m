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
        
        % Handle DUR field if it exists
        if isfield(h_temp.EVENT, 'DUR')
            if isfield(EVENT_concat, 'DUR')
                EVENT_concat.DUR = [EVENT_concat.DUR; h_temp.EVENT.DUR];
            else
                EVENT_concat.DUR = h_temp.EVENT.DUR;
            end
        end
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
fixation_off_code = FIXATION_CODE + EVENT_OFF;  % OFF event = ON + 32768

for i = 1:length(fixation_on_events)
    fix_start = EVENT.POS(fixation_on_events(i));
    % Find the corresponding OFF event
    fix_off_idx = find(EVENT.POS > fix_start & EVENT.TYP == fixation_off_code, 1, 'first');
    if ~isempty(fix_off_idx)
        fix_end = EVENT.POS(fix_off_idx) - 1;
    else
        fix_end = min(fix_start + round(2*fs), total_samples);
    end
    Fk(fix_start:min(fix_end, total_samples)) = FIXATION_CODE;
end

% c. Create Ak [cue periods] (0 or event value)
Ak = zeros(1, total_samples);
cue_events = find(EVENT.TYP == CUE_BOTH_HAND | EVENT.TYP == CUE_BOTH_FEET | EVENT.TYP == REST_CODE);

for i = 1:length(cue_events)
    cue_start = EVENT.POS(cue_events(i));
    cue_type = EVENT.TYP(cue_events(i));
    cue_off_code = cue_type + EVENT_OFF;  % OFF event = ON + 32768
    
    % Find the corresponding OFF event
    cue_off_idx = find(EVENT.POS > cue_start & EVENT.TYP == cue_off_code, 1, 'first');
    if ~isempty(cue_off_idx)
        cue_end = EVENT.POS(cue_off_idx) - 1;
    else
        cue_end = min(cue_start + round(3*fs), total_samples);
    end
    Ak(cue_start:min(cue_end, total_samples)) = cue_type;
end

% d. Create CFk [continuous feedback periods] (0 or event value)
CFk = zeros(1, total_samples);
feedback_on_events = find(EVENT.TYP == FEEDBACK_CODE);
feedback_off_code = FEEDBACK_CODE + EVENT_OFF;  % OFF event = ON + 32768

for i = 1:length(feedback_on_events)
    fb_start = EVENT.POS(feedback_on_events(i));
    % Find the corresponding OFF event
    fb_off_idx = find(EVENT.POS > fb_start & EVENT.TYP == feedback_off_code, 1, 'first');
    if ~isempty(fb_off_idx)
        fb_end = EVENT.POS(fb_off_idx) - 1;
    else
        fb_end = min(fb_start + round(4*fs), total_samples);
    end
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

% a. Common Average Reference (CAR) filter
fprintf('Applying CAR filter...\n');
% CAR: subtract the average of all channels from each channel
% Apply CAR only to first 16 channels (EEG channels)
s_car = s(:, 1:16) - mean(s(:, 1:16), 2);
fprintf('CAR filter applied to first 16 EEG channels.\n');

% b. Laplacian filter
fprintf('Loading Laplacian filter mask...\n');
% Load the Laplacian mask
load('laplacian16.mat');
fprintf('Laplacian mask size: [%d x %d]\n', size(lap, 1), size(lap, 2));

fprintf('Applying Laplacian filter...\n');
% Apply Laplacian spatial filter only to first 16 channels (EEG channels)
% s_lap = s_eeg * lap where s_eeg = s(:, 1:16)
s_eeg = s(:, 1:16);  % Extract only EEG channels
s_lap = s_eeg * lap;  % [samples x 16] * [16 x 16] = [samples x 16]
fprintf('Laplacian filter applied to first 16 EEG channels.\n');
fprintf('Laplacian filtered data size: [%d x %d]\n', size(s_lap, 1), size(s_lap, 2));

fprintf('Spatial filtering complete!\n');
fprintf('  Original data: [%d x %d]\n', size(s, 1), size(s, 2));
fprintf('  Original EEG only: [%d x %d]\n', size(s_eeg, 1), size(s_eeg, 2));
fprintf('  CAR filtered: [%d x %d]\n', size(s_car, 1), size(s_car, 2));
fprintf('  Laplacian filtered: [%d x %d]\n', size(s_lap, 1), size(s_lap, 2));
fprintf('===========================\n');

%% Signal Processing - Logarithmic Band Power (No Spatial Filter)
fprintf('\n=== Computing Logarithmic Band Power (No Spatial Filter) ===\n'); 

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
ma_filter = ones(window_size, 1) / window_size;

% Initialize matrices for smoothed signals
s_mu_smooth = zeros(size(s_mu_squared));
s_beta_smooth = zeros(size(s_beta_squared));

% Apply moving average to each channel
for ch = 1:num_channels
    % μ band - use 'same' to maintain length
    s_mu_smooth(:, ch) = conv(s_mu_squared(:, ch), ma_filter, 'same');
    
    % β band
    s_beta_smooth(:, ch) = conv(s_beta_squared(:, ch), ma_filter, 'same');
end

% d. Apply a logarithmic transform to the result
fprintf('Applying logarithmic transform...\n');

% Add small epsilon to avoid log(0)
epsilon = 1e-10;
s_mu_logbp = log(s_mu_smooth + epsilon);
s_beta_logbp = log(s_beta_smooth + epsilon);

fprintf('Logarithmic band power computation complete!\n');
fprintf('  μ band log BP size: [%d x %d]\n', size(s_mu_logbp, 1), size(s_mu_logbp, 2));
fprintf('  β band log BP size: [%d x %d]\n', size(s_beta_logbp, 1), size(s_beta_logbp, 2));
fprintf('===========================\n');

% % Visualize the processing pipeline for one channel
% sample_channel = 7;  % C3 channel (or modify as needed)
% sample_start = round(60 * fs);  % Start at 60 seconds
% sample_duration = round(10 * fs);  % Show 10 seconds
% sample_end = min(sample_start + sample_duration, total_samples);
% sample_indices = sample_start:sample_end;
% time_samples = (sample_indices - 1) / fs;
% 
% figure('Position', [100, 100, 1400, 900]);
% 
% % Original signal
% subplot(5, 2, 1);
% plot(time_samples, s(sample_indices, sample_channel), 'k', 'LineWidth', 1);
% ylabel('Amplitude');
% title(sprintf('Original Signal (Channel %d)', sample_channel));
% grid on;
% 
% subplot(5, 2, 2);
% plot(time_samples, s(sample_indices, sample_channel), 'k', 'LineWidth', 1);
% ylabel('Amplitude');
% title(sprintf('Original Signal (Channel %d)', sample_channel));
% grid on;
% 
% % Filtered signals
% subplot(5, 2, 3);
% plot(time_samples, s_mu(sample_indices, sample_channel), 'b', 'LineWidth', 1);
% ylabel('Amplitude');
% title(sprintf('μ Band Filtered (%.0f-%.0f Hz)', mu_band(1), mu_band(2)));
% grid on;
% 
% subplot(5, 2, 4);
% plot(time_samples, s_beta(sample_indices, sample_channel), 'r', 'LineWidth', 1);
% ylabel('Amplitude');
% title(sprintf('β Band Filtered (%.0f-%.0f Hz)', beta_band(1), beta_band(2)));
% grid on;
% 
% % Squared signals
% subplot(5, 2, 5);
% plot(time_samples, s_mu_squared(sample_indices, sample_channel), 'b', 'LineWidth', 1);
% ylabel('Power');
% title('μ Band Squared');
% grid on;
% 
% subplot(5, 2, 6);
% plot(time_samples, s_beta_squared(sample_indices, sample_channel), 'r', 'LineWidth', 1);
% ylabel('Power');
% title('β Band Squared');
% grid on;
% 
% % Smoothed signals
% subplot(5, 2, 7);
% plot(time_samples, s_mu_smooth(sample_indices, sample_channel), 'b', 'LineWidth', 1);
% ylabel('Power');
% title('μ Band Smoothed (1s MA)');
% grid on;
% 
% subplot(5, 2, 8);
% plot(time_samples, s_beta_smooth(sample_indices, sample_channel), 'r', 'LineWidth', 1);
% ylabel('Power');
% title('β Band Smoothed (1s MA)');
% grid on;
% 
% % Log BP
% subplot(5, 2, 9);
% plot(time_samples, s_mu_logbp(sample_indices, sample_channel), 'b', 'LineWidth', 1.5);
% ylabel('Log Power');
% title('μ Band Log BP');
% xlabel('Time (s)');
% grid on;
% 
% subplot(5, 2, 10);
% plot(time_samples, s_beta_logbp(sample_indices, sample_channel), 'r', 'LineWidth', 1.5);
% ylabel('Log Power');
% title('β Band Log BP');
% xlabel('Time (s)');
% grid on;
% 
% sgtitle(sprintf('Signal Processing Pipeline - Channel %d (10s sample)', sample_channel));

%% Trial Extraction
fprintf('\n=== Extracting Trials ===\n');

% Trial lasts from fixation cross to end of continuous feedback period
% Get all fixation events (these mark trial starts)
fixation_event_indices = find(EVENT.TYP == FIXATION_CODE);
num_trials_extracted = length(fixation_event_indices);

fprintf('Found %d trials (fixation events)\n', num_trials_extracted);

% Determine trial length by analyzing the structure
% Trial structure: Fixation -> Cue -> Feedback
trial_durations = [];
for i = 1:min(10, num_trials_extracted)  % Check first few trials
    fix_pos = EVENT.POS(fixation_event_indices(i));
    
    % Find end of feedback for this trial
    fb_off_events_after_fix = find(EVENT.POS > fix_pos & EVENT.TYP == (FEEDBACK_CODE + EVENT_OFF));
    if ~isempty(fb_off_events_after_fix)
        trial_end = EVENT.POS(fb_off_events_after_fix(1));
        trial_durations = [trial_durations, trial_end - fix_pos];
    end
end

% Use maximum trial duration to ensure we capture complete trials
if ~isempty(trial_durations)
    trial_length = round(max(trial_durations));
else
    % Default: 10 seconds (typical for MI-BCI tasks with fixation+cue+feedback)
    trial_length = round(10 * fs);
end

fprintf('Trial length: %d samples (%.2f seconds)\n', trial_length, trial_length/fs);
fprintf('  (from fixation cross to end of continuous feedback)\n');

% Initialize trial matrices for different signal types
% Raw signal trials
trials_raw = zeros(trial_length, num_channels, num_trials_extracted);
% μ band log BP trials
trials_mu_logbp = zeros(trial_length, num_channels, num_trials_extracted);
% β band log BP trials
trials_beta_logbp = zeros(trial_length, num_channels, num_trials_extracted);

% Create a vector Ck with cue information regarding each trial [trials x 1]
Ck = zeros(num_trials_extracted, 1);

% Extract each trial
valid_trials = 0;
for i = 1:num_trials_extracted
    fix_event_idx = fixation_event_indices(i);
    trial_start = EVENT.POS(fix_event_idx);  % Start at fixation cross
    trial_end = trial_start + trial_length - 1;
    
    % Find the cue event for this trial (occurs after fixation)
    cue_events_after_fix = find(EVENT.POS > trial_start & ...
                                 (EVENT.TYP == CUE_BOTH_HAND | EVENT.TYP == CUE_BOTH_FEET));
    
    % Check if trial fits within data bounds and has a cue
    if trial_end <= total_samples && ~isempty(cue_events_after_fix)
        % Extract trial data (from fixation to end of feedback)
        trials_raw(:, :, valid_trials + 1) = s(trial_start:trial_end, :);
        trials_mu_logbp(:, :, valid_trials + 1) = s_mu_logbp(trial_start:trial_end, :);
        trials_beta_logbp(:, :, valid_trials + 1) = s_beta_logbp(trial_start:trial_end, :);
        
        % Store cue type for this trial
        cue_idx = cue_events_after_fix(1);
        Ck(valid_trials + 1) = EVENT.TYP(cue_idx);
        
        valid_trials = valid_trials + 1;
    else
        if trial_end > total_samples
            fprintf('Warning: Trial %d extends beyond data bounds, skipping.\n', i);
        else
            fprintf('Warning: Trial %d has no cue event, skipping.\n', i);
        end
    end
end

% Trim arrays to valid trials only
trials_raw = trials_raw(:, :, 1:valid_trials);
trials_mu_logbp = trials_mu_logbp(:, :, 1:valid_trials);
trials_beta_logbp = trials_beta_logbp(:, :, 1:valid_trials);
Ck = Ck(1:valid_trials);

fprintf('Successfully extracted %d trials\n', valid_trials);
fprintf('Trial matrices size: [%d samples x %d channels x %d trials]\n', ...
        size(trials_raw, 1), size(trials_raw, 2), size(trials_raw, 3));

% Display cue distribution
fprintf('\nCue distribution:\n');
num_both_hand = sum(Ck == CUE_BOTH_HAND);
num_both_feet = sum(Ck == CUE_BOTH_FEET);
fprintf('  Both Hand (773): %d trials\n', num_both_hand);
fprintf('  Both Feet (771): %d trials\n', num_both_feet);
fprintf('===========================\n');

%% Visualization
fprintf('\n=== Visualization ===\n');

% Display all available channel labels
fprintf('\nAvailable EEG channels:\n');
if isfield(h, 'Label') && ~isempty(h.Label)
    for ch_idx = 1:min(length(h.Label), num_channels)
        fprintf('  Channel %2d: %s\n', ch_idx, h.Label{ch_idx});
    end
else
    fprintf('  No channel labels found in header.\n');
    for ch_idx = 1:num_channels
        fprintf('  Channel %2d\n', ch_idx);
    end
end

% a. Select three meaningful EEG channels
% For motor imagery, select channels over motor cortex areas
% Common choices: C3, Cz, C4 (central channels)
% Try to automatically find C3, Cz, C4 channels
if isfield(h, 'Label') && ~isempty(h.Label)
    % Search for C3, Cz, C4 in channel labels
    c3_idx = find(strcmpi(h.Label, 'C3'), 1);
    cz_idx = find(strcmpi(h.Label, 'Cz'), 1);
    c4_idx = find(strcmpi(h.Label, 'C4'), 1);
    
    if ~isempty(c3_idx) && ~isempty(cz_idx) && ~isempty(c4_idx)
        selected_channels = [c3_idx, cz_idx, c4_idx];
        channel_names = {'C3', 'Cz', 'C4'};
        fprintf('\nAutomatically found motor cortex channels:\n');
        fprintf('  C3: Channel %d\n', c3_idx);
        fprintf('  Cz: Channel %d\n', cz_idx);
        fprintf('  C4: Channel %d\n', c4_idx);
    else
        % Fallback to default channels
        selected_channels = [7, 9, 11];
        channel_names = {sprintf('Ch %d', 7), sprintf('Ch %d', 9), sprintf('Ch %d', 11)};
        fprintf('\nC3, Cz, C4 not found. Using default channels: %d, %d, %d\n', ...
                selected_channels(1), selected_channels(2), selected_channels(3));
        fprintf('Please verify these are appropriate motor cortex channels.\n');
    end
else
    % No labels available, use default
    selected_channels = [7, 9, 11];
    channel_names = {sprintf('Ch %d', 7), sprintf('Ch %d', 9), sprintf('Ch %d', 11)};
    fprintf('\nUsing default channels: %d, %d, %d\n', selected_channels(1), selected_channels(2), selected_channels(3));
end

fprintf('Selected channels for visualization: %d, %d, %d\n', selected_channels(1), selected_channels(2), selected_channels(3));

%% Visualization - Single Trial
% Plot trial 65 (Both Hand)
sample_trial_idx = 65;
fprintf('\nVisualizing trial #%d (cue type: %d)\n', sample_trial_idx, Ck(sample_trial_idx));

% Create time vector
time_trial = (0:trial_length-1) / fs;

% Create figure with 3 subplots in a single row
figure('Position', [100, 100, 1500, 400]);

% Subplot 1: Raw signal
subplot(1, 3, 1);
hold on;
for i = 1:3
    ch = selected_channels(i);
    plot(time_trial, trials_raw(:, ch, sample_trial_idx), 'DisplayName', channel_names{i});
end
hold off;
ylabel('Amplitude [\muV]');
xlabel('Time [s]');
title('Raw Signal - Trial 65');
legend(channel_names, 'Location', 'best');
grid on;
ylim(ylim + [-1 1] * range(ylim) * 0.1);

% Subplot 2: Filtered Mu band
subplot(1, 3, 2);
hold on;
for i = 1:3
    ch = selected_channels(i);
    % Extract filtered mu band signal for this trial
    trial_start = EVENT.POS(fixation_event_indices(sample_trial_idx));
    trial_end = trial_start + trial_length - 1;
    if trial_end <= total_samples
        s_mu_trial = s_mu(trial_start:trial_end, ch);
        plot(time_trial, s_mu_trial, 'DisplayName', channel_names{i});
    end
end
hold off;
ylabel('Amplitude [\muV]');
xlabel('Time [s]');
title('\mu Band (8-12 Hz) - Trial 65');
legend(channel_names, 'Location', 'best');
grid on;
ylim(ylim + [-1 1] * range(ylim) * 0.1);

% Subplot 3: Filtered Beta band
subplot(1, 3, 3);
hold on;
for i = 1:3
    ch = selected_channels(i);
    % Extract filtered beta band signal for this trial
    trial_start = EVENT.POS(fixation_event_indices(sample_trial_idx));
    trial_end = trial_start + trial_length - 1;
    if trial_end <= total_samples
        s_beta_trial = s_beta(trial_start:trial_end, ch);
        plot(time_trial, s_beta_trial, 'DisplayName', channel_names{i});
    end
end
hold off;
ylabel('Amplitude [\muV]');
xlabel('Time [s]');
title('\beta Band (18-22 Hz) - Trial 65');
legend(channel_names, 'Location', 'best');
grid on;
ylim(ylim + [-1 1] * range(ylim) * 0.1);

sgtitle('Single Trial Visualization - Channels 7, 9, 11');

% c. Compute the power of the three channels after averaging across all trials for each class
fprintf('Computing averaged power across trials for each class...\n');

% Get trial indices for each class
hand_trials = find(Ck == CUE_BOTH_HAND);
feet_trials = find(Ck == CUE_BOTH_FEET);

% Extract log band power for each class
trials_mu_hand = trials_mu_logbp(:, :, hand_trials);
trials_mu_feet = trials_mu_logbp(:, :, feet_trials);
trials_beta_hand = trials_beta_logbp(:, :, hand_trials);
trials_beta_feet = trials_beta_logbp(:, :, feet_trials);

% Compute average power across trials for each class
avg_mu_hand = mean(trials_mu_hand, 3);    % [samples x channels]
avg_mu_feet = mean(trials_mu_feet, 3);
avg_beta_hand = mean(trials_beta_hand, 3);
avg_beta_feet = mean(trials_beta_feet, 3);

fprintf('Averaged across %d Both Hand trials and %d Both Feet trials\n', ...
        length(hand_trials), length(feet_trials));

%% Visualization - Averaged Band Power across Classes
% Create figure with 6 subplots (2 rows x 3 columns)
% First row: Mu band averaged power for each channel (both classes)
% Second row: Beta band averaged power for each channel (both classes)
figure('Position', [100, 100, 1500, 600]);

% First row: Mu band averaged power for each of the 3 channels
for i = 1:3
    ch = selected_channels(i);
    subplot(2, 3, i);
    hold on;
    plot(time_trial, avg_mu_hand(:, ch), 'r', 'DisplayName', 'Both Hands');
    plot(time_trial, avg_mu_feet(:, ch), 'b', 'DisplayName', 'Both Feet');
    hold off;
    ylabel('Log Power');
    xlabel('Time [s]');
    title(sprintf('%s - \\mu Band', channel_names{i}));
    legend('Location', 'best');
    grid on;
    ylim(ylim + [-1 1] * range(ylim) * 0.1);
end

% Second row: Beta band averaged power for each of the 3 channels
for i = 1:3
    ch = selected_channels(i);
    subplot(2, 3, i + 3);
    hold on;
    plot(time_trial, avg_beta_hand(:, ch), 'r', 'DisplayName', 'Both Hands');
    plot(time_trial, avg_beta_feet(:, ch), 'b', 'DisplayName', 'Both Feet');
    hold off;
    ylabel('Log Power');
    xlabel('Time [s]');
    title(sprintf('%s - \\beta Band', channel_names{i}));
    legend('Location', 'best');
    grid on;
    ylim(ylim + [-1 1] * range(ylim) * 0.1);
end

sgtitle('Averaged Logarithmic Band Power - Both Hand vs Both Feet');

% % Original commented plotting code below
% % figure('Position', [100, 100, 1400, 900]);
% 
% for i = 1:3
%     ch = selected_channels(i);
%     
%     % μ band - Both Hand
%     subplot(3, 4, (i-1)*4 + 1);
%     plot(time_trial, avg_mu_hand(:, ch), 'b', 'LineWidth', 2);
%     ylabel('Log Power');
%     title(sprintf('%s - μ Band Both Hand', channel_names{i}));
%     grid on;
%     
%     % μ band - Both Feet
%     subplot(3, 4, (i-1)*4 + 2);
%     plot(time_trial, avg_mu_feet(:, ch), 'r', 'LineWidth', 2);
%     ylabel('Log Power');
%     title(sprintf('%s - μ Band Both Feet', channel_names{i}));
%     grid on;
%     
%     % β band - Both Hand
%     subplot(3, 4, (i-1)*4 + 3);
%     plot(time_trial, avg_beta_hand(:, ch), 'b', 'LineWidth', 2);
%     ylabel('Log Power');
%     title(sprintf('%s - β Band Both Hand', channel_names{i}));
%     grid on;
%     
%     % β band - Both Feet
%     subplot(3, 4, (i-1)*4 + 4);
%     plot(time_trial, avg_beta_feet(:, ch), 'r', 'LineWidth', 2);
%     ylabel('Log Power');
%     title(sprintf('%s - β Band Both Feet', channel_names{i}));
%     grid on;
% end
% 
% % Add x-labels to bottom row
% for j = 9:12
%     subplot(3, 4, j);
%     xlabel('Time (s)');
% end
% 
% sgtitle(sprintf('Averaged Log Band Power Across Trials (n_hand=%d, n_feet=%d)', ...
%         length(hand_trials), length(feet_trials)));
% 
% % Create comparison plots for each channel (overlay both classes)
% figure('Position', [100, 100, 1400, 600]);
% 
% for i = 1:3
%     ch = selected_channels(i);
%     
%     % μ band comparison
%     subplot(2, 3, i);
%     plot(time_trial, avg_mu_hand(:, ch), 'b', 'LineWidth', 2);
%     hold on;
%     plot(time_trial, avg_mu_feet(:, ch), 'r', 'LineWidth', 2);
%     hold off;
%     ylabel('Log Power');
%     xlabel('Time (s)');
%     title(sprintf('%s - μ Band (10-12 Hz)', channel_names{i}));
%     legend('Both Hand', 'Both Feet', 'Location', 'best');
%     grid on;
%     
%     % β band comparison
%     subplot(2, 3, i + 3);
%     plot(time_trial, avg_beta_hand(:, ch), 'b', 'LineWidth', 2);
%     hold on;
%     plot(time_trial, avg_beta_feet(:, ch), 'r', 'LineWidth', 2);
%     hold off;
%     ylabel('Log Power');
%     xlabel('Time (s)');
%     title(sprintf('%s - β Band (18-24 Hz)', channel_names{i}));
%     legend('Both Hand', 'Both Feet', 'Location', 'best');
%     grid on;
% end
% 
% sgtitle('Comparison of Averaged Log Band Power Between Classes (No Spatial Filter)');

fprintf('Visualization complete!\n');
fprintf('===========================\n');

%% Process CAR-Filtered Data
fprintf('\n\n==========================================================\n');
fprintf('=== PROCESSING WITH CAR FILTER ===\n');
fprintf('==========================================================\n\n');

% Use CAR-filtered data
s_process = s_car;
num_channels_process = size(s_process, 2);

%% Signal Processing - Logarithmic Band Power (CAR)
fprintf('\n=== Computing Logarithmic Band Power (CAR Filter) ===\n');

% Apply bandpass filtering
s_mu_car = zeros(size(s_process));
s_beta_car = zeros(size(s_process));

fprintf('Applying zero-phase filtering to all channels (CAR)...\n');
for ch = 1:num_channels_process
    s_mu_car(:, ch) = filtfilt(b_mu, a_mu, s_process(:, ch));
    s_beta_car(:, ch) = filtfilt(b_beta, a_beta, s_process(:, ch));
end

% Rectify
s_mu_squared_car = s_mu_car .^ 2;
s_beta_squared_car = s_beta_car .^ 2;

% Moving average
s_mu_smooth_car = zeros(size(s_mu_squared_car));
s_beta_smooth_car = zeros(size(s_beta_squared_car));
for ch = 1:num_channels_process
    s_mu_smooth_car(:, ch) = conv(s_mu_squared_car(:, ch), ma_filter, 'same');
    s_beta_smooth_car(:, ch) = conv(s_beta_squared_car(:, ch), ma_filter, 'same');
end

% Log transform
s_mu_logbp_car = log(s_mu_smooth_car + epsilon);
s_beta_logbp_car = log(s_beta_smooth_car + epsilon);

fprintf('Log band power computation complete (CAR)!\n');

%% Trial Extraction (CAR)
fprintf('\n=== Extracting Trials (CAR Filter) ===\n');

trials_raw_car = zeros(trial_length, num_channels_process, num_trials_extracted);
trials_mu_logbp_car = zeros(trial_length, num_channels_process, num_trials_extracted);
trials_beta_logbp_car = zeros(trial_length, num_channels_process, num_trials_extracted);
Ck_car = zeros(num_trials_extracted, 1);

valid_trials_car = 0;
for i = 1:num_trials_extracted
    fix_event_idx = fixation_event_indices(i);
    trial_start = EVENT.POS(fix_event_idx);
    trial_end = trial_start + trial_length - 1;
    
    cue_events_after_fix = find(EVENT.POS > trial_start & ...
                                 (EVENT.TYP == CUE_BOTH_HAND | EVENT.TYP == CUE_BOTH_FEET));
    
    if trial_end <= total_samples && ~isempty(cue_events_after_fix)
        trials_raw_car(:, :, valid_trials_car + 1) = s_process(trial_start:trial_end, :);
        trials_mu_logbp_car(:, :, valid_trials_car + 1) = s_mu_logbp_car(trial_start:trial_end, :);
        trials_beta_logbp_car(:, :, valid_trials_car + 1) = s_beta_logbp_car(trial_start:trial_end, :);
        
        cue_idx = cue_events_after_fix(1);
        Ck_car(valid_trials_car + 1) = EVENT.TYP(cue_idx);
        
        valid_trials_car = valid_trials_car + 1;
    end
end

trials_raw_car = trials_raw_car(:, :, 1:valid_trials_car);
trials_mu_logbp_car = trials_mu_logbp_car(:, :, 1:valid_trials_car);
trials_beta_logbp_car = trials_beta_logbp_car(:, :, 1:valid_trials_car);
Ck_car = Ck_car(1:valid_trials_car);

fprintf('Successfully extracted %d trials (CAR)\n', valid_trials_car);

%% Visualization (CAR)
fprintf('\n=== Visualization (CAR Filter) ===\n');

hand_trials_car = find(Ck_car == CUE_BOTH_HAND);
feet_trials_car = find(Ck_car == CUE_BOTH_FEET);

% Averaged power
trials_mu_hand_car = trials_mu_logbp_car(:, :, hand_trials_car);
trials_mu_feet_car = trials_mu_logbp_car(:, :, feet_trials_car);
trials_beta_hand_car = trials_beta_logbp_car(:, :, hand_trials_car);
trials_beta_feet_car = trials_beta_logbp_car(:, :, feet_trials_car);

avg_mu_hand_car = mean(trials_mu_hand_car, 3);
avg_mu_feet_car = mean(trials_mu_feet_car, 3);
avg_beta_hand_car = mean(trials_beta_hand_car, 3);
avg_beta_feet_car = mean(trials_beta_feet_car, 3);

%% Visualization - CAR Filtered Data
% Figure 3: Single Trial with CAR filter
figure('Position', [100, 100, 1500, 400]);

% Subplot 1: Raw signal (CAR filtered)
subplot(1, 3, 1);
hold on;
for i = 1:3
    ch = selected_channels(i);
    plot(time_trial, trials_raw_car(:, ch, sample_trial_idx), 'DisplayName', channel_names{i});
end
hold off;
ylabel('Amplitude [\muV]');
xlabel('Time [s]');
title('Raw Signal (CAR) - Trial 65');
legend(channel_names, 'Location', 'best');
grid on;
ylim(ylim + [-1 1] * range(ylim) * 0.1);

% Subplot 2: Filtered Mu band (CAR)
subplot(1, 3, 2);
hold on;
for i = 1:3
    ch = selected_channels(i);
    % Extract filtered mu band signal for this trial
    trial_start = EVENT.POS(fixation_event_indices(sample_trial_idx));
    trial_end = trial_start + trial_length - 1;
    if trial_end <= total_samples
        s_mu_trial_car = s_mu_car(trial_start:trial_end, ch);
        plot(time_trial, s_mu_trial_car, 'DisplayName', channel_names{i});
    end
end
hold off;
ylabel('Amplitude [\muV]');
xlabel('Time [s]');
title('\mu Band (8-12 Hz) (CAR) - Trial 65');
legend(channel_names, 'Location', 'best');
grid on;
ylim(ylim + [-1 1] * range(ylim) * 0.1);

% Subplot 3: Filtered Beta band (CAR)
subplot(1, 3, 3);
hold on;
for i = 1:3
    ch = selected_channels(i);
    % Extract filtered beta band signal for this trial
    trial_start = EVENT.POS(fixation_event_indices(sample_trial_idx));
    trial_end = trial_start + trial_length - 1;
    if trial_end <= total_samples
        s_beta_trial_car = s_beta_car(trial_start:trial_end, ch);
        plot(time_trial, s_beta_trial_car, 'DisplayName', channel_names{i});
    end
end
hold off;
ylabel('Amplitude [\muV]');
xlabel('Time [s]');
title('\beta Band (18-22 Hz) (CAR) - Trial 65');
legend(channel_names, 'Location', 'best');
grid on;
ylim(ylim + [-1 1] * range(ylim) * 0.1);

sgtitle('Single Trial Visualization (CAR Filter) - Channels 7, 9, 11');

% Figure 4: Averaged Band Power with CAR filter
figure('Position', [100, 100, 1500, 600]);

% First row: Mu band averaged power for each of the 3 channels
for i = 1:3
    ch = selected_channels(i);
    subplot(2, 3, i);
    hold on;
    plot(time_trial, avg_mu_hand_car(:, ch), 'r', 'DisplayName', 'Both Hands');
    plot(time_trial, avg_mu_feet_car(:, ch), 'b', 'DisplayName', 'Both Feet');
    hold off;
    ylabel('Log Power');
    xlabel('Time [s]');
    title(sprintf('%s - \\mu Band (CAR)', channel_names{i}));
    legend('Location', 'best');
    grid on;
    ylim(ylim + [-1 1] * range(ylim) * 0.1);
end

% Second row: Beta band averaged power for each of the 3 channels
for i = 1:3
    ch = selected_channels(i);
    subplot(2, 3, i + 3);
    hold on;
    plot(time_trial, avg_beta_hand_car(:, ch), 'r', 'DisplayName', 'Both Hands');
    plot(time_trial, avg_beta_feet_car(:, ch), 'b', 'DisplayName', 'Both Feet');
    hold off;
    ylabel('Log Power');
    xlabel('Time [s]');
    title(sprintf('%s - \\beta Band (CAR)', channel_names{i}));
    legend('Location', 'best');
    grid on;
    ylim(ylim + [-1 1] * range(ylim) * 0.1);
end

sgtitle('Averaged Logarithmic Band Power (CAR Filter) - Both Hand vs Both Feet');

fprintf('CAR filter processing complete!\n');

%% Process Laplacian-Filtered Data
fprintf('\n\n==========================================================\n');
fprintf('=== PROCESSING WITH LAPLACIAN FILTER ===\n');
fprintf('==========================================================\n\n');

s_process = s_lap;
num_channels_process = size(s_process, 2);

%% Signal Processing - Logarithmic Band Power (Laplacian)
fprintf('\n=== Computing Logarithmic Band Power (Laplacian Filter) ===\n');

s_mu_lap = zeros(size(s_process));
s_beta_lap = zeros(size(s_process));

fprintf('Applying zero-phase filtering to all channels (Laplacian)...\n');
for ch = 1:num_channels_process
    s_mu_lap(:, ch) = filtfilt(b_mu, a_mu, s_process(:, ch));
    s_beta_lap(:, ch) = filtfilt(b_beta, a_beta, s_process(:, ch));
end

s_mu_squared_lap = s_mu_lap .^ 2;
s_beta_squared_lap = s_beta_lap .^ 2;

s_mu_smooth_lap = zeros(size(s_mu_squared_lap));
s_beta_smooth_lap = zeros(size(s_beta_squared_lap));
for ch = 1:num_channels_process
    s_mu_smooth_lap(:, ch) = conv(s_mu_squared_lap(:, ch), ma_filter, 'same');
    s_beta_smooth_lap(:, ch) = conv(s_beta_squared_lap(:, ch), ma_filter, 'same');
end

s_mu_logbp_lap = log(s_mu_smooth_lap + epsilon);
s_beta_logbp_lap = log(s_beta_smooth_lap + epsilon);

fprintf('Log band power computation complete (Laplacian)!\n');

%% Trial Extraction (Laplacian)
fprintf('\n=== Extracting Trials (Laplacian Filter) ===\n');

trials_raw_lap = zeros(trial_length, num_channels_process, num_trials_extracted);
trials_mu_logbp_lap = zeros(trial_length, num_channels_process, num_trials_extracted);
trials_beta_logbp_lap = zeros(trial_length, num_channels_process, num_trials_extracted);
Ck_lap = zeros(num_trials_extracted, 1);

valid_trials_lap = 0;
for i = 1:num_trials_extracted
    fix_event_idx = fixation_event_indices(i);
    trial_start = EVENT.POS(fix_event_idx);
    trial_end = trial_start + trial_length - 1;
    
    cue_events_after_fix = find(EVENT.POS > trial_start & ...
                                 (EVENT.TYP == CUE_BOTH_HAND | EVENT.TYP == CUE_BOTH_FEET));
    
    if trial_end <= total_samples && ~isempty(cue_events_after_fix)
        trials_raw_lap(:, :, valid_trials_lap + 1) = s_process(trial_start:trial_end, :);
        trials_mu_logbp_lap(:, :, valid_trials_lap + 1) = s_mu_logbp_lap(trial_start:trial_end, :);
        trials_beta_logbp_lap(:, :, valid_trials_lap + 1) = s_beta_logbp_lap(trial_start:trial_end, :);
        
        cue_idx = cue_events_after_fix(1);
        Ck_lap(valid_trials_lap + 1) = EVENT.TYP(cue_idx);
        
        valid_trials_lap = valid_trials_lap + 1;
    end
end

trials_raw_lap = trials_raw_lap(:, :, 1:valid_trials_lap);
trials_mu_logbp_lap = trials_mu_logbp_lap(:, :, 1:valid_trials_lap);
trials_beta_logbp_lap = trials_beta_logbp_lap(:, :, 1:valid_trials_lap);
Ck_lap = Ck_lap(1:valid_trials_lap);

fprintf('Successfully extracted %d trials (Laplacian)\n', valid_trials_lap);

%% Visualization (Laplacian)
fprintf('\n=== Visualization (Laplacian Filter) ===\n');

hand_trials_lap = find(Ck_lap == CUE_BOTH_HAND);
feet_trials_lap = find(Ck_lap == CUE_BOTH_FEET);

% Averaged power
trials_mu_hand_lap = trials_mu_logbp_lap(:, :, hand_trials_lap);
trials_mu_feet_lap = trials_mu_logbp_lap(:, :, feet_trials_lap);
trials_beta_hand_lap = trials_beta_logbp_lap(:, :, hand_trials_lap);
trials_beta_feet_lap = trials_beta_logbp_lap(:, :, feet_trials_lap);

avg_mu_hand_lap = mean(trials_mu_hand_lap, 3);
avg_mu_feet_lap = mean(trials_mu_feet_lap, 3);
avg_beta_hand_lap = mean(trials_beta_hand_lap, 3);
avg_beta_feet_lap = mean(trials_beta_feet_lap, 3);

%% Visualization - Laplacian Filtered Data
% Figure 5: Single Trial with Laplacian filter
figure('Position', [100, 100, 1500, 400]);

% Subplot 1: Raw signal (Laplacian filtered)
subplot(1, 3, 1);
hold on;
for i = 1:3
    ch = selected_channels(i);
    plot(time_trial, trials_raw_lap(:, ch, sample_trial_idx), 'DisplayName', channel_names{i});
end
hold off;
ylabel('Amplitude [\muV]');
xlabel('Time [s]');
title('Raw Signal (Laplacian) - Trial 65');
legend(channel_names, 'Location', 'best');
grid on;
ylim(ylim + [-1 1] * range(ylim) * 0.1);

% Subplot 2: Filtered Mu band (Laplacian)
subplot(1, 3, 2);
hold on;
for i = 1:3
    ch = selected_channels(i);
    % Extract filtered mu band signal for this trial
    trial_start = EVENT.POS(fixation_event_indices(sample_trial_idx));
    trial_end = trial_start + trial_length - 1;
    if trial_end <= total_samples
        s_mu_trial_lap = s_mu_lap(trial_start:trial_end, ch);
        plot(time_trial, s_mu_trial_lap, 'DisplayName', channel_names{i});
    end
end
hold off;
ylabel('Amplitude [\muV]');
xlabel('Time [s]');
title('\mu Band (8-12 Hz) (Laplacian) - Trial 65');
legend(channel_names, 'Location', 'best');
grid on;
ylim(ylim + [-1 1] * range(ylim) * 0.1);

% Subplot 3: Filtered Beta band (Laplacian)
subplot(1, 3, 3);
hold on;
for i = 1:3
    ch = selected_channels(i);
    % Extract filtered beta band signal for this trial
    trial_start = EVENT.POS(fixation_event_indices(sample_trial_idx));
    trial_end = trial_start + trial_length - 1;
    if trial_end <= total_samples
        s_beta_trial_lap = s_beta_lap(trial_start:trial_end, ch);
        plot(time_trial, s_beta_trial_lap, 'DisplayName', channel_names{i});
    end
end
hold off;
ylabel('Amplitude [\muV]');
xlabel('Time [s]');
title('\beta Band (18-22 Hz) (Laplacian) - Trial 65');
legend(channel_names, 'Location', 'best');
grid on;
ylim(ylim + [-1 1] * range(ylim) * 0.1);

sgtitle('Single Trial Visualization (Laplacian Filter) - Channels 7, 9, 11');

% Figure 6: Averaged Band Power with Laplacian filter
figure('Position', [100, 100, 1500, 600]);

% First row: Mu band averaged power for each of the 3 channels
for i = 1:3
    ch = selected_channels(i);
    subplot(2, 3, i);
    hold on;
    plot(time_trial, avg_mu_hand_lap(:, ch), 'r', 'DisplayName', 'Both Hands');
    plot(time_trial, avg_mu_feet_lap(:, ch), 'b', 'DisplayName', 'Both Feet');
    hold off;
    ylabel('Log Power');
    xlabel('Time [s]');
    title(sprintf('%s - \\mu Band (Laplacian)', channel_names{i}));
    legend('Location', 'best');
    grid on;
    ylim(ylim + [-1 1] * range(ylim) * 0.1);
end

% Second row: Beta band averaged power for each of the 3 channels
for i = 1:3
    ch = selected_channels(i);
    subplot(2, 3, i + 3);
    hold on;
    plot(time_trial, avg_beta_hand_lap(:, ch), 'r', 'DisplayName', 'Both Hands');
    plot(time_trial, avg_beta_feet_lap(:, ch), 'b', 'DisplayName', 'Both Feet');
    hold off;
    ylabel('Log Power');
    xlabel('Time [s]');
    title(sprintf('%s - \\beta Band (Laplacian)', channel_names{i}));
    legend('Location', 'best');
    grid on;
    ylim(ylim + [-1 1] * range(ylim) * 0.1);
end

sgtitle('Averaged Logarithmic Band Power (Laplacian Filter) - Both Hand vs Both Feet');

fprintf('Laplacian filter processing complete!\n');

fprintf('\n\n==========================================================\n');
fprintf('=== ALL PROCESSING COMPLETE ===\n');
fprintf('Processed data with:\n');
fprintf('  1. No spatial filter\n');
fprintf('  2. CAR filter\n');
fprintf('  3. Laplacian filter\n');
fprintf('==========================================================\n\n');

