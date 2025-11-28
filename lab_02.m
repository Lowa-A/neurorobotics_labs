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

% Get sampling rate from header
fs = h.SampleRate;  % Sampling frequency in Hz

% Select 3 EEG channels
channels = [1, 2, 3];  % Change these to plot different channels

% Calculate number of samples for 5 seconds
num_samples = 5 * fs;

% Extract 5 seconds of data (starting from sample 1)
if size(s, 1) >= num_samples
    data_5sec = s(1:num_samples, channels);
else
    data_5sec = s(:, channels);  % Use all data if less than 5 seconds
    num_samples = size(data_5sec, 1);
end

% Create time vector in seconds
time = (0:num_samples-1) / fs;

% Calculate common y-axis limits for all channels
y_min = min(data_5sec(:));
y_max = max(data_5sec(:));
y_range = y_max - y_min;
y_limits = [y_min - 0.1*y_range, y_max + 0.1*y_range];

% Create figure with 3 subplots
figure;

for i = 1:3
    subplot(3, 1, i);
    plot(time, data_5sec(:, i), 'b', 'LineWidth', 1);
    xlabel('Time (seconds)');
    ylabel('Amplitude (\muV)');
    
    % Set common y-axis limits
    ylim(y_limits);
    
    % Add channel label if available
    if isfield(h, 'Label')
        title(sprintf('Channel %d: %s - 5 seconds', channels(i), h.Label{channels(i)}));
    else
        title(sprintf('EEG Channel %d - 5 seconds', channels(i)));
    end
    
    grid on;
end

% Adjust spacing between subplots
sgtitle('5 Seconds of EEG Data - 3 Channels');

%% Lab 02 - Creation of Label Vectors
% Extract event information from the GDF file
EVENT = h.EVENT;  % Event structure from header

% Get total number of samples in the recording
total_samples = size(s, 1);

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

% f. Plot the label vectors
% figure;
% time_vector = (0:total_samples-1) / fs;
% 
% % Plot 1: Run index (Rk)
% subplot(6, 1, 1);
% plot(time_vector, Rk, 'LineWidth', 1.5);
% ylabel('Run #');
% title('Rk: Run Index Vector');
% grid on;
% ylim([0 max(Rk)+1]);
% 
% % Plot 2: Trial vector (Tk)
% subplot(6, 1, 2);
% plot(time_vector, Tk, 'LineWidth', 1.5);
% ylabel('Trial #');
% title('Tk: Trial Vector');
% grid on;
% ylim([0 max(Tk)+1]);
% 
% % Plot 3: Fixation periods (Fk)
% subplot(6, 1, 3);
% plot(time_vector, Fk, 'LineWidth', 1.5);
% ylabel('Event Code');
% title('Fk: Fixation Periods');
% grid on;
% 
% % Plot 4: Cue periods (Ak)
% subplot(5, 1, 3);
% subplot(6, 1, 4);
% plot(time_vector, Ak, 'LineWidth', 1.5);
% ylabel('Event Code');
% title('Ak: Cue Periods (773=Both Hand, 771=Both Feet, 783=Rest)');
% grid on;
% 
% % Plot 5: Continuous feedback periods (CFk)
% subplot(6, 1, 5);
% plot(time_vector, CFk, 'LineWidth', 1.5);
% ylabel('Event Code');
% title('CFk: Continuous Feedback Periods');
% grid on;
% 
% % Plot 6: Hit/Miss periods (Xk)
% subplot(6, 1, 6);
% plot(time_vector, Xk, 'LineWidth', 1.5);
% ylabel('Event Code');
% title('Xk: Hit/Miss Periods (897=Hit, 898=Miss)');
% xlabel('Time (seconds)');
% grid on;
% 
% sgtitle('Label Vectors for Concatenated GDF Files');

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

%% Trial Extraction
fprintf('\n=== Extracting Trials ===\n');

% Extract trials corresponding to the two motor imagery (MI) tasks
% a. Use the event types (EVENT.TYP) to identify which tasks correspond to which class
fprintf('Identifying MI task classes from events...\n');
fprintf('  Class 1: Both Hand (773)\n');
fprintf('  Class 2: Both Feet (771)\n');

% b. Trial lasts from fixation cross to end of continuous feedback period
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

% Get number of EEG channels
num_channels = size(s, 2);

% Initialize trial matrices [samples x channels x trials]
trials = zeros(trial_length, num_channels, num_trials_extracted);

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
        trials(:, :, valid_trials + 1) = s(trial_start:trial_end, :);
        
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
trials = trials(:, :, 1:valid_trials);
Ck = Ck(1:valid_trials);

fprintf('Successfully extracted %d trials\n', valid_trials);
fprintf('Trial matrix size: [%d samples x %d channels x %d trials]\n', ...
        size(trials, 1), size(trials, 2), size(trials, 3));

% Display cue distribution
fprintf('\nCue distribution:\n');
num_both_hand = sum(Ck == CUE_BOTH_HAND);
num_both_feet = sum(Ck == CUE_BOTH_FEET);
fprintf('  Both Hand (773): %d trials\n', num_both_hand);
fprintf('  Both Feet (771): %d trials\n', num_both_feet);
fprintf('===========================\n');

%% Signal Processing - Logarithmic Band Power
fprintf('\n=== Computing Logarithmic Band Power ===\n');

% Define frequency bands
mu_band = [8 12];     % μ (mu) band: 10-12 Hz
beta_band = [18 22];   % β (beta) band: 18-24 Hz

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

% Initialize matrices for filtered trials
trials_mu = zeros(size(trials));
trials_beta = zeros(size(trials));

fprintf('Applying zero-phase filtering to all trials...\n');

% Apply zero-phase filtering using filtfilt for each trial and channel
for trial = 1:size(trials, 3)
    for ch = 1:size(trials, 2)
        % μ band
        trials_mu(:, ch, trial) = filtfilt(b_mu, a_mu, trials(:, ch, trial));
        
        % β band
        trials_beta(:, ch, trial) = filtfilt(b_beta, a_beta, trials(:, ch, trial));
    end
end

fprintf('Filtering complete.\n');

% b. Rectify the signal (square it)
fprintf('Rectifying signals (squaring)...\n');
trials_mu_squared = trials_mu .^ 2;
trials_beta_squared = trials_beta .^ 2;

% c. Apply a moving average using a 1-second window
window_size = round(fs);  % 1 second window
fprintf('Applying moving average (window size: %d samples = 1 second)...\n', window_size);

% Create moving average filter
ma_filter = ones(window_size, 1) / window_size;

% Initialize matrices for smoothed signals
trials_mu_smooth = zeros(size(trials_mu_squared));
trials_beta_smooth = zeros(size(trials_beta_squared));

% Apply moving average to each trial and channel
for trial = 1:size(trials, 3)
    for ch = 1:size(trials, 2)
        % μ band - use 'same' to maintain length
        trials_mu_smooth(:, ch, trial) = conv(trials_mu_squared(:, ch, trial), ma_filter, 'same');
        
        % β band
        trials_beta_smooth(:, ch, trial) = conv(trials_beta_squared(:, ch, trial), ma_filter, 'same');
    end
end

% d. Apply a logarithmic transform to the result
fprintf('Applying logarithmic transform...\n');

% Add small epsilon to avoid log(0)
epsilon = 1e-10;
trials_mu_logbp = log(trials_mu_smooth + epsilon);
trials_beta_logbp = log(trials_beta_smooth + epsilon);

fprintf('Logarithmic band power computation complete!\n');
fprintf('  μ band log BP size: [%d x %d x %d]\n', size(trials_mu_logbp, 1), size(trials_mu_logbp, 2), size(trials_mu_logbp, 3));
fprintf('  β band log BP size: [%d x %d x %d]\n', size(trials_beta_logbp, 1), size(trials_beta_logbp, 2), size(trials_beta_logbp, 3));
fprintf('===========================\n');

% Visualize the processing pipeline for one trial and channel
sample_trial = 1;
sample_channel = 3;

figure;
time_trial = (0:trial_length-1) / fs;

% Original signal
subplot(5, 2, 1);
plot(time_trial, trials(:, sample_channel, sample_trial), 'k', 'LineWidth', 1);
ylabel('Amplitude');
title(sprintf('Original Signal (Ch %d, Trial %d)', sample_channel, sample_trial));
grid on;

subplot(5, 2, 2);
plot(time_trial, trials(:, sample_channel, sample_trial), 'k', 'LineWidth', 1);
ylabel('Amplitude');
title(sprintf('Original Signal (Ch %d, Trial %d)', sample_channel, sample_trial));
grid on;

% Filtered signals
subplot(5, 2, 3);
plot(time_trial, trials_mu(:, sample_channel, sample_trial), 'b', 'LineWidth', 1);
ylabel('Amplitude');
title(sprintf('μ Band Filtered (%.0f-%.0f Hz)', mu_band(1), mu_band(2)));
grid on;

subplot(5, 2, 4);
plot(time_trial, trials_beta(:, sample_channel, sample_trial), 'r', 'LineWidth', 1);
ylabel('Amplitude');
title(sprintf('β Band Filtered (%.0f-%.0f Hz)', beta_band(1), beta_band(2)));
grid on;

% Squared signals
subplot(5, 2, 5);
plot(time_trial, trials_mu_squared(:, sample_channel, sample_trial), 'b', 'LineWidth', 1);
ylabel('Power');
title('μ Band Squared');
grid on;

subplot(5, 2, 6);
plot(time_trial, trials_beta_squared(:, sample_channel, sample_trial), 'r', 'LineWidth', 1);
ylabel('Power');
title('β Band Squared');
grid on;

% Smoothed signals
subplot(5, 2, 7);
plot(time_trial, trials_mu_smooth(:, sample_channel, sample_trial), 'b', 'LineWidth', 1);
ylabel('Power');
title('μ Band Smoothed (1s MA)');
grid on;

subplot(5, 2, 8);
plot(time_trial, trials_beta_smooth(:, sample_channel, sample_trial), 'r', 'LineWidth', 1);
ylabel('Power');
title('β Band Smoothed (1s MA)');
grid on;

% Log BP
subplot(5, 2, 9);
plot(time_trial, trials_mu_logbp(:, sample_channel, sample_trial), 'b', 'LineWidth', 1.5);
ylabel('Log Power');
title('μ Band Log BP');
xlabel('Time (s)');
grid on;

subplot(5, 2, 10);
plot(time_trial, trials_beta_logbp(:, sample_channel, sample_trial), 'r', 'LineWidth', 1.5);
ylabel('Log Power');
title('β Band Log BP');
xlabel('Time (s)');
grid on;

sgtitle(sprintf('Signal Processing Pipeline - Channel %d, Trial %d', sample_channel, sample_trial));

% Visualize a sample trial from each class
figure;
time_trial = (0:trial_length-1) / fs;

% Find first trial of each class
trial_hand_idx = find(Ck == CUE_BOTH_HAND, 1);
trial_feet_idx = find(Ck == CUE_BOTH_FEET, 1);

% Plot 3 channels for Both Hand trial
if ~isempty(trial_hand_idx)
    for ch = 1:3
        subplot(3, 2, ch*2-1);
        plot(time_trial, trials(:, ch, trial_hand_idx), 'b', 'LineWidth', 1);
        ylabel(sprintf('Ch %d (\\muV)', ch));
        if ch == 1
            title('Both Hand Trial (773)');
        end
        grid on;
        if ch == 3
            xlabel('Time (s)');
        end
    end
end

% Plot 3 channels for Both Feet trial
if ~isempty(trial_feet_idx)
    for ch = 1:3
        subplot(3, 2, ch*2);
        plot(time_trial, trials(:, ch, trial_feet_idx), 'r', 'LineWidth', 1);
        ylabel(sprintf('Ch %d (\\muV)', ch));
        if ch == 1
            title('Both Feet Trial (771)');
        end
        grid on;
        if ch == 3
            xlabel('Time (s)');
        end
    end
end

sgtitle('Sample Trials - First 3 Channels');

%% Visualization
fprintf('\n=== Visualization ===\n');

% Select channel to visualize (you can change this)
channel_to_plot = 3;  % Channel 3
fprintf('Selected channel for visualization: %d\n', channel_to_plot);

% a. Select a trial for each cue and plot a channel of your choice
% Get indices for each cue type
hand_trials = find(Ck == CUE_BOTH_HAND);
feet_trials = find(Ck == CUE_BOTH_FEET);

% Select one trial from each class (e.g., first trial)
selected_hand_trial = hand_trials(1);
selected_feet_trial = feet_trials(1);

% Create time vector for trials
time_trial = (0:trial_length-1) / fs;

% Plot individual trials for selected channel
figure;

subplot(2, 1, 1);
plot(time_trial, trials(:, channel_to_plot, selected_hand_trial), 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
title(sprintf('Single Trial - Both Hand (Channel %d, Trial #%d)', channel_to_plot, selected_hand_trial));
grid on;

subplot(2, 1, 2);
plot(time_trial, trials(:, channel_to_plot, selected_feet_trial), 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
title(sprintf('Single Trial - Both Feet (Channel %d, Trial #%d)', channel_to_plot, selected_feet_trial));
grid on;

sgtitle(sprintf('Individual Trials - Channel %d', channel_to_plot));

% b. Compute the grand averages for each cue and plot a channel of your choice
fprintf('Computing grand averages...\n');

% Extract trials for each class
trials_hand = trials(:, :, hand_trials);  % [samples x channels x trials]
trials_feet = trials(:, :, feet_trials);  % [samples x channels x trials]

% Compute grand average (mean across trials)
grand_avg_hand = mean(trials_hand, 3);  % [samples x channels]
grand_avg_feet = mean(trials_feet, 3);  % [samples x channels]

% Compute standard error for visualization (optional)
std_hand = std(trials_hand, 0, 3);
std_feet = std(trials_feet, 0, 3);
sem_hand = std_hand / sqrt(length(hand_trials));
sem_feet = std_feet / sqrt(length(feet_trials));

fprintf('Grand averages computed:\n');
fprintf('  Both Hand: averaged over %d trials\n', length(hand_trials));
fprintf('  Both Feet: averaged over %d trials\n', length(feet_trials));

% Plot grand averages for selected channel
figure;

subplot(2, 1, 1);
plot(time_trial, grand_avg_hand(:, channel_to_plot), 'b', 'LineWidth', 2);
hold on;
% Add shaded error region (standard error)
fill([time_trial, fliplr(time_trial)], ...
     [grand_avg_hand(:, channel_to_plot)' + sem_hand(:, channel_to_plot)', ...
      fliplr(grand_avg_hand(:, channel_to_plot)' - sem_hand(:, channel_to_plot)')], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold off;
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
title(sprintf('Grand Average - Both Hand (Channel %d, n=%d trials)', channel_to_plot, length(hand_trials)));
grid on;
legend('Mean', 'SEM', 'Location', 'best');

subplot(2, 1, 2);
plot(time_trial, grand_avg_feet(:, channel_to_plot), 'r', 'LineWidth', 2);
hold on;
% Add shaded error region (standard error)
fill([time_trial, fliplr(time_trial)], ...
     [grand_avg_feet(:, channel_to_plot)' + sem_feet(:, channel_to_plot)', ...
      fliplr(grand_avg_feet(:, channel_to_plot)' - sem_feet(:, channel_to_plot)')], ...
     'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold off;
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
title(sprintf('Grand Average - Both Feet (Channel %d, n=%d trials)', channel_to_plot, length(feet_trials)));
grid on;
legend('Mean', 'SEM', 'Location', 'best');

sgtitle(sprintf('Grand Averages - Channel %d', channel_to_plot));

% Overlay comparison plot
figure;
plot(time_trial, grand_avg_hand(:, channel_to_plot), 'b', 'LineWidth', 2);
hold on;
plot(time_trial, grand_avg_feet(:, channel_to_plot), 'r', 'LineWidth', 2);
hold off;
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
title(sprintf('Grand Averages Comparison - Channel %d', channel_to_plot));
legend('Both Hand', 'Both Feet', 'Location', 'best');
grid on;

fprintf('Visualization complete!\n');
fprintf('========================\n');

%% Enhanced Visualization for Processed Signals
fprintf('\n=== Enhanced Visualization ===\n');

% a. Select three meaningful EEG channels
% For motor imagery, select channels over motor cortex areas
% Common choices: C3, Cz, C4 (central channels)
% Assuming standard 10-20 system, we'll select channels 7, 9, 11
% You can modify these based on your actual channel layout
selected_channels = [7, 9, 11];  % Modify based on your channel layout
channel_names = {'C3', 'Cz', 'C4'}; 

fprintf('Selected channels for visualization: %d, %d, %d\n', selected_channels(1), selected_channels(2), selected_channels(3));

% b. Plot the raw and filtered signals for a given trial and for the selected channels
% Use trial 65 (Both Hand - 773)
sample_trial_idx = 65;
fprintf('Visualizing trial #%d (Both Hand)\n', sample_trial_idx);

% Create time vector
time_trial = (0:trial_length-1) / fs;

% Plot raw vs filtered signals - 3 plots (raw, mu, beta), each with 3 channels
figure('Position', [100, 100, 1200, 800]);

% Plot 1: Raw signals for all 3 channels
subplot(3, 1, 1);
hold on;
for i = 1:3
    ch = selected_channels(i);
    plot(time_trial, trials(:, ch, sample_trial_idx), 'LineWidth', 1.5, 'DisplayName', channel_names{i});
end
hold off;
ylabel('Amplitude (\muV)');
title(sprintf('Raw Signal - Trial %d (Both Hand)', sample_trial_idx));
legend('Location', 'best');
grid on;

% Plot 2: μ band filtered signals for all 3 channels
subplot(3, 1, 2);
hold on;
for i = 1:3
    ch = selected_channels(i);
    plot(time_trial, trials_mu(:, ch, sample_trial_idx), 'LineWidth', 1.5, 'DisplayName', channel_names{i});
end
hold off;
ylabel('Amplitude (\muV)');
title(sprintf('μ Band Filtered (10-12 Hz) - Trial %d (Both Hand)', sample_trial_idx));
legend('Location', 'best');
grid on;

% Plot 3: β band filtered signals for all 3 channels
subplot(3, 1, 3);
hold on;
for i = 1:3
    ch = selected_channels(i);
    plot(time_trial, trials_beta(:, ch, sample_trial_idx), 'LineWidth', 1.5, 'DisplayName', channel_names{i});
end
hold off;
ylabel('Amplitude (\muV)');
xlabel('Time (s)');
title(sprintf('β Band Filtered (18-24 Hz) - Trial %d (Both Hand)', sample_trial_idx));
legend('Location', 'best');
grid on;

sgtitle(sprintf('Raw vs Filtered Signals - Trial %d (Both Hand) - All Selected Channels', sample_trial_idx));

% c. Plot the power of the three channels after averaging across all trials for each class
fprintf('Computing averaged power across trials for each class...\n');

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

% Plot averaged power for the three selected channels
figure('Position', [100, 100, 1400, 900]);

for i = 1:3
    ch = selected_channels(i);
    
    % μ band - Both Hand
    subplot(3, 4, (i-1)*4 + 1);
    plot(time_trial, avg_mu_hand(:, ch), 'b', 'LineWidth', 2);
    ylabel('Log Power');
    title(sprintf('%s - μ Band Both Hand', channel_names{i}));
    grid on;
    
    % μ band - Both Feet
    subplot(3, 4, (i-1)*4 + 2);
    plot(time_trial, avg_mu_feet(:, ch), 'r', 'LineWidth', 2);
    ylabel('Log Power');
    title(sprintf('%s - μ Band Both Feet', channel_names{i}));
    grid on;
    
    % β band - Both Hand
    subplot(3, 4, (i-1)*4 + 3);
    plot(time_trial, avg_beta_hand(:, ch), 'b', 'LineWidth', 2);
    ylabel('Log Power');
    title(sprintf('%s - β Band Both Hand', channel_names{i}));
    grid on;
    
    % β band - Both Feet
    subplot(3, 4, (i-1)*4 + 4);
    plot(time_trial, avg_beta_feet(:, ch), 'r', 'LineWidth', 2);
    ylabel('Log Power');
    title(sprintf('%s - β Band Both Feet', channel_names{i}));
    grid on;
end

% Add x-labels to bottom row
for j = 9:12
    subplot(3, 4, j);
    xlabel('Time (s)');
end

sgtitle(sprintf('Averaged Log Band Power Across Trials (n_hand=%d, n_feet=%d)', ...
        length(hand_trials), length(feet_trials)));

% Create comparison plots for each channel (overlay both classes)
figure('Position', [100, 100, 1400, 600]);

for i = 1:3
    ch = selected_channels(i);
    
    % μ band comparison
    subplot(2, 3, i);
    plot(time_trial, avg_mu_hand(:, ch), 'b', 'LineWidth', 2);
    hold on;
    plot(time_trial, avg_mu_feet(:, ch), 'r', 'LineWidth', 2);
    hold off;
    ylabel('Log Power');
    xlabel('Time (s)');
    title(sprintf('%s - μ Band', channel_names{i}));
    legend('Both Hand', 'Both Feet', 'Location', 'best');
    grid on;
    
    % β band comparison
    subplot(2, 3, i + 3);
    plot(time_trial, avg_beta_hand(:, ch), 'b', 'LineWidth', 2);
    hold on;
    plot(time_trial, avg_beta_feet(:, ch), 'r', 'LineWidth', 2);
    hold off;
    ylabel('Log Power');
    xlabel('Time (s)');
    title(sprintf('%s - ', channel_names{i}));
    legend('Both Hand', 'Both Feet', 'Location', 'best');
    grid on;
end

sgtitle('Comparison of Averaged Log Band Power Between Classes');

fprintf('Enhanced visualization complete!\n');
fprintf('===================================\n');

