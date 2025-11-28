function PSD_files(gdf_filename)
% PSD_files - Load and process GDF file with PSD computation
%
% Usage:
%   PSD_files('filename.gdf')
%
% Input:
%   gdf_filename - Path to the GDF file to process
%
% This function:
%   1. Loads the GDF file
%   2. Applies Laplacian spatial filter
%   3. Computes PSD over time using proc_spectrogram
%   4. Selects meaningful frequency subset (4-48 Hz, step 2 Hz)
%   5. Recomputes EVENT.POS and .DUR for PSD windows
%   6. Saves results to .mat file

    % Check if file exists
    if ~exist(gdf_filename, 'file')
        error('File not found: %s', gdf_filename);
    end
    
    % Load GDF file using BIOSIG toolbox
    fprintf('Loading GDF file: %s\n', gdf_filename);
    [s, h] = sload(gdf_filename);
    
    % Extract metadata
    samplerate = h.SampleRate;
    num_channels = size(s, 2);
    num_samples = size(s, 1);
    
    fprintf('Sampling rate: %.1f Hz\n', samplerate);
    fprintf('Number of channels: %d\n', num_channels);
    fprintf('Number of samples: %d\n', num_samples);
    fprintf('Duration: %.2f seconds\n', num_samples/samplerate);
    
    % Load and apply Laplacian spatial filter
    fprintf('Applying Laplacian spatial filter...\n');
    load('laplacian16.mat', 'lap');
    % Use only first 16 channels for Laplacian filter
    data = s(:, 1:16) * lap;  % Apply Laplacian filter
    
    % Compute PSD over time using proc_spectrogram
    fprintf('Computing PSD over time...\n');
    wlength = 0.5;    % seconds. Length of the external window
    pshift = 0.25;    % seconds. Shift of the internal windows
    wshift = 0.0625;  % seconds. Shift of the external window
    mlength = 1;      % seconds
    
    [PSD, f] = proc_spectrogram(data, wlength, wshift, pshift, samplerate, mlength);
    
    fprintf('PSD size: [%d x %d x %d] (windows x frequencies x channels)\n', ...
            size(PSD, 1), size(PSD, 2), size(PSD, 3));
    fprintf('Frequency range: %.2f - %.2f Hz\n', f(1), f(end));
    
    % Select meaningful frequency subset (4 Hz to 48 Hz, step 2 Hz)
    fprintf('Selecting frequency subset (4-48 Hz, step 2 Hz)...\n');
    target_freqs = 4:2:48;  % 4, 6, 8, 10, ..., 48 Hz
    
    % Find closest indices in frequency grid
    freq_indices = zeros(length(target_freqs), 1);
    selected_freqs = zeros(length(target_freqs), 1);
    
    for i = 1:length(target_freqs)
        [~, idx] = min(abs(f - target_freqs(i)));
        freq_indices(i) = idx;
        selected_freqs(i) = f(idx);
    end
    
    % Extract PSD at selected frequencies
    PSD = PSD(:, freq_indices, :);
    f = selected_freqs;
    
    fprintf('Selected %d frequencies\n', length(f));
    
    % Recompute EVENT.POS and .DUR with respect to PSD windows
    fprintf('Recomputing EVENT positions for PSD windows...\n');
    winconv = 'backward';
    
    events = struct();
    events.TYP = h.EVENT.TYP;
    events.POS = proc_pos2win(h.EVENT.POS, wshift * samplerate, winconv, wlength * samplerate);
    events.DUR = h.EVENT.DUR;
    events.SampleRate = samplerate;
    
    % Adjust DUR if needed (convert from samples to windows)
    if isfield(h.EVENT, 'DUR') && ~isempty(h.EVENT.DUR)
        % DUR is in samples, convert to number of windows
        events.DUR = round(h.EVENT.DUR / (wshift * samplerate));
    end
    
    fprintf('Number of events: %d\n', length(events.TYP));
    
    % Generate output filename
    [filepath, name, ~] = fileparts(gdf_filename);
    output_filename = fullfile(filepath, [name '.mat']);
    
    % Save processed data
    fprintf('Saving processed data to: %s\n', output_filename);
    
    % Store all relevant information
    channels = h.Label;
    
    save(output_filename, 'PSD', 'f', 'events', 'samplerate', 'channels', ...
         'wlength', 'wshift', 'pshift', 'mlength', '-v7.3');
    
    fprintf('Processing complete!\n');
    fprintf('Output saved to: %s\n', output_filename);
    fprintf('PSD dimensions: [%d windows x %d frequencies x %d channels]\n', ...
            size(PSD, 1), size(PSD, 2), size(PSD, 3));
    
end
