%% Setup and Configuration
clc; clear; close all;

% ---------------------------------------------------------
% PATH CONFIGURATION
% ---------------------------------------------------------
toolbox_dir = 'YOUR/PATH/TO/nirs-toolbox';
data_dir    = 'YOUR/PATH/TO/RAW_DATA';   % Folder containing NIRx data folders
output_dir  = 'YOUR/PATH/TO/RESULTS';    % Folder to save CSVs and MAT files

% Add Toolboxes to Path
addpath(genpath(toolbox_dir));

%% 1. Data Loading
raw = nirs.io.loadDirectory(data_dir, {'subject', 'group'}); % Load all NIRx data

%% 2. Preprocessing Pipeline

% Resample
job = nirs.modules.Resample();
job.Fs = 2;
raw_down = job.run(raw);

% Trim Baseline and Adjust Durations
job = nirs.modules.TrimBaseline();
job.preBaseline  = 10;
job.postBaseline = 20; % Note: Changed to 20s as per your script (was 30s in previous ones)
raw_trim = job.run(raw);

% Standardize stimulus duration
raw_trim = nirs.design.change_stimulus_duration(raw_trim, 'trig1', 10);
raw_trim = nirs.design.change_stimulus_duration(raw_trim, 'trig2', 10);
raw_trim = nirs.design.change_stimulus_duration(raw_trim, 'trig3', 10);

% Quality Control (Bad Channel Detection)
% Detect bad channels based on SCI threshold (0.6)
bad_chan_table = bad_channels_indices(raw_trim, 0.6);

% Signal Processing (OD, TDDR, Filter, MBLL)
job = nirs.modules.OpticalDensity();
OD = job.run(raw_trim);

job = nirs.modules.TDDR();               % Motion correction
job = eeg.modules.BandPassFilter(job);   % Filtering
job.lowpass  = 0.12;
job.highpass = 0.01;
job = nirs.modules.BeerLambertLaw(job);  % Convert to HbO/HbR
Hb = job.run(OD);

%% 3. GLM Analysis (OLS with Built-in SC Regression)

job = nirs.modules.GLM();
job.type = 'OLS';               % Ordinary Least Squares
job.AddShortSepRegressors = 1;  % 1 = Automatically include Short Channels as regressors

% Run GLM
sub_stats = job.run(Hb);

% Visualization of Design Matrix (First subject)
dm = nirs.util.getdesign_matrix(Hb(1));
figure;
imagesc(dm);
title('Design Matrix (GLM OLS)');
colorbar;

%% 4. Export Results

% Create output directory if it doesn't exist
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

disp('Exporting Results...');

%  Export Individual Subject Stats to CSV (for R analysis)
% Loop through each subject and save their stats table separately
for i = 1:numel(sub_stats)
    stats_table = sub_stats(i).table;
    
    % Construct filename: suj_1_all_sc_v2.csv, suj_2..., etc.
    filename = sprintf('suj_%d_all_sc_v2.csv', i);
    full_path = fullfile(output_dir, filename);
    
    writetable(stats_table, full_path);
end

disp('Done.');