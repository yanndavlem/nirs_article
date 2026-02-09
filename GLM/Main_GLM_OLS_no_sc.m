%% Setup and Configuration
clc; clear; close all;

% ---------------------------------------------------------
% PATH CONFIGURATION
% ---------------------------------------------------------
toolbox_dir = 'YOUR/PATH/TO/nirs-toolbox';
data_dir    = 'YOUR/PATH/TO/RAW_DATA';   % Folder containing NIRx data folders
output_dir  = 'YOUR/PATH/TO/RESULTS';    % Folder to save CSVs and MAT files

% Add Toolboxes
addpath(genpath(toolbox_dir));

%% 1. Data Loading
raw = nirs.io.loadDirectory(data_dir, {'subject', 'group'});

%% 2. Preprocessing Pipeline

% Resample
job = nirs.modules.Resample();
job.Fs = 2;
raw_down = job.run(raw);

% Trim Baseline and Adjust Durations
job = nirs.modules.TrimBaseline();
job.preBaseline  = 10;
job.postBaseline = 20;
raw_trim = job.run(raw_down);

% Standardize stimulus duration
raw_trim = nirs.design.change_stimulus_duration(raw_trim, 'trig1', 10);
raw_trim = nirs.design.change_stimulus_duration(raw_trim, 'trig2', 10);
raw_trim = nirs.design.change_stimulus_duration(raw_trim, 'trig3', 10);

% Signal Processing (OD, TDDR, Filter, MBLL)
job = nirs.modules.OpticalDensity();
OD = job.run(raw_trim);

job = nirs.modules.TDDR();               % Motion correction
job = eeg.modules.BandPassFilter(job);   % Filtering
job.do_downsample = 0;
job.lowpass  = 0.12;
job.highpass = 0.01;
job = nirs.modules.BeerLambertLaw(job);  % Convert to HbO/HbR
Hb = job.run(OD);

%% 3. Analysis: GLM OLS (No Short Channel Regression)
% Logic: Standard GLM without using Short Channels as regressors.


job = [];
job = nirs.modules.GLM();
job.type = 'OLS';
job.AddShortSepRegressors = 0; % 0 = Do NOT include Short Channels

sub_stats = job.run(Hb);

% Visualization of Design Matrix (First subject)
dm = nirs.util.getdesign_matrix(Hb(1));
figure; imagesc(dm); title('Design Matrix (No SC)'); colorbar;


%% 4. Export Results
if ~exist(output_dir, 'dir'); mkdir(output_dir); end
disp('Exporting Results...');

% A. Save the full stats object
save(fullfile(output_dir, 'stats_trial_block_OLS_no_sc.mat'), 'sub_stats');

% B. Export Individual Subject Stats to CSV
for i = 1:length(sub_stats)
    stats_table = sub_stats(i).table;
    
    % Construct filename: suj_1_bp001_02_OLS_no_sc.csv, etc.
    filename = sprintf('suj_%d_bp001_02_OLS_no_sc.csv', i);
    full_path = fullfile(output_dir, filename);
    
    writetable(stats_table, full_path);
end

disp('Done.');