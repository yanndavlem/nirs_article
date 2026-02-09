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

%% 3. Analysis: GLM with Three Specific SC Pairs
% Logic: Use a custom GLM module that takes a specific list of Short Channels
% (indices 5,6, 25,26, 43,44) and uses them as regressors.

% Define the specific SCs to use
good_SC = [5, 6, 25, 26, 43, 44];
SC_associate_LC = good_SC;


% NOTE: 'nirs.modules.GLM_three_SC' is a CUSTOM module available on our GITHUB. 
% It must be present in the toolbox path nirs-toolbox-master/nirs/modules/GLM_three_SC.m
job = [];
if exist('nirs.modules.GLM_three_SC', 'class')
    job = nirs.modules.GLM_three_SC();
else
    error('Custom module nirs.modules.GLM_three_SC not found.');
end

job.type = 'OLS';
job.AddShortSepRegressors_three_SC = 1; % Activate specific logic
job.whichSC = SC_associate_LC;          % Pass the list of SC indices

sub_stats = job.run(Hb);

% Visualization of Design Matrix (First subject)
dm = nirs.util.getdesign_matrix(Hb(1));
figure; imagesc(dm); title('Design Matrix (Three SC)'); colorbar;


%% 4. Export Results
if ~exist(output_dir, 'dir'); mkdir(output_dir); end
disp('Exporting Results...');

% Save the full stats object
save(fullfile(output_dir, 'stats_trial_block_OLS_ThreeSC.mat'), 'sub_stats');

% Export Individual Subject Stats to CSV
for i = 1:length(sub_stats)
    stats_table = sub_stats(i).table;
    
    % Construct filename: suj_1_Three_SC_NEW.csv, etc.
    filename = sprintf('suj_%d_Three_SC_NEW.csv', i);
    full_path = fullfile(output_dir, filename);
    
    writetable(stats_table, full_path);
end

disp('Done.');