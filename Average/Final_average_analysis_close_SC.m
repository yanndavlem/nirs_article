%% Setup and Configuration
clc; clear; close all;

% ---------------------------------------------------------
% PATH CONFIGURATION
% ---------------------------------------------------------
% Define your local paths here
toolbox_dir = 'YOUR/PATH/TO/nirs-toolbox';
homer2_dir  = 'YOUR/PATH/TO/homer2';
data_dir    = 'YOUR/PATH/TO/RAW_DATA';   % Folder containing NIRx data folders
output_dir  = 'YOUR/PATH/TO/RESULTS';    % Folder to save CSVs

% Add Toolboxes to Path
addpath(genpath(toolbox_dir));
addpath(genpath(homer2_dir));

%% 1. Data Loading
% Load all NIRx data from the directory automatically
raw = nirs.io.loadDirectory(data_dir, {'subject', 'group'});

%% 2. Preprocessing Pipeline

% Trim Baseline and Adjust Durations
job = nirs.modules.TrimBaseline();
job.preBaseline  = 10;
job.postBaseline = 30;
raw_trim = job.run(raw);

% Standardize stimulus duration (e.g., set all to 10s)
% Note: Ensure 'trig1', 'trig2', etc. match the trigger names in your raw data
raw_trim = nirs.design.change_stimulus_duration(raw_trim, 'trig1', 10);
raw_trim = nirs.design.change_stimulus_duration(raw_trim, 'trig2', 10);
raw_trim = nirs.design.change_stimulus_duration(raw_trim, 'trig3', 10);

% Store stimulus info for later reinjection
tbl_stims = nirs.createStimulusTable(raw_trim);

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

%% 3. Analysis: "Close" Short Channel Regression
% Regression is performed using the SC sharing the same source as the LC

% Identify Channel Indices
lc_idx = find(Hb(1).probe.link.ShortSeperation == 0); % Long Channels
sc_idx = find(Hb(1).probe.link.ShortSeperation == 1); % Short Channels

% Prepare data: Remove stims to treat as resting state for regression
Hb_nostim = discard_all_stims(Hb);
data_cor = Hb_nostim;

disp('Starting Regression Loop (Nearest SC)...');

for nxn = 1:length(lc_idx)
    % 1. Identify the specific Short Channel associated with this Long Channel
    % (Matches Source index and ShortSeperation flag)
    tbl_probe = Hb(1).probe.link;
    good_SC = find(tbl_probe.source == tbl_probe.source(lc_idx(nxn)) & ...
                   strcmp(tbl_probe.type, tbl_probe.type(lc_idx(nxn))) & ...
                   tbl_probe.ShortSeperation == 1);
               
    % 2. Add specific SC regressor
    % Note: 'add_SC_regressors_PCA_close_SC' is a custom function available
    % on our GITHUB
    Hb_nostim_sc = add_SC_regressors_PCA_close_SC(Hb_nostim, good_SC);
    
    % 3. Run GLM (OLS)
    job = nirs.modules.GLM();
    job.type = 'OLS';
    job.AddShortSepRegressors = 0; % Manual handling
    sub_stats = job.run(Hb_nostim_sc);
    
    % 4. Subtract the regressed signal from the data
    for i = 1:size(sub_stats, 2)
        probe_tbl = data_cor(i).probe.link;
        
        % Extract data for the specific SC
        data_sc_PCA = data_cor(i).data(:, good_SC);
        
        % Extract betas
        stats_tbl = sub_stats(i).table;
        beta_vec_sc = stats_tbl.beta(stats_tbl.source == probe_tbl.source(lc_idx(nxn)) & ...
                                     stats_tbl.detector == probe_tbl.detector(lc_idx(nxn)) & ...
                                     strcmp(stats_tbl.type, probe_tbl.type(lc_idx(nxn))));
        
        % Reconstruct artifact signal
        beta_matrix_sc = repmat(beta_vec_sc, [1, size(data_sc_PCA, 1)])';
        regressors_sum = data_sc_PCA .* beta_matrix_sc;
        
        % Subtract artifact from the Long Channel
        data_cor(i).data(:, lc_idx(nxn)) = data_cor(i).data(:, lc_idx(nxn)) - regressors_sum;
    end
    
end
disp('Regression Loop Completed.');

%% 4. Reconstruct and Average

% Remove regressors
data_cor_nostim = discard_all_stims(data_cor);

% Reinject original triggers
job = nirs.modules.ChangeStimulusInfo();
job.ChangeTable = tbl_stims;
data_cor_stims = job.run(data_cor_nostim);

% Block Averaging (Homer2)
job = nirs.modules.Run_HOMER2();
job.fcn = 'hmrBlockAvg';
job.vars.trange = [-5 25];
average_screg = job.run(data_cor_stims);

% Export to R-compatible table
if exist('export_avg_data_to_R', 'file')
    r_table_average_close_sc = export_avg_data_to_R(average_screg);
end

% Plotting results
if exist('plots_average', 'file')
    plots_average(average_screg, 'hbo', [-20 25], "HbO Average (Close SC Correction)");
end

%% 5. Export Results to CSV
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

disp('Exporting data...');
if exist('r_table_average_close_sc', 'var')
    writetable(r_table_average_close_sc, fullfile(output_dir, 'r_table_average_close_sc.csv'));
end
disp('Done.');