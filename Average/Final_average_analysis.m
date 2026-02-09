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
% Assumes standard folder structure. 
raw = nirs.io.loadDirectory(data_dir, {'subject', 'group'});

%% 2. Preprocessing Pipeline

% Trim Baseline and Adjust Durations
job = nirs.modules.TrimBaseline();
job.preBaseline  = 10;
job.postBaseline = 30;
raw_trim = job.run(raw);

% Standardize stimulus duration (e.g., set all to 10s)
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

%% 3. Analysis WITHOUT Short Channel (SC) Correction

% Block Averaging using Homer2 module
job = nirs.modules.Run_HOMER2();
job.fcn = 'hmrBlockAvg';
job.vars.trange = [-5 25];
average_no_sc = job.run(Hb);

% Plotting (Requires custom function 'plots_average')
if exist('plots_average', 'file')
    plots_average(average_no_sc, 'hbo', [-5 25], "HbO Average (No SC Correction)");
end

% Export to R-ready table (Requires custom function 'export_avg_data_to_R')
if exist('export_avg_data_to_R', 'file')
    r_table_no_sc = export_avg_data_to_R(average_no_sc);
end

%% 4. Analysis WITH Short Channel (SC) Correction (PCA Regression)

% Prepare Data
% Remove stimuli to treat everything as resting state for regression
Hb_nostim = discard_all_stims(Hb); 

% Add PCA of Short Channels as regressors
Hb_nostim_sc = add_SC_regressors_PCA(Hb_nostim); % Function available on our GITHUB

% GLM to regress Short Channels
job = nirs.modules.GLM();
job.type = 'OLS';
job.AddShortSepRegressors = 0; % We manage SC manually via PCA regressors
sub_stats = job.run(Hb_nostim_sc);

% Subtract Short Channel contributions from the signal
data_cor = Hb_nostim_sc;

for i = 1:size(sub_stats, 2)
    % Link table to identify Short (SC) vs Long Channels (LC)
    probe_tbl = data_cor(i).probe.link;
    lc_idx = find(probe_tbl.ShortSeperation == 0);
    
    % Extract PCA SC data and GLM Betas
    data_sc_PCA = extract_sc_PCA(data_cor(i));
    stats_tbl = sub_stats(i).table;
    
    for j = 1:length(lc_idx)
        % Find betas corresponding to the current long channel and SC regressors
        mask = stats_tbl.source == probe_tbl.source(lc_idx(j)) & ...
               stats_tbl.detector == probe_tbl.detector(lc_idx(j)) & ...
               strcmp(stats_tbl.type, probe_tbl.type(lc_idx(j))) & ...
               contains(stats_tbl.cond, 'SC');
           
        beta_vec_sc = stats_tbl.beta(mask);
        
        % Reconstruct the SC component signal
        beta_matrix_sc = repmat(beta_vec_sc, [1, size(data_sc_PCA, 1)])';
        regressors_sum = sum(data_sc_PCA .* beta_matrix_sc, 2);
        
        % Subtract SC component from the Long Channel data
        data_cor(i).data(:, lc_idx(j)) = data_cor(i).data(:, lc_idx(j)) - regressors_sum;
    end
end

% Reconstruct Structure for Averaging
% Remove the PCA regressors now that correction is applied
data_cor_nostim = discard_all_stims(data_cor);

% Reinject the original stimulus timing
job = nirs.modules.ChangeStimulusInfo();
job.ChangeTable = tbl_stims;
data_cor_stims = job.run(data_cor_nostim);

% Block Averaging on Corrected Data
job = nirs.modules.Run_HOMER2();
job.fcn = 'hmrBlockAvg';
job.vars.trange = [-5 25];
average_sc = job.run(data_cor_stims);

% Plotting
if exist('plots_average', 'file')
    plots_average(average_sc, 'hbo', [-20 20], "HbO Average (With SC Correction)");
end

% Export to R table
if exist('export_avg_data_to_R', 'file')
    r_table_sc = export_avg_data_to_R(average_sc);
end

%% 5. Export Results to CSV
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

disp('Exporting data...');
writetable(r_table_sc, fullfile(output_dir, 'r_table_average_with_sc.csv'));
writetable(r_table_no_sc, fullfile(output_dir, 'r_table_average_no_sc.csv'));
disp('Done.');