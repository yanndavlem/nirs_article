%% Setup and Configuration
clc; clear; close all;

% ---------------------------------------------------------
% PATH CONFIGURATION
% ---------------------------------------------------------
toolbox_dir = 'YOUR/PATH/TO/nirs-toolbox';
homer2_dir  = 'YOUR/PATH/TO/homer2';
data_dir    = 'YOUR/PATH/TO/RAW_DATA';   % Folder containing NIRx data folders
output_dir  = 'YOUR/PATH/TO/RESULTS';    % Folder to save CSVs

% Add Toolboxes to Path
addpath(genpath(toolbox_dir));
addpath(genpath(homer2_dir));

%% 1. Data Loading
raw = nirs.io.loadDirectory(data_dir, {'subject', 'group'}); % Load all NIRx data

%% 2. Preprocessing Pipeline

job = nirs.modules.TrimBaseline();
job.preBaseline  = 10;
job.postBaseline = 30;
raw_trim = job.run(raw);

% Standardize stimulus duration
raw_trim = nirs.design.change_stimulus_duration(raw_trim, 'trig1', 10);
raw_trim = nirs.design.change_stimulus_duration(raw_trim, 'trig2', 10);
raw_trim = nirs.design.change_stimulus_duration(raw_trim, 'trig3', 10);

tbl_stims = nirs.createStimulusTable(raw_trim); % Store triggers

% Signal Processing (OD, TDDR, Filter, MBLL)
job = nirs.modules.OpticalDensity();
OD = job.run(raw_trim);

job = nirs.modules.TDDR();
job = eeg.modules.BandPassFilter(job);
job.do_downsample = 0;
job.lowpass  = 0.12;
job.highpass = 0.01;
job = nirs.modules.BeerLambertLaw(job);
Hb = job.run(OD);

%% 3. Analysis: Regression using "Three SC" Strategy
% Logic: Uses custom functions to add regressors based on 3 specific SCs.

% 1. Prepare Data
Hb_nostim = discard_all_stims(Hb); 

% 2. Add Specific Regressors (Custom function available on our GITHUB)
Hb_nostim_sc = add_SC_regressors_three_sc(Hb_nostim); 

% 3. Run GLM (OLS)
job = nirs.modules.GLM();
job.type = 'OLS';
job.AddShortSepRegressors = 0; % Manual handling
sub_stats = job.run(Hb_nostim_sc);

% 4. Correction Loop
disp('Applying Three SC Regression Correction...');
data_cor = Hb_nostim_sc; 

for i = 1:length(data_cor)
    probe_tbl = data_cor(i).probe.link;
    lc_indices = find(probe_tbl.ShortSeperation == 0); % Long Channels
    
    % Extract the regressor data using your custom function
    % (Should return a matrix of size Time x NumRegressors)
    data_sc_PCA = recuperer_regresseurs_sc(data_cor(i));
    
    % Extract Betas
    a = sub_stats(i).table;
    
    for j = 1:length(lc_indices)
        chan_idx = lc_indices(j);
        
        % Filter betas for this specific channel and SC regressors
        % Note: 'contains(..., 'SC')' catches all regressors with 'SC' in name
        beta_vec_sc = a.beta(a.source == probe_tbl.source(chan_idx) & ...
                             a.detector == probe_tbl.detector(chan_idx) & ...
                             strcmp(a.type, probe_tbl.type(chan_idx)) & ...
                             contains(a.cond, 'SC'));
        
        % Safety Check: Ensure dimensions match
        if length(beta_vec_sc) ~= size(data_sc_PCA, 2)
            % This handles cases where betas might be missing or mismatched
            warning('Dimension mismatch: Betas (%d) vs Regressors (%d) for LC %d', ...
                    length(beta_vec_sc), size(data_sc_PCA, 2), chan_idx);
            continue; 
        end
        
        % Reconstruct Noise (Matrix Multiplication)
        % Dimensions: [Time x Regressors] * [Regressors x 1] = [Time x 1]
        % Note: I adjusted the calculation to be cleaner than repmat+sum
        regressors_sum = data_sc_PCA * beta_vec_sc;
        
        % Subtract Noise
        data_cor(i).data(:, chan_idx) = data_cor(i).data(:, chan_idx) - regressors_sum;
    end
end
disp('Correction Loop Completed.');

%% 4. Reconstruct and Average

data_cor_nostim = discard_all_stims(data_cor); % Remove regressors

% Reinject triggers
job = nirs.modules.ChangeStimulusInfo();
job.ChangeTable = tbl_stims;
data_cor_stims = job.run(data_cor_nostim);

% Block Averaging
job = nirs.modules.Run_HOMER2();
job.fcn = 'hmrBlockAvg';
job.vars.trange = [-5 25];
average_screg = job.run(data_cor_stims);

% Export and Plot
if exist('plots_average', 'file')
    plots_average(average_screg, 'hbo', [-20 20], "HbO Average (Three SC Correction)");
end

if exist('export_avg_data_to_R', 'file')
    r_table_three_sc = export_avg_data_to_R(average_screg);
end

%% 5. Export Results to CSV
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

disp('Exporting data...');
if exist('r_table_three_sc', 'var')
    writetable(r_table_three_sc, fullfile(output_dir, 'r_table_average_three_sc_v2.csv'));
end
disp('Done.');