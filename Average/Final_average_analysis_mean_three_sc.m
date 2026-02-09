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
raw = nirs.io.loadDirectory(data_dir, {'subject', 'group'});

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

% Signal Processing
job = nirs.modules.OpticalDensity();
OD = job.run(raw_trim);

job = nirs.modules.TDDR();
job = eeg.modules.BandPassFilter(job);
job.do_downsample = 0;
job.lowpass  = 0.12;
job.highpass = 0.01;
job = nirs.modules.BeerLambertLaw(job);
Hb = job.run(OD);

% Resample (if needed)
job = nirs.modules.Resample();
job.Fs = 2; 
raw = job.run(raw); 

%% 3. Analysis: Regression using Specific SCs (PCA/Orthogonalized)
% Logic: Uses 'add_SC_regressors_three_sc.m' to add orthogonalized regressors
% from channels [5, 6, 25, 26, 43, 44].

% 1. Prepare Data
Hb_nostim = discard_all_stims(Hb); 

% 2. Add Specific Regressors
%It adds SC_1_PCA to SC_6_PCA. Custom function available on our GITHUB
Hb_nostim_sc = add_SC_regressors_three_sc(Hb_nostim); 

% 3. Run GLM (OLS)
job = nirs.modules.GLM();
job.type = 'OLS';
job.AddShortSepRegressors = 0; 
sub_stats = job.run(Hb_nostim_sc);

% 4. Correction Loop
disp('Applying Orthogonalized SC Regression Correction...');
data_cor = Hb_nostim_sc; 

for i = 1:length(data_cor)
    beta_table = sub_stats(i).table;
    probe_tbl = data_cor(i).probe.link;
    lc_indices = find(probe_tbl.ShortSeperation == 0); % Long Channels only
    
    % Loop through Long Channels
    for j = 1:length(lc_indices)
        chan_idx = lc_indices(j);
        
        chan_source = probe_tbl.source(chan_idx);
        chan_detector = probe_tbl.detector(chan_idx);
        chan_type = probe_tbl.type{chan_idx}; 
        
        % Filter betas for this specific channel
        chan_betas = beta_table(beta_table.source == chan_source & ...
                                beta_table.detector == chan_detector & ...
                                strcmp(beta_table.type, chan_type), :);
        
        % Identify SC Regressors (Any condition starting with 'SC_')
        % This will catch 'SC_1_PCA', 'SC_2_PCA', etc. created by your function.
        sc_rows = startsWith(chan_betas.cond, 'SC_');
        
        correction_term = 0;
        
        if any(sc_rows)
            % Extract betas and condition names for SCs
            sc_betas = chan_betas.beta(sc_rows);
            sc_conds = chan_betas.cond(sc_rows);
            
            % Sum the noise from all 6 orthogonal components
            for k = 1:length(sc_conds)
                reg_vector = data_cor(i).stimulus(sc_conds{k}).vector;
                correction_term = correction_term + (sc_betas(k) * reg_vector);
            end
            
            % Apply Subtraction
            data_cor(i).data(:, chan_idx) = data_cor(i).data(:, chan_idx) - correction_term;
        end
    end
end
disp('Correction Loop Completed.');

%% 4. Reconstruct and Average

data_cor_nostim = discard_all_stims(data_cor); 

% Reinject triggers
job = nirs.modules.ChangeStimulusInfo();
job.ChangeTable = tbl_stims;
data_cor_stims = job.run(data_cor_nostim);

% Block Averaging
job = nirs.modules.Run_HOMER2();
job.fcn = 'hmrBlockAvg';
job.vars.trange = [-5 25];
average_screg = job.run(data_cor_stims);

% Export
if exist('export_avg_data_to_R', 'file')
    r_table_three_sc = export_avg_data_to_R(average_screg);
end

if exist('plots_average', 'file')
    plots_average(average_screg, 'hbo', [-20 25], "HbO Average (Three SC PCA Correction)");
end

%% 5. Export Results
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

disp('Exporting data...');
if exist('r_table_three_sc', 'var')
    writetable(r_table_three_sc, fullfile(output_dir, 'r_table_average_three_sc_pca.csv'));
end
disp('Done.');