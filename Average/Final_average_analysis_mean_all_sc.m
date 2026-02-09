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
% Load all NIRx data from the directory automatically
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

tbl_stims = nirs.createStimulusTable(raw_trim); % Store stimulus info for later reinjection

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

% Resample (Optional - verify if needed for your specific dataset)
% job = nirs.modules.Resample();
% job.Fs = 2; 
% Hb = job.run(Hb); 

%% 3. Analysis: Regression using Mean of ALL Short Channels
% Logic: Calculate mean SC (HbO) and mean SC (HbR), regress them out of Long Channels.

% 1. Prepare Data: Remove triggers to treat as resting state
Hb_nostim = discard_all_stims(Hb); 

% 2. Add Mean SC Regressors (Custom function required available on our GITHUB)
% Adds 'SC_mean_hbo' and 'SC_mean_hbr' as stimulus vectors
Hb_nostim_sc = add_SC_regressors_mean_all_sc(Hb_nostim); 

% 3. Run GLM (OLS) to get betas for the mean SCs
job = nirs.modules.GLM();
job.type = 'OLS';
job.AddShortSepRegressors = 0; % Manual handling
sub_stats = job.run(Hb_nostim_sc);

% 4. Correction Loop: Subtract estimated noise
disp('Applying Mean SC Regression Correction...');
data_cor = Hb_nostim_sc; 

for i = 1:length(data_cor)
    % Extract Beta Table and Probe Info
    beta_table = sub_stats(i).table;
    probe_tbl = data_cor(i).probe.link;
    lc_indices = find(probe_tbl.ShortSeperation == 0); % Long Channels only

    % Extract the Regressor Vectors created in Step 2
    regressor_hbo = data_cor(i).stimulus('SC_mean_hbo').vector;
    regressor_hbr = data_cor(i).stimulus('SC_mean_hbr').vector;
    
    % Loop through Long Channels to apply correction
    for j = 1:length(lc_indices)
        chan_idx = lc_indices(j);
        
        % Channel Props
        chan_source = probe_tbl.source(chan_idx);
        chan_detector = probe_tbl.detector(chan_idx);
        chan_type = probe_tbl.type{chan_idx}; % 'hbo' or 'hbr'
        
        correction_term = 0;
        
        % Identify correct regressor and beta based on type (HbO vs HbR)
        if strcmp(chan_type, 'hbo')
            target_cond = 'SC_mean_hbo';
            target_regressor = regressor_hbo;
        elseif strcmp(chan_type, 'hbr')
            target_cond = 'SC_mean_hbr';
            target_regressor = regressor_hbr;
        else
            continue; % Skip if unknown type
        end
        
        % Find the specific beta value
        beta_mask = beta_table.source == chan_source & ...
                    beta_table.detector == chan_detector & ...
                    strcmp(beta_table.type, chan_type) & ...
                    strcmp(beta_table.cond, target_cond);
        
        if any(beta_mask)
            beta_val = beta_table.beta(beta_mask);
            % Calculate Noise: Beta * Regressor Signal
            correction_term = beta_val * target_regressor;
            
            % Apply Subtraction
            data_cor(i).data(:, chan_idx) = data_cor(i).data(:, chan_idx) - correction_term;
        end
    end
end
disp('Correction Loop Completed.');

%% 4. Reconstruct and Average

data_cor_nostim = discard_all_stims(data_cor); % Remove regressors

% Reinject original triggers
job = nirs.modules.ChangeStimulusInfo();
job.ChangeTable = tbl_stims;
data_cor_stims = job.run(data_cor_nostim);

% Block Averaging (Homer2)
job = nirs.modules.Run_HOMER2();
job.fcn = 'hmrBlockAvg';
job.vars.trange = [-5 25];
average_screg = job.run(data_cor_stims);

% Export and Plot
if exist('export_avg_data_to_R', 'file')
    r_table_mean_sc = export_avg_data_to_R(average_screg);
end

if exist('plots_average', 'file')
    plots_average(average_screg, 'hbo', [-20 25], "HbO Average (Mean SC Correction)");
end

%% 5. Export Results to CSV
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

disp('Exporting data...');
if exist('r_table_mean_sc', 'var')
    writetable(r_table_mean_sc, fullfile(output_dir, 'r_table_average_with_mean_all_sc.csv'));
end
disp('Done.');