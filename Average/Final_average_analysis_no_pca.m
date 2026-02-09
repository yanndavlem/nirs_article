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

tbl_stims = nirs.createStimulusTable(raw_trim); % Store triggers for later

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

% Resample
job = nirs.modules.Resample();
job.Fs = 2; 
raw = job.run(raw); 

%% 3. Analysis Part A: NO Short Channel Correction

% Block Averaging (Homer2)
job = nirs.modules.Run_HOMER2();
job.fcn = 'hmrBlockAvg';
job.vars.trange = [-5 25];
average_no_sc = job.run(Hb);

% Export & Plot (No SC)
if exist('plots_average', 'file')
    plots_average(average_no_sc, 'hbo', [-20 20], "HbO Average (No SC Correction)");
end

if exist('export_avg_data_to_R', 'file')
    r_table_no_sc = export_avg_data_to_R(average_no_sc);
    if ~exist(output_dir, 'dir'); mkdir(output_dir); end
    writetable(r_table_no_sc, fullfile(output_dir, 'r_table_average_no_sc.csv'));
end


%% 4. Analysis Part B: WITH Short Channel Correction (NO PCA)
% Logic: Use every SC as an individual regressor without dimensionality reduction.

% 1. Prepare Data
Hb_nostim = discard_all_stims(Hb); 

% 2. Add Regressors (Custom function: No PCA, adds all SCs. Function available on our GITHUB)
Hb_nostim_sc = add_SC_regressors_no_PCA(Hb_nostim); 

% 3. Run GLM (OLS)
job = nirs.modules.GLM();
job.type = 'OLS';
job.AddShortSepRegressors = 0; % Manual handling
sub_stats = job.run(Hb_nostim_sc);

% 4. Correction Loop
disp('Applying SC Regression (No PCA) Correction...');
data_cor = Hb_nostim_sc; 

for i = 1:length(data_cor)
    beta_table = sub_stats(i).table;
    probe_tbl = data_cor(i).probe.link;
    lc_indices = find(probe_tbl.ShortSeperation == 0); % Long Channels only
    
    for j = 1:length(lc_indices)
        chan_idx = lc_indices(j);
        
        % Channel Info
        chan_source = probe_tbl.source(chan_idx);
        chan_detector = probe_tbl.detector(chan_idx);
        chan_type = probe_tbl.type{chan_idx};
        
        % Filter betas for this specific channel
        chan_betas = beta_table(beta_table.source == chan_source & ...
                                beta_table.detector == chan_detector & ...
                                strcmp(beta_table.type, chan_type), :);
        
        % Identify SC Regressors (Conditions starting with 'SC_')
        sc_rows = startsWith(chan_betas.cond, 'SC_');
        
        correction_term = 0;
        
        if any(sc_rows)
            % Extract betas and condition names
            sc_betas = chan_betas.beta(sc_rows);
            sc_conds = chan_betas.cond(sc_rows);
            
            % Sum noise from all SC regressors found
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

%% 5. Reconstruct and Average (Corrected Data)

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
    plots_average(average_screg, 'hbo', [-20 20], "HbO Average (All SCs - No PCA)");
end

if exist('export_avg_data_to_R', 'file')
    r_table_no_pca = export_avg_data_to_R(average_screg);
end

%% 6. Export Results to CSV
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

disp('Exporting data...');
if exist('r_table_no_pca', 'var')
    writetable(r_table_no_pca, fullfile(output_dir, 'r_table_average_with_sc_no_pca.csv'));
end
disp('Done.');