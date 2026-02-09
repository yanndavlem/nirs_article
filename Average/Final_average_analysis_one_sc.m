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

% Resample 
 job = nirs.modules.Resample();
 job.Fs = 2; 
 raw = job.run(raw); 

%% 3. Analysis: Regression using Specific Assigned Short Channels
% Logic: Manually assigns specific SCs to specific groups of Long Channels.

lc_idx = find(Hb(1).probe.link.ShortSeperation == 0); % Indices of Long Channels
sc_idx = find(Hb(1).probe.link.ShortSeperation == 1); % Indices of Short Channels

% Prepare Data
Hb_nostim = discard_all_stims(Hb); 
data_cor = Hb_nostim;

disp('Starting Assigned SC Regression Loop...');

for nxn = 1:length(lc_idx)
    current_lc = lc_idx(nxn);
    
    % --- ASSIGNMENT LOGIC ---
    % Define which SC pair to use based on Long Channel Index
    if ismember(nxn, 1:8)
        good_SC = [5, 6];
    elseif ismember(nxn, 9:22)
        good_SC = [25, 26];
    elseif ismember(nxn, 23:48)
        good_SC = [43, 44];
    else
        warning('Long Channel index %d is outside defined ranges. Skipping.', nxn);
        continue;
    end
    
    % --- REGRESSION ---
    
    % Add specific SC regressors (Custom function: add_SC_regressors_one_sc. Function available on our GITHUB)
    Hb_nostim_sc = add_SC_regressors_one_sc(Hb_nostim, good_SC); 
    
    % Run GLM (OLS)
    job = nirs.modules.GLM();
    job.type = 'OLS';
    job.AddShortSepRegressors = 0; 
    sub_stats = job.run(Hb_nostim_sc);
    
    % Apply Correction
    for i = 1:size(sub_stats, 2)
        probe_tbl = data_cor(i).probe.link;
        
        % Extract SC Data (Regressor Signal)
        data_sc_PCA = data_cor(i).data(:, good_SC);
        
        % Extract Betas
        a = sub_stats(i).table;
        
        % Find betas for the current Long Channel (nxn)
        % Note: We look for betas associated with the SC regressors added above.
        % The 'add_SC_regressors_one_sc' function typically names them 'SC_1_PCA', etc.
        % However, since we are doing channel-by-channel subtraction here, we need
        % to match the specific SC regressors.
        
        % Filter table for current LC
        lc_rows = (a.source == probe_tbl.source(current_lc) & ...
                   a.detector == probe_tbl.detector(current_lc) & ...
                   strcmp(a.type, probe_tbl.type(current_lc)));
        
        a_lc = a(lc_rows, :);
        
        % Extract SC betas (Assumes custom function adds conditions starting with 'SC_')
        beta_vec_sc = a_lc.beta(startsWith(a_lc.cond, 'SC_'));
        
        % Safety Check
        if length(beta_vec_sc) ~= length(good_SC)
             % Fallback logic if names differ slightly or if assignment failed
             warning('Beta mismatch for LC %d. Check regressor names.', current_lc);
             beta_vec_sc = zeros(length(good_SC), 1);
        end

        % Reconstruct Noise
        beta_matrix_sc = repmat(beta_vec_sc, [1, size(data_sc_PCA, 1)])';
        regressors_sum = sum(data_sc_PCA .* beta_matrix_sc, 2);
        
        % Subtract Noise
        data_cor(i).data(:, current_lc) = data_cor(i).data(:, current_lc) - regressors_sum;
    end
end
disp('Regression Loop Completed.');

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
    plots_average(average_screg, 'hbo', [-20 20], "HbO Average (Assigned SC Correction)");
end

if exist('export_avg_data_to_R', 'file')
    r_table_one_sc = export_avg_data_to_R(average_screg);
end

%% 5. Export Results to CSV
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

disp('Exporting data...');
if exist('r_table_one_sc', 'var')
    writetable(r_table_one_sc, fullfile(output_dir, 'r_table_average_one_sc.csv'));
end
disp('Done.');