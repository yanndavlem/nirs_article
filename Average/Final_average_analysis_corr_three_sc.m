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
% Note: Ensure 'trig1', 'trig2', etc. match the trigger names in your raw data
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

%% 3. Analysis: Regression using Best of 3 Specific SC Pairs
% Logic: Restricted search. Find the best correlated SC pair among 
% 3 specific pairs (indices 5-6, 25-26, 43-44).

lc_idx = find(Hb(1).probe.link.ShortSeperation == 0); % Identify Long Channels

% Define the specific Target Short Channels (Indices in Hb.data)
% Pair 1: 5,6 | Pair 2: 25,26 | Pair 3: 43,44
target_sc_indices = [5, 6, 25, 26, 43, 44]; 

% Prepare data: Remove stims to treat as resting state for regression
Hb_nostim = discard_all_stims(Hb);
data_cor = Hb_nostim;

disp('Starting Regression Loop (Best of 3 Target Pairs)...');

for nxn = 1:length(lc_idx)
    current_lc_global_idx = lc_idx(nxn);
    
    % --- STEP A: Find the Best Correlated SC among Targets ---
    all_corrs_this_lc = zeros(length(Hb_nostim), length(target_sc_indices)); % Init matrix
    
    for i = 1:length(Hb_nostim)
        current_lc_data = Hb_nostim(i).data(:, current_lc_global_idx);
        
        % Safety check for indices
        if max(target_sc_indices) > size(Hb_nostim(i).data, 2)
             error('Target SC index (%d) out of bounds.', max(target_sc_indices));
        end
        
        target_sc_data = Hb_nostim(i).data(:, target_sc_indices);
        
        corrs_vector = corr(current_lc_data, target_sc_data); % Calculate correlation
        corrs_vector(isnan(corrs_vector)) = 0;                % Handle NaNs
        all_corrs_this_lc(i, :) = corrs_vector;
    end
    
    avg_abs_corrs = mean(abs(all_corrs_this_lc), 1);     % Avg absolute correlation
    [~, best_sc_local_idx] = max(avg_abs_corrs);         % Find best SC (Index 1-6)
    best_sc_global_idx = target_sc_indices(best_sc_local_idx); % Convert to Global Index
    
    % --- STEP B: Identify the HbO/HbR Pair ---
    if mod(best_sc_global_idx, 2) == 1 
        pair_global_idx = best_sc_global_idx + 1; % Odd index -> Pair is Next
    else 
        pair_global_idx = best_sc_global_idx - 1; % Even index -> Pair is Previous
    end
    
    good_SC = sort([best_sc_global_idx, pair_global_idx]); % Define pair for regression
    
    % --- STEP C: Perform Regression with Specific Pair ---
    Hb_nostim_sc = add_SC_regressors_one_sc(Hb_nostim, good_SC); % Add specific SC regressors from a custom function available on our GITHUB
    
    % Run GLM (OLS)
    job = nirs.modules.GLM();
    job.type = 'OLS';
    job.AddShortSepRegressors = 0; % Manual handling
    sub_stats = job.run(Hb_nostim_sc);
    
    % Subtract the regressed signal from the data
    for i = 1:size(sub_stats, 2)
        probe_tbl = data_cor(i).probe.link;
        
        % Extract SC data and Betas
        data_sc_PCA = data_cor(i).data(:, good_SC);
        stats_tbl = sub_stats(i).table;
        
        % Filter table for current Long Channel
        lc_rows_idx = (stats_tbl.source == probe_tbl.source(lc_idx(nxn)) & ...
                       stats_tbl.detector == probe_tbl.detector(lc_idx(nxn)) & ...
                       strcmp(stats_tbl.type, probe_tbl.type(lc_idx(nxn))));
        a_lc_table = stats_tbl(lc_rows_idx, :);
        
        beta_vec_sc = a_lc_table.beta(startsWith(a_lc_table.cond, 'SC_')); % Identify SC betas
        
        if length(beta_vec_sc) ~= length(good_SC) % Safety Check
            warning('Beta mismatch for LC %d, Subject %d.', lc_idx(nxn), i);
            beta_vec_sc = zeros(length(good_SC), 1);
        end
        
        % Reconstruct and subtract artifact
        beta_matrix_sc = repmat(beta_vec_sc, [1, size(data_sc_PCA, 1)])';
        regressors_sum = sum(data_sc_PCA .* beta_matrix_sc, 2);
        data_cor(i).data(:, lc_idx(nxn)) = data_cor(i).data(:, lc_idx(nxn)) - regressors_sum;
    end
end
disp('Regression Loop Completed.');

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
    r_table_corr_three_sc = export_avg_data_to_R(average_screg);
end

if exist('plots_average', 'file')
    plots_average(average_screg, 'hbo', [-20 25], "HbO Average (Best of 3 Target SCs)");
end

%% 5. Export Results to CSV
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

disp('Exporting data...');
if exist('r_table_corr_three_sc', 'var')
    writetable(r_table_corr_three_sc, fullfile(output_dir, 'r_table_average_corr_three_sc.csv'));
end
disp('Done.');