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

%% 3. Analysis: GLM with Most Correlated SC Pair
% Logic: Loop LC -> Find Max Corr SC (across subjects) -> Select Pair -> Run Custom GLM

lc_idx = find(Hb(1).probe.link.ShortSeperation == 0); % Long Channels
sc_idx = find(Hb(1).probe.link.ShortSeperation == 1); % Short Channels
conditions = {'trig1', 'trig2', 'trig3'};

% Initialize Results Table
num_rows = length(Hb) * length(lc_idx) * length(conditions); 
varNames = {'source', 'detector', 'type', 'cond', 'beta', 'sub'};
varTypes = {'double', 'double', 'string', 'string', 'double', 'double'};
table_results = table('Size', [num_rows, 6], 'VariableTypes', varTypes, 'VariableNames', varNames);

row_counter = 0; 
disp('Starting GLM Loop (Most Correlated SC Pair)...');

for nxn = 1:length(lc_idx)
    current_lc_global_idx = lc_idx(nxn);
    
    % --- STEP A: Find Best Correlated SC (HbO+HbR Logic) ---
    all_corrs = zeros(length(Hb), length(sc_idx));
    
    for i = 1:length(Hb)
        lc_data = Hb(i).data(:, current_lc_global_idx);
        sc_data = Hb(i).data(:, sc_idx);
        
        c = corr(lc_data, sc_data);
        c(isnan(c)) = 0;
        all_corrs(i, :) = c;
    end
    
    % Average Absolute Correlation
    [~, best_sc_local_idx] = max(mean(abs(all_corrs), 1));
    best_sc_global_idx = sc_idx(best_sc_local_idx);
    
    % Find Pair (Even/Odd Logic)
    if mod(best_sc_global_idx, 2) == 1 
        pair_global_idx = best_sc_global_idx + 1; % Odd -> Next
    else 
        pair_global_idx = best_sc_global_idx - 1; % Even -> Prev
    end
    SC_associate_LC = sort([best_sc_global_idx, pair_global_idx]);
    
    
    % --- STEP B: Run Custom GLM ---
    % NOTE: 'nirs.modules.GLM_close_SC' is a CUSTOM module available on our GITHUB. 
    % It must be present in the toolbox path nirs-toolbox-master/nirs/modules/GLM.m
    if exist('nirs.modules.GLM_corr_three_SC', 'class')
        job = nirs.modules.GLM_corr_three_SC();
    else
        error('Custom module nirs.modules.GLM_corr_three_SC not found.');
    end
    
    job.type = 'OLS';
    job.AddShortSepRegressors_corr_three_SC = 1; 
    job.whichSC = SC_associate_LC; 
    
    sub_stats = job.run(Hb);
    
    
    % --- STEP C: Extract Results for Current LC ---
    tbl_probe = Hb(1).probe.link;
    
    for jx = 1:length(sub_stats)
        stats_tbl = sub_stats(jx).table;
        
        for cond = 1:length(conditions)
            ab_idx = find(stats_tbl.source == tbl_probe.source(current_lc_global_idx) & ...
                          stats_tbl.detector == tbl_probe.detector(current_lc_global_idx) & ...
                          strcmp(stats_tbl.type, tbl_probe.type(current_lc_global_idx)) & ...
                          strcmp(stats_tbl.cond, conditions{cond}));
            
            if ~isempty(ab_idx)
                row_counter = row_counter + 1;
                table_results.source(row_counter)   = stats_tbl.source(ab_idx);
                table_results.detector(row_counter) = stats_tbl.detector(ab_idx);
                table_results.type(row_counter)     = stats_tbl.type(ab_idx);
                table_results.cond(row_counter)     = stats_tbl.cond(ab_idx);
                table_results.beta(row_counter)     = stats_tbl.beta(ab_idx);
                table_results.sub(row_counter)      = jx;
            end
        end
    end
    % fprintf('Processed LC %d/%d\n', nxn, length(lc_idx));
end

% Cleanup empty rows
table_results = table_results(1:row_counter, :);

%% 4. Export Results
if ~exist(output_dir, 'dir'); mkdir(output_dir); end
disp('Exporting Results...');

writetable(table_results, fullfile(output_dir, 'beta_table_GLM_MostCorrelated_SC_NEW.csv'));
save(fullfile(output_dir, 'stats_trial_block_OLS_MostCorrelated.mat'), 'sub_stats');

disp('Done.');