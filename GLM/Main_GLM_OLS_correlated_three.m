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

%% 3. Analysis: GLM with Best Correlated SC (from 3 Target Pairs)
% Logic: Check correlation with 3 specific pairs -> Select Best Pair -> Run Custom GLM

lc_idx = find(Hb(1).probe.link.ShortSeperation == 0); % Long Channels
conditions = {'trig1', 'trig2', 'trig3'};

% Define Target SC Pairs
target_sc_pairs = {[5, 6], [25, 26], [43, 44]};
target_sc_indices = [target_sc_pairs{:}]; % Flattened: [5, 6, 25, 26, 43, 44]

% Initialize Results Table
num_rows = length(Hb) * length(lc_idx) * length(conditions); 
varNames = {'source', 'detector', 'type', 'cond', 'beta', 'sub'};
varTypes = {'double', 'double', 'string', 'string', 'double', 'double'};
table_results = table('Size', [num_rows, 6], 'VariableTypes', varTypes, 'VariableNames', varNames);

row_counter = 0; 
disp('Starting GLM Loop (Best of 3 Target Pairs)...');

for nxn = 1:length(lc_idx)
    current_lc_global_idx = lc_idx(nxn);
    
    % --- STEP A: Find Best Correlated Pair among Targets ---
    all_corrs = zeros(length(Hb), length(target_sc_indices));
    
    for i = 1:length(Hb)
        lc_data = Hb(i).data(:, current_lc_global_idx);
        
        if max(target_sc_indices) > size(Hb(i).data, 2)
            error('Target SC index out of bounds for Subject %d', i);
        end
        sc_data = Hb(i).data(:, target_sc_indices);
        
        c = corr(lc_data, sc_data);
        c(isnan(c)) = 0;
        all_corrs(i, :) = c;
    end
    
    % Average Absolute Correlation across subjects
    avg_abs_corrs = mean(abs(all_corrs), 1);
    
    % Find best SC index (1 to 6)
    [~, best_sc_local_idx] = max(avg_abs_corrs);
    
    % Map local index (1-6) to Pair Index (1, 2, or 3)
    pair_index = ceil(best_sc_local_idx / 2);
    
    % Get the actual indices for the selected pair
    SC_associate_LC = target_sc_pairs{pair_index};
    
    
    % --- STEP B: Run Custom GLM ---
    % NOTE: 'nirs.modules.GLM_corr_three_SC' is a CUSTOM module available on our GITHUB. 
    % It must be present in the toolbox path nirs-toolbox-master/nirs/modules/GLM.m
    job = [];
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
end

% Cleanup empty rows
table_results = table_results(1:row_counter, :);
disp('GLM Loop Completed.');

%% 4. Export Results
if ~exist(output_dir, 'dir'); mkdir(output_dir); end
disp('Exporting Results...');

writetable(table_results, fullfile(output_dir, 'beta_table_GLM_MostCorrelated_three_SC_NEW.csv'));
save(fullfile(output_dir, 'stats_trial_block_OLS_MostCorrelated_Three.mat'), 'sub_stats');

disp('Done.');