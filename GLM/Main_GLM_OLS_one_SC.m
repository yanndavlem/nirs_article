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

%% 3. Analysis: GLM with ROI-Assigned Short Channel
% Logic: Assign a specific SC to each LC based on its ROI (Temporal Left, Occipital, Temporal Right)
% and run a GLM using that specific SC.

lc_idx = find(Hb(1).probe.link.ShortSeperation == 0); % Long Channels
sc_idx = find(Hb(1).probe.link.ShortSeperation == 1); % Short Channels
lc_roi = tbl_roi(tbl_roi.SC == 0,:); % Filter ROI table for LCs only

conditions = {'trig1', 'trig2', 'trig3'};

% Initialize Results Table
num_rows = length(Hb) * height(lc_roi) * length(conditions); 
varNames = {'source', 'detector', 'type', 'cond', 'beta', 'sub'};
varTypes = {'double', 'double', 'string', 'string', 'double', 'double'};
table_results = table('Size', [num_rows, 6], 'VariableTypes', varTypes, 'VariableNames', varNames);

row_counter = 0; 
tbl_probe = Hb(1).probe.link;

disp('Starting GLM Loop (ROI-Assigned SC)...');

for nxn = 1:height(lc_roi)
    
    % --- Determine "Good SC" based on ROI & Type ---
    current_roi = lc_roi.ROI{nxn};
    current_type = lc_roi.Type{nxn}; % Expected format: "'hbo'" or "'hbr'"
    
    good_SC_global_idx = NaN;
    
    % ROI Logic
    if strcmp(current_roi, 'temp_g')
        if contains(current_type, 'hbo')
            good_SC_global_idx = 5;
        else
            good_SC_global_idx = 6;
        end
    elseif strcmp(current_roi, 'occi')
        if contains(current_type, 'hbo')
            good_SC_global_idx = 25;
        else
            good_SC_global_idx = 26;
        end
    elseif strcmp(current_roi, 'temp_d')
        if contains(current_type, 'hbo')
            good_SC_global_idx = 43;
        else
            good_SC_global_idx = 44;
        end
    end
    
    if isnan(good_SC_global_idx)
        warning('Could not determine SC for ROI: %s, Type: %s. Skipping.', current_roi, current_type);
        continue;
    end

    % Get the relative index of this SC among Short Channels
    SC_associate_LC = find(sc_idx == good_SC_global_idx);
    
    if isempty(SC_associate_LC)
        error('Assigned SC index (%d) is not a valid Short Channel in the data.', good_SC_global_idx);
    end

    
    % --- Run Custom GLM ---
    % NOTE: 'nirs.modules.GLM_close_SC' is a CUSTOM module available on our GITHUB. 
    % It must be present in the toolbox path nirs-toolbox-master/nirs/modules/GLM_close_SC.m
    job = [];
    if exist('nirs.modules.GLM_close_SC', 'class')
        job = nirs.modules.GLM_close_SC();
    else
        error('Custom module nirs.modules.GLM_close_SC not found.');
    end
    
    job.type = 'OLS';
    job.AddShortSepRegressors_close_SC = 1; 
    job.whichSC = SC_associate_LC; 
    
    sub_stats = job.run(Hb);
    
    
    % --- Extract Results ---
    % Note: lc_idx(nxn) assumes lc_roi is sorted same as Hb.probe.link LC entries.
    % To be safer, we use source/detector from lc_roi if available, 
    % but here we stick to your logic mapping nxn to lc_idx.
    
    current_lc_idx = lc_idx(nxn); 
    
    for jx = 1:length(sub_stats)
        stats_tbl = sub_stats(jx).table;
        
        for cond = 1:length(conditions)
            ab_idx = find(stats_tbl.source == tbl_probe.source(current_lc_idx) & ...
                          stats_tbl.detector == tbl_probe.detector(current_lc_idx) & ...
                          strcmp(stats_tbl.type, tbl_probe.type(current_lc_idx)) & ...
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

writetable(table_results, fullfile(output_dir, 'beta_table_GLM_one_SC_NEW.csv'));
save(fullfile(output_dir, 'stats_trial_block_OLS_AssignedROI.mat'), 'sub_stats');

disp('Done.');