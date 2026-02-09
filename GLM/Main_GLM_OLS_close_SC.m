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
raw = nirs.io.loadDirectory(data_dir, {'subject', 'group'}); % Load NIRx

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

%% 3. Analysis: GLM with "Close" Short Channel (Iterative)
% Logic: Loop each LC -> Find closest SC -> Run custom GLM (GLM_close_SC)

lc_idx = find(Hb(1).probe.link.ShortSeperation == 0); % Long Channels
sc_idx = find(Hb(1).probe.link.ShortSeperation == 1); % Short Channels
conditions = {'trig1', 'trig2', 'trig3'};

% Initialize Results Table
num_rows = length(Hb) * length(lc_idx) * length(conditions); 
varNames = {'source', 'detector', 'type', 'cond', 'beta', 'sub'};
varTypes = {'double', 'double', 'string', 'string', 'double', 'double'};
table_close_SC = table('Size', [num_rows, 6], 'VariableTypes', varTypes, 'VariableNames', varNames);

row_counter = 0; 
disp('Starting GLM Loop (One GLM per Long Channel)...');

for nxn = 1:length(lc_idx)
    % 1. Identify current LC info
    tbl_probe = Hb(1).probe.link;
    current_source = tbl_probe.source(lc_idx(nxn));
    current_type   = tbl_probe.type(lc_idx(nxn));
    
    % 2. Find closest SC (Same Source, Same Type)
    good_SC = find(tbl_probe.source == current_source & ...
                   strcmp(tbl_probe.type, current_type) & ...
                   tbl_probe.ShortSeperation == 1);
               
    SC_associate_LC = find(sc_idx == good_SC); % Get SC index
    
    if isempty(SC_associate_LC)
        warning('No SC found for LC index %d. Skipping.', lc_idx(nxn)); continue;
    end
    
    % 3. Run Custom GLM Module (Must be in toolbox path!)
    % NOTE: 'nirs.modules.GLM_close_SC' is a CUSTOM module available on our GITHUB. 
    % It must be present in the toolbox path nirs-toolbox-master/nirs/modules/GLM.m
    if exist('nirs.modules.GLM_close_SC', 'class')
        job = nirs.modules.GLM_close_SC();
    else
        error('Custom module nirs.modules.GLM_close_SC not found.');
    end
    
    job.type = 'OLS';
    job.AddShortSepRegressors_close_SC = 1; 
    job.whichSC = SC_associate_LC; 
    
    sub_stats = job.run(Hb);
    
    % 4. Extract Results for this specific LC
    for jx = 1:length(sub_stats) 
        stats_tbl = sub_stats(jx).table;
        
        for cond = 1:length(conditions)
            ab_idx = find(stats_tbl.source == tbl_probe.source(lc_idx(nxn)) & ...
                          stats_tbl.detector == tbl_probe.detector(lc_idx(nxn)) & ...
                          strcmp(stats_tbl.type, tbl_probe.type(lc_idx(nxn))) & ...
                          strcmp(stats_tbl.cond, conditions{cond}));
            
            if ~isempty(ab_idx)
                row_counter = row_counter + 1;
                table_close_SC.source(row_counter)   = stats_tbl.source(ab_idx);
                table_close_SC.detector(row_counter) = stats_tbl.detector(ab_idx);
                table_close_SC.type(row_counter)     = stats_tbl.type(ab_idx);
                table_close_SC.cond(row_counter)     = stats_tbl.cond(ab_idx);
                table_close_SC.beta(row_counter)     = stats_tbl.beta(ab_idx);
                table_close_SC.sub(row_counter)      = jx;
            end
        end
    end
end
disp('GLM Loop Completed.');

%% 4. Export Results
if ~exist(output_dir, 'dir'); mkdir(output_dir); end
disp('Exporting Results...');

% Save table and last run stats
writetable(table_close_SC, fullfile(output_dir, 'beta_table_GLM_close_SC_NEW.csv'));
save(fullfile(output_dir, 'stats_trial_block_OLS_last_run.mat'), 'sub_stats');

disp('Done.');