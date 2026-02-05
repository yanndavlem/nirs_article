% L'export des données vers un tableau utilisable par R long. Si ce n'est
% pas utile, il faut commenter %"r_table_average_no_sc" et
% r_table_average_with_sc


clc;clear;close all
cd('H:/Autres ordinateurs/Mon ordinateur/Yann/Data Analyse/1. Analyse NIRS Toolbox/Analyse_TPP/Analyse CRNL/Scripts_Matlab/')
% ajoute au path MATLAB la NIRSToolbox (mettre l'accès à la place des trois petits
% points
dir_nirtoolbox = 'C:/fNIRS AnalyzIR/nirs-toolbox-master/';
dir_hmr2 = 'H:/Autres ordinateurs/Mon ordinateur/Yann/Mat�riels utiles/Toolbox/homer2/';
dir_toolbox = 'H:/Autres ordinateurs/Mon ordinateur/Yann/Mat�riels utiles/Toolbox/';
addpath(genpath(dir_hmr2));
addpath(genpath(dir_toolbox));
addpath(genpath(dir_nirtoolbox));


% importation des donnees brutes
folders = dir('H:/Autres ordinateurs/Mon ordinateur/Yann/Data Analyse/1. Analyse NIRS Toolbox/Analyse_TPP/Analyse CRNL/fichiers_nirx/*2024*'); % faire un dossier avec tous tes dossiers NIRX et mettre dans les trois petits points le chemin vers ce dossier. La fonction dir va lister tous les noms de dossier
for i = 1:length(folders)
    raw(i,1) =  nirs.io.loadNIRx(['H:/Autres ordinateurs/Mon ordinateur/Yann/Data Analyse/1. Analyse NIRS Toolbox/Analyse_TPP/Analyse CRNL/fichiers_nirx/' folders(i).name],true); % importation, mettre dans les trois petits points le chemin vers le dossier contenat les données brutes
end
type_stim_par_sujet = readtable('H:/Autres ordinateurs/Mon ordinateur/Yann/Data Analyse/1. Analyse NIRS Toolbox/Analyse_TPP/Analyse CRNL/Script_R/type_stim_par_sujet.csv');

% Renomme les triggers
j = [];
j = nirs.modules.RenameStims();
j.listOfChanges = { ...
    'channel_1' 'trig1'
    'channel_2' 'trig2'
    'channel_4' 'trig3'
    };
raw = j.run(raw);



%% Check et corrige les triggers

%Récupère les triggers manquants
tbl = nirs.createStimulusTable(raw); %Recupere les triggers

%Vérifie les liens où les triggers sont absents
cell_trig = {};
trig_ok = [];
for i = 1:size(tbl,1) %Additionne et verifie le nombre de trigger
    num_trig_1 = length(tbl.trig1(i).onset);
    num_trig_2 = length(tbl.trig2(i).onset);
    num_trig_3 = length(tbl.trig3(i).onset);
    num_total_trig = num_trig_1+num_trig_2+num_trig_3;
    
    cell_trig(i) = {[num_trig_1 num_trig_2 num_trig_3]}; %met 1 si tous les triggers / 0 s'il manque 1 trigger
    if num_total_trig <45
        trig_ok = [trig_ok 1];
    else
        trig_ok = [trig_ok 0];
    end
end

%Ajoute les triggers manquant grâce au fichier datalog txt qui a été préalablement
%converti en fichier csv
for j = find(trig_ok)
    time_trig_1 = [tbl.trig1(j).onset repmat(1,[length(tbl.trig1(j).onset) 1])];
    time_trig_2 = [tbl.trig2(j).onset repmat(2,[length(tbl.trig2(j).onset) 1])];
    time_trig_3 = [tbl.trig3(j).onset repmat(3,[length(tbl.trig3(j).onset) 1])];
    
    time_trig_all = sortrows([time_trig_1;time_trig_2;time_trig_3],1);
    trig_miss = find(diff(time_trig_all(:,1)) >30);
    true_trig_miss = [];
    
    for k =1:length(trig_miss)
        %true_trig_miss = [true_trig_miss trig_miss(k)+k];
        true_trig_miss = trig_miss(k)+k;
        tbl_good_suj = type_stim_par_sujet(type_stim_par_sujet.Sujet == j,:);
        new_trig_onset = time_trig_all(trig_miss(k),1) + 10 + tbl_good_suj.Duree_Baseline(true_trig_miss);
        trig_name = tbl_good_suj.triggers(true_trig_miss);
        
        raw(j).stimulus(trig_name).onset = [raw(j).stimulus(trig_name).onset; new_trig_onset];
        raw(j).stimulus(trig_name).dur = repmat(1,length(raw(j).stimulus(trig_name).onset),1);
        raw(j).stimulus(trig_name).amp = repmat(1,length(raw(j).stimulus(trig_name).onset),1);
    end
    
    %Ajoute une ligne spécial pour le 7e sujet, il manque le dernier trigger
    if j == 7
        new_trig_special_7 = time_trig_all(end,1) + 10 + 10;
        raw(j).stimulus('trig1').onset = [raw(j).stimulus('trig1').onset; new_trig_special_7];
        raw(j).stimulus('trig1').dur = repmat(1,length(raw(j).stimulus('trig1').onset),1);
        raw(j).stimulus('trig1').amp = repmat(1,length(raw(j).stimulus('trig1').onset),1);
    end
end

%% Recupere le signal des 15 premiers stims
% tbl = nirs.createStimulusTable(raw); %Recupere les triggers
%
% for m = 1:size(raw,1)
%     time_trig_1 = [tbl.trig1(m).onset repmat(1,[length(tbl.trig1(m).onset) 1])];
%     time_trig_2 = [tbl.trig2(m).onset repmat(2,[length(tbl.trig2(m).onset) 1])];
%     time_trig_3 = [tbl.trig3(m).onset repmat(3,[length(tbl.trig3(m).onset) 1])];
%
%     time_trig_all = sortrows([time_trig_1;time_trig_2;time_trig_3],1);
%     time_trig_15 = time_trig_all(1:15,:);
%
%     for n = 1:3
%         raw(m).stimulus(['trig' num2str(n)]).onset = time_trig_15(time_trig_15(:,2)== n ,1);
%         raw(m).stimulus(['trig' num2str(n)]).dur = repmat(1,length(raw(m).stimulus(['trig' num2str(n)]).onset),1);
%         raw(m).stimulus(['trig' num2str(n)]).amp = repmat(1,length(raw(m).stimulus(['trig' num2str(n)]).onset),1);
%     end
% end
%% Preprocess
% Decoupe le signal pour ne garder que le signal intéressant
j=[];
j = nirs.modules.TrimBaseline( j );
j.preBaseline  = 10;
j.postBaseline = 30;
raw_trim = j.run(raw);

%change stimulus durations
raw_trim = nirs.design.change_stimulus_duration(raw_trim,'trig1',10);
raw_trim = nirs.design.change_stimulus_duration(raw_trim,'trig2',10);
raw_trim = nirs.design.change_stimulus_duration(raw_trim,'trig3',10);

% bad_chan = bad_channels_indices(raw_trim);

%pretraitement des données par la NIRS Toolbox

% Stocke pour plus tard la table de triggers à réinjecter à la fin
tbl_stims = nirs.createStimulusTable(raw_trim); %Recupere les triggers


% % Resample
% j = [];
% j = nirs.modules.Resample(j);
% j.Fs = 2;
% raw_trim = j.run(raw_trim);


%pretraitement des données par la NIRS Toolbox

jobs = [];
jobs = nirs.modules.OpticalDensity(jobs);
OD = jobs.run(raw_trim); % transforme tes données en densité optiques

jobs = [];
jobs = nirs.modules.TDDR(); % ajoute aux "jobs" la correction TDDR
jobs = eeg.modules.BandPassFilter(jobs); % ajoute aux jobs un filtre bandpass
jobs.do_downsample = 0;
jobs.lowpass = 0.12; % choisis ton passe-bas
jobs.highpass = 0.01; % choisis ton passe-haut
jobs = nirs.modules.BeerLambertLaw(jobs); % ajoute aux jobs la transformation de BeerLambert
Hb = jobs.run(OD); % applique tous les jobs précédent aux densités optiques, tu obtiens des données d'hémoglobine oxy/desoxy



%% Fait l'analyse pour ne prendre que la r�gression d'un short channel

lc_idx = find(Hb(1).probe.link.ShortSeperation == 0); % recupere les numeros de canaux des long channels
sc_idx = find(Hb(1).probe.link.ShortSeperation == 1); % recupere les numeros de canaux des shorts channels
conditions = {'trig1','trig2','trig3'};
sz = [size(Hb,1)*length(lc_idx)*3,6];
varNames = {'source','detector','type','cond','beta','sub'};
varTypes = {'double','double','string','string','double','double'};
table_close_SC = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
flag = 0;

% Gets rid of all stims
Hb_nostim = discard_all_stims(Hb);
data_cor = Hb_nostim;

for nxn = 1:length(lc_idx)
    disp(['long channel ' num2str(lc_idx(nxn)) ' start'])
    
    %Simplifier la lecture
    tbl_probe = Hb(1).probe.link;
    %Recupere les SC par rapport aux longs channels correspondant
    good_SC = find(tbl_probe.source == tbl_probe.source(lc_idx(nxn)) & strcmp(tbl_probe.type,tbl_probe.type(lc_idx(nxn)))& tbl_probe.ShortSeperation == 1);
    %Attribue le SC associe au long channel
    SC_associate_LC = find(sc_idx == good_SC);
    
    
    % Adds the PCA of SC as regressors
    Hb_nostim_sc = add_SC_regressors_PCA_close_SC(Hb_nostim,good_SC);
    
    % GLM
    job = [];
    job = nirs.modules.GLM;
    job.type = 'OLS'; % Choisi le type de GLM
    job.AddShortSepRegressors= 0; % D�cide de prendre ou non le signal des shorts channels (0 : non / 1 : oui)
    
    sub_stats = job.run(Hb_nostim_sc);
    
    dm = nirs.util.getdesign_matrix(Hb_nostim_sc(1));
    
    for i = 1:size(sub_stats,2)
        % extracts probe object to obtain channel numbers
        probe_tbl = data_cor(i).probe.link;
        % extracts index of short-channels in probe_tbl
        sc_idx = find(probe_tbl.ShortSeperation == 1);
        % extracts index of long-channels in probe_tbl
        lc_idx = find(probe_tbl.ShortSeperation == 0);
        % extracts PCA short-channel data (should be time*16)
        data_sc_PCA = data_cor(i).data(:,good_SC);
        % extracts all beta values from ols regression
        a = sub_stats(i).table;
        
        % creates a vector of beta corresponding to the right long-channel
        % indexed by lc_idx(j)
        beta_vec_sc = a.beta(a.source == probe_tbl.source(lc_idx(nxn)) & a.detector == probe_tbl.detector(lc_idx(nxn)) & strcmp(a.type,probe_tbl.type(lc_idx(nxn))));
        % creates a matrix of beta values of the size of data_sc to
        % multiply properly
        beta_matrix_sc = repmat(beta_vec_sc, [1, size(data_sc_PCA, 1)])';
        % sums all right side of the equation
        regressors_sum = data_sc_PCA.*beta_matrix_sc;
        % substracts, for each long-channel, the regressors_sum
        data_cor(i).data(:,lc_idx(nxn)) = data_cor(i).data(:,lc_idx(nxn)) - regressors_sum;
    end
    
    %disp(['long channel ' num2str(lc_idx(nxn)) ' done'])
end


dm = nirs.util.getdesign_matrix(Hb(1));
imagesc(dm);
% save


%% Données avec shorts channels


% Gets rid of all SC_PCA regressors
data_cor_nostim = discard_all_stims(data_cor);

% reinjects trial triggers for averaging
j = [];
j = nirs.modules.ChangeStimulusInfo();
j.ChangeTable = tbl_stims;
data_cor_stims = j.run(data_cor_nostim);

% export data to .nirs format
nirs.io.saveDotNirs(OD(1))

% Make average
j = [];
j = nirs.modules.Run_HOMER2();
j.fcn = 'hmrBlockAvg';
j.vars.trange = [-5 25];
average_screg = j.run(data_cor_stims);

% export data to R
r_table_average_close_sc = export_avg_data_to_R(average_screg);


table_test = r_table_average_close_sc(r_table_average_close_sc.type == "hbo" & r_table_average_close_sc.ID == "sub_1" & r_table_average_close_sc.channel == "1-1" & r_table_average_close_sc.stim == "trig1",:);
plot(table_test.data)
plots_average(average_screg,'hbo',[-20 25],"HbO average close SC correction")

%% Export fichier CSV

%Export les données avec short channel dans un tableau avec le en csv
writetable(r_table_average_close_sc,'H:/Autres ordinateurs/Mon ordinateur/Yann/Data Analyse/1. Analyse NIRS Toolbox/Analyse_TPP/Analyse CRNL/1. Script pour figure papier/V vs A R Analyses/Average_trial_block_close_sc/r_table_average_close_sc.csv','Delimiter',',')
