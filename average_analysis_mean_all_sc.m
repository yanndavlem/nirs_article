clc;clear;close all
cd('H:/Autres ordinateurs/Mon ordinateur/Yann/Data Analyse/1. Analyse NIRS Toolbox/Analyse_TPP/Analyse CRNL/Scripts_Matlab/')
% ajoute au path MATLAB la NIRSToolbox (mettre l'accès à la place des trois petits
% points
dir_nirtoolbox = 'C:/fNIRS AnalyzIR/nirs-toolbox-master/';
addpath(genpath(dir_nirtoolbox));
dir_hmr2 = 'H:/Autres ordinateurs/Mon ordinateur/Yann/Matériels utiles/Toolbox/homer2/';
addpath(genpath(dir_hmr2));

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
%% Preprocess
% Decoupe le signal pour ne garder que le signal intÃ©ressant
j=[];
j = nirs.modules.TrimBaseline( j );
j.preBaseline  = 10;
j.postBaseline = 30;
raw_trim = j.run(raw);

%change stimulus durations
raw_trim = nirs.design.change_stimulus_duration(raw_trim,'trig1',10);
raw_trim = nirs.design.change_stimulus_duration(raw_trim,'trig2',10);
raw_trim = nirs.design.change_stimulus_duration(raw_trim,'trig3',10);

bad_chan = bad_channels_indices(raw_trim,0.7);

%pretraitement des donnÃ©es par la NIRS Toolbox

% Stocke pour plus tard la table de triggers Ã  rÃ©injecter Ã  la fin
tbl_stims = nirs.createStimulusTable(raw_trim); %Recupere les triggers


%pretraitement des donnÃ©es par la NIRS Toolbox

jobs = [];
jobs = nirs.modules.OpticalDensity(jobs);
OD = jobs.run(raw_trim); % transforme tes donnÃ©es en densitÃ© optiques

jobs = [];
jobs = nirs.modules.TDDR(); % ajoute aux "jobs" la correction TDDR
jobs = eeg.modules.BandPassFilter(jobs); % ajoute aux jobs un filtre bandpass
jobs.do_downsample = 0;
jobs.lowpass = 0.12; % choisis ton passe-bas
jobs.highpass = 0.01; % choisis ton passe-haut
jobs = nirs.modules.BeerLambertLaw(jobs); % ajoute aux jobs la transformation de BeerLambert
Hb = jobs.run(OD); % applique tous les jobs prÃ©cÃ©dent aux densitÃ©s optiques, tu obtiens des donnÃ©es d'hÃ©moglobine oxy/desoxy

% Resample
j = [];
j = nirs.modules.Resample( j );
j.Fs = 2; 
raw = j.run(raw);



%% DonnÃ©es sans shorts channels 

% % Make average
% j = [];
% j = nirs.modules.Run_HOMER2();
% j.fcn = 'hmrBlockAvg';
% j.vars.trange = [-5 25];
% average = j.run(Hb);
% 
% plots_average(average,'hbo',[-5 25],"HbO average without SC correction")
% 
% % exports to a table ready for R
% r_table_average_no_sc = export_avg_data_to_R(average);
% 
% %r_table(r_table.type == "hbo" & r_table.ID == "sub_1" & r_table.channel == "1-1" & r_table.stim == "trig1",:);


%% Données avec shorts channels 

% 1. Supprime les triggers d'essais pour ne pas qu'ils influencent le calcul des bêtas
Hb_nostim = discard_all_stims(Hb);

% 2. Ajoute les deux régresseurs moyens (SC_mean_hbo et SC_mean_hbr) aux données
Hb_nostim_sc = add_SC_regressors_mean_all_sc(Hb_nostim);

% 3. Exécute le GLM pour obtenir les poids bêta de chaque régresseur sur chaque canal
job = [];
job = nirs.modules.GLM;
job.type = 'OLS';
% Cette option est importante, elle dit au GLM d'utiliser les régresseurs que NOUS avons créés
job.AddShortSepRegressors =  0; 
sub_stats = job.run(Hb_nostim_sc);

% (Optionnel) Visualisation de la matrice de design pour un sujet
% Vous devriez voir des colonnes pour SC_mean_hbo et SC_mean_hbr
dm = nirs.util.getdesign_matrix(Hb_nostim_sc(1));
figure;
imagesc(dm);
title('Matrice de Design (avec régresseurs SC)');
colorbar;


% 4. NOUVELLE BOUCLE DE CORRECTION : Calcule et soustrait le bruit estimé
% =========================================================================
disp('Application de la correction basée sur les régresseurs SC moyens...');
data_cor = Hb_nostim_sc; % On part des données qui contiennent les régresseurs

for i = 1:length(data_cor)
    % Récupère la table des bêtas calculés pour ce sujet
    beta_table = sub_stats(i).table;
    % Récupère la description des canaux (source, détecteur, type)
    probe_tbl = data_cor(i).probe.link;
    % Récupère les indices de tous les canaux longs
    lc_indices = find(probe_tbl.ShortSeperation == 0);

    % Extrait les vecteurs de nos deux régresseurs pour ce sujet
    regressor_hbo = data_cor(i).stimulus('SC_mean_hbo').vector;
    regressor_hbr = data_cor(i).stimulus('SC_mean_hbr').vector;
    
    % Boucle sur chaque canal long pour lui appliquer la correction
    for j = 1:length(lc_indices)
        chan_idx = lc_indices(j);

        % Informations sur le canal long actuel
        chan_source = probe_tbl.source(chan_idx);
        chan_detector = probe_tbl.detector(chan_idx);
        chan_type = probe_tbl.type{chan_idx}; % 'hbo' ou 'hbr'

        correction_term = 0; % Initialise le terme de correction

        % La correction dépend du type de canal (HbO ou HbR)
        if strcmp(chan_type, 'hbo')
            % Trouve la ligne de bêta qui correspond à CE canal et au régresseur HBO
            beta_row = beta_table.source == chan_source & ...
                       beta_table.detector == chan_detector & ...
                       strcmp(beta_table.type, 'hbo') & ...
                       strcmp(beta_table.cond, 'SC_mean_hbo');
            
            if any(beta_row)
                beta_val = beta_table.beta(beta_row);
                % Le bruit à soustraire = bêta * signal du régresseur
                correction_term = beta_val * regressor_hbo;
            end

        elseif strcmp(chan_type, 'hbr')
            % Trouve la ligne de bêta qui correspond à CE canal et au régresseur HBR
            beta_row = beta_table.source == chan_source & ...
                       beta_table.detector == chan_detector & ...
                       strcmp(beta_table.type, 'hbr') & ...
                       strcmp(beta_table.cond, 'SC_mean_hbr');
            
            if any(beta_row)
                beta_val = beta_table.beta(beta_row);
                % Le bruit à soustraire = bêta * signal du régresseur
                correction_term = beta_val * regressor_hbr;
            end
        end

        % Applique la correction en soustrayant le bruit estimé
        data_cor(i).data(:, chan_idx) = data_cor(i).data(:, chan_idx) - correction_term;
    end   
end
disp('Correction terminée.');
% =========================================================================


% 5. Supprime TOUS les "stimulus" (y compris nos régresseurs SC) des données corrigées
data_cor_nostim = discard_all_stims(data_cor);

% 6. Réinjecte les VRAIS triggers d'essais pour pouvoir faire la moyenne par blocs
j = [];
j = nirs.modules.ChangeStimulusInfo();
j.ChangeTable = tbl_stims;
data_cor_stims = j.run(data_cor_nostim);

% 7. Fait la moyenne par blocs sur les données maintenant corrigées
j = [];
j = nirs.modules.Run_HOMER2();
j.fcn = 'hmrBlockAvg';
j.vars.trange = [-5 25];
average_screg = j.run(data_cor_stims);

% 8. Exporte les données pour R et visualise les résultats
r_table_average_with_sc = export_avg_data_to_R(average_screg);
plots_average(average_screg,'hbo',[-5 25],"HbO average AVEC correction SC (moyenne)")

%% Export fichier CSV

% Export des données AVEC correction par la moyenne des SC
writetable(r_table_average_with_sc, 'H:/Autres ordinateurs/Mon ordinateur/Yann/Data Analyse/1. Analyse NIRS Toolbox/Analyse_TPP/Analyse CRNL/1. Script pour figure papier/V vs A R Analyses/Average_trial_block_with_mean_all_sc/r_table_average_with_mean_all_sc.csv', 'Delimiter', ',');