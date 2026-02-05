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

%% Recupere le signal des 15 premiers stims
% (Section commentée dans votre original, laissée telle quelle)
% tbl = nirs.createStimulusTable(raw); %Recupere les triggers
%
% for m = 1:size(raw,1)
% ...
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

% % Resample
% (Section commentée dans votre original, laissée telle quelle)
% j = [];
% j = nirs.modules.Resample( j );
% j.Fs = 2;
% raw = j.run(raw);


%% Données avec shorts channels (LOGIQUE MODIFIÉE)

lc_idx = find(Hb(1).probe.link.ShortSeperation == 0); % recupere les numeros de canaux des long channels
sc_idx = find(Hb(1).probe.link.ShortSeperation == 1); % recupere les numeros de canaux des shorts channels
conditions = {'trig1','trig2','trig3'};

% Note : Les variables 'sz', 'varNames', 'table_close_SC', 'flag' ont été
% supprimées car elles n'étaient pas utilisées dans ce script d'averaging.

% Gets rid of all stims
Hb_nostim = discard_all_stims(Hb);
data_cor = Hb_nostim;


for nxn = 1:length(lc_idx)
    disp(['long channel ' num2str(lc_idx(nxn)) ' start'])
    
    % ====================================================================
    % DÉBUT DE LA LOGIQUE : TROUVER LE MEILLEUR SC (PARMI TOUS)
    % ====================================================================
    
    current_lc_global_idx = lc_idx(nxn); % Index global du LC courant
    
    % On va stocker les corrélations pour ce LC (nxn) à travers tous les sujets
    % Lignes: Sujets, Colonnes: nombre total de SCs
    all_corrs_this_lc = zeros(length(Hb_nostim), length(sc_idx));
    
    for i = 1:length(Hb_nostim) % Boucle sur chaque sujet
        % Récupère la série temporelle pour le LC courant
        current_lc_data = Hb_nostim(i).data(:, current_lc_global_idx);
        
        % Récupère les séries temporelles pour TOUS les SCs de ce sujet
        all_sc_data = Hb_nostim(i).data(:, sc_idx);
        
        % Calcule la corrélation
        corrs_vector = corr(current_lc_data, all_sc_data);
        corrs_vector(isnan(corrs_vector)) = 0; % Gère les NaNs
        
        all_corrs_this_lc(i, :) = corrs_vector;
    end
    
    % Calcule la moyenne des corrélations *absolues*
    avg_abs_corrs = mean(abs(all_corrs_this_lc), 1);
    
    % Trouve l'index *local* (parmi les SCs) du SC avec la plus haute corrélation
    [~, best_sc_local_idx] = max(avg_abs_corrs);
    
    % Récupère l'index *global* de ce SC (l'indice réel dans Hb.data)
    best_sc_global_idx = sc_idx(best_sc_local_idx);
    
    % ====================================================================
    % FIN DE LA LOGIQUE DE CORRÉLATION
    % ====================================================================
    
    % ====================================================================
    % DÉBUT : TROUVE LA PAIRE DU MEILLEUR SC
    % ====================================================================
    
    % Trouve l'indice de la paire (HbO/HbR) du meilleur SC
    if mod(best_sc_global_idx, 2) == 1 % Si l'indice est impair (ex: 5, 11, ...)
        pair_global_idx = best_sc_global_idx + 1; % Sa paire est l'indice pair suivant
    else % Si l'indice est pair (ex: 6, 12, ...)
        pair_global_idx = best_sc_global_idx - 1; % Sa paire est l'indice impair précédent
    end
    
    % Définit la variable 'good_SC' comme étant le VECTEUR de la paire
    good_SC = sort([best_sc_global_idx, pair_global_idx]);
    
    % fprintf('  -> LC %d: Best correlated global SC index %d. Using pair [%d, %d].\n', ...
    %          current_lc_global_idx, best_sc_global_idx, good_SC(1), good_SC(2));
    
    % ====================================================================
    % FIN DE LA LOGIQUE DE PAIRING
    % (L'ancien bloc "if ismember..." est remplacé par ce qui précède)
    % ====================================================================

    
    % Adds the SC signal (basé sur la fonction add_SC_regressors_one_sc.m)
    % Utilise maintenant la paire 'good_SC' dynamique
    Hb_nostim_sc = add_SC_regressors_one_sc(Hb_nostim,good_SC);
    
    % GLM
    job = [];
    job = nirs.modules.GLM;
    job.type = 'OLS'; % Choisi le type de GLM
    job.AddShortSepRegressors= 0; % On le fait manuellement
    
    sub_stats = job.run(Hb_nostim_sc);
    
    dm = nirs.util.getdesign_matrix(Hb_nostim_sc(1));
    
    for i = 1:size(sub_stats,2)
        % extracts probe object to obtain channel numbers
        probe_tbl = data_cor(i).probe.link;
        % extracts index of short-channels in probe_tbl
        % sc_idx = find(probe_tbl.ShortSeperation == 1); % (Déjà défini en dehors de la boucle)
        % extracts index of long-channels in probe_tbl
        % lc_idx = find(probe_tbl.ShortSeperation == 0); % (Déjà défini en dehors de la boucle)

        % extracts PCA short-channel data (should be time*2)
        data_sc_PCA = data_cor(i).data(:,good_SC);
        % extracts all beta values from ols regression
        a = sub_stats(i).table;
        
        % ====================================================================
        % DÉBUT DU BLOC D'EXTRACTION DE BÊTA CORRIGÉ
        % ====================================================================

        % Sépare la table 'a' pour ne garder que le LC courant
        lc_rows_idx = (a.source == probe_tbl.source(lc_idx(nxn)) & ...
                       a.detector == probe_tbl.detector(lc_idx(nxn)) & ...
                       strcmp(a.type, probe_tbl.type(lc_idx(nxn))));
        a_lc_table = a(lc_rows_idx, :);

        % Extrait UNIQUEMENT les bêtas des régresseurs SC
        % (Basé sur la fct add_SC_regressors_one_sc.m qui nomme 'SC_1_PCA', 'SC_2_PCA')
        beta_vec_sc = a_lc_table.beta(startsWith(a_lc_table.cond, 'SC_'));
        
        % Vérification : s'assure qu'on a le bon nombre de bêtas
        if length(beta_vec_sc) ~= length(good_SC)
           warning('Problème de bêta pour LC %d, Sujet %d. Bêtas trouvés: %d, attendus: %d', ...
                   lc_idx(nxn), i, length(beta_vec_sc), length(good_SC));
           % Remplace par des zéros pour éviter un crash
           beta_vec_sc = zeros(length(good_SC), 1);
        end
        % ====================================================================
        % FIN DU BLOC D'EXTRACTION DE BÊTA CORRIGÉ
        % ====================================================================

        % creates a matrix of beta values of the size of data_sc to
        % multiply properly
        beta_matrix_sc = repmat(beta_vec_sc, [1, size(data_sc_PCA, 1)])';
        % sums all right side of the equation
        regressors_sum = sum([data_sc_PCA.*beta_matrix_sc],2);
        % substracts, for each long-channel, the regressors_sum
        data_cor(i).data(:,lc_idx(nxn)) = data_cor(i).data(:,lc_idx(nxn)) - regressors_sum;
    end
    
    %disp(['long channel ' num2str(lc_idx(nxn)) ' done'])
end

% saves regressed data for all participants
%save('../results/average_scregPCA_delay/data/data_after_reg.mat', 'data_cor')

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
r_table_average_with_sc = export_avg_data_to_R(average_screg);



%% Export fichier CSV

%Export les données avec short channel dans un tableau avec le en csv
writetable(r_table_average_with_sc,'H:/Autres ordinateurs/Mon ordinateur/Yann/Data Analyse/1. Analyse NIRS Toolbox/Analyse_TPP/Analyse CRNL/1. Script pour figure papier/V vs A R Analyses/Average_trial_corr_all_sc/r_table_average_most_correlated_pair.csv','Delimiter',',')
% ATTENTION: J'ai changé le nom du fichier de sortie pour ne pas écraser l'ancien
% 'r_table_average_one_sc.csv' est devenu 'r_table_average_most_correlated_pair.csv'