% SUR DE CETTE ANALYSE, MOINS DE L'ANALYSE ALL CORRELATED

clc;clear;close all
cd('H:/Autres ordinateurs/Mon ordinateur/Yann/Data Analyse/1. Analyse NIRS Toolbox/Analyse_TPP/Analyse CRNL/Scripts_Matlab/')
% ajoute au path MATLAB la NIRSToolbox
dir_nirtoolbox = 'C:/fNIRS AnalyzIR/nirs-toolbox-master/';
addpath(genpath(dir_nirtoolbox));

tbl_roi = readtable('H:/Autres ordinateurs/Mon ordinateur/Yann/Data Analyse/1. Analyse NIRS Toolbox/Analyse_TPP/Analyse CRNL/Scripts_Matlab/correspondance_one_channel.csv');

% importation des donnees brutes
folders = dir('H:/Autres ordinateurs/Mon ordinateur/Yann/Data Analyse/1. Analyse NIRS Toolbox/Analyse_TPP/Analyse CRNL/fichiers_nirx/*2024*'); 
for i = 1:length(folders)
    raw(i,1) =  nirs.io.loadNIRx(['H:/Autres ordinateurs/Mon ordinateur/Yann/Data Analyse/1. Analyse NIRS Toolbox/Analyse_TPP/Analyse CRNL/fichiers_nirx/' folders(i).name],true); 
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

% Récupère les triggers manquants
tbl = nirs.createStimulusTable(raw); %Recupere les triggers

% Vérifie les liens où les triggers sont absents
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

% Ajoute les triggers manquant grâce au fichier datalog txt qui a été préalablement
% converti en fichier csv
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
    
    % Ajoute une ligne spécial pour le 7e sujet, il manque le dernier trigger
    if j == 7
        new_trig_special_7 = time_trig_all(end,1) + 10 + 10;
        raw(j).stimulus('trig1').onset = [raw(j).stimulus('trig1').onset; new_trig_special_7];
        raw(j).stimulus('trig1').dur = repmat(1,length(raw(j).stimulus('trig1').onset),1);
        raw(j).stimulus('trig1').amp = repmat(1,length(raw(j).stimulus('trig1').onset),1);
    end
end

%% Preprocess

% Resample
j = [];
j = nirs.modules.Resample( j );
j.Fs = 2;
raw = j.run(raw);

% Decoupe le signal pour ne garder que le signal intéressant
j=[];
j = nirs.modules.TrimBaseline( j );
j.preBaseline  = 10;
j.postBaseline = 20;
raw_trim = j.run(raw);

%change stimulus durations
raw_trim = nirs.design.change_stimulus_duration(raw_trim,'trig1',10);
raw_trim = nirs.design.change_stimulus_duration(raw_trim,'trig2',10);
raw_trim = nirs.design.change_stimulus_duration(raw_trim,'trig3',10);

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

one_HB = [];
full_sub_stats = [];

%% Fait l'analyse pour ne prendre que la régression du SC LE PLUS CORRÉLÉ

lc_idx = find(Hb(1).probe.link.ShortSeperation == 0); % recupere les numeros de canaux des long channels
sc_idx = find(Hb(1).probe.link.ShortSeperation == 1); % recupere les numeros de canaux des shorts channels
conditions = {'trig1','trig2','trig3'};
sz = [size(Hb,1)*length(lc_idx)*3,6];
varNames = {'source','detector','type','cond','beta','sub'};
varTypes = {'double','double','string','string','double','double'};
table_close_SC = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
tbl_probe = Hb(1).probe.link;

flag = 0;
% full_SC_associate_LC = []; % On n'utilise plus ça

for nxn = 1:length(lc_idx) % Boucle sur chaque Long Channel
    
    current_lc_global_idx = lc_idx(nxn);
    fprintf('Processing LC %d/%d (Global Index: %d)\n', nxn, length(lc_idx), current_lc_global_idx);

% ====================================================================
% DÉBUT DE LA LOGIQUE MODIFIÉE : TROUVER ET PASSER LA PAIRE D'INDICES
% ====================================================================

% Définit les paires de SC cibles dans une "cell array"
target_sc_pairs = {[5, 6], [25, 26], [43, 44]};

% Aplatit la liste pour le calcul de corrélation
target_sc_indices = [target_sc_pairs{:}]; % Devient [5, 6, 25, 26, 43, 44]
num_target_scs = length(target_sc_indices);

% On va stocker les corrélations pour ce LC (nxn) à travers tous les sujets
% Lignes: Sujets, Colonnes: 6 SCs cibles
all_corrs_this_lc = zeros(length(Hb), num_target_scs);

for i = 1:length(Hb) % Boucle sur chaque sujet
    % Récupère la série temporelle pour le LC courant et le sujet courant
    current_lc_data = Hb(i).data(:, current_lc_global_idx);
    
    % Récupère les séries temporelles pour les 6 SCs cibles
    if max(target_sc_indices) > size(Hb(i).data, 2)
         error('Indice SC cible (%d) hors limites pour Hb.data (taille: %d)', max(target_sc_indices), size(Hb(i).data, 2));
    end
    target_sc_data = Hb(i).data(:, target_sc_indices);
    
    % Calcule la corrélation entre ce LC et les 6 SCs cibles
    corrs_vector = corr(current_lc_data, target_sc_data);
    corrs_vector(isnan(corrs_vector)) = 0; % Gère les NaNs
    
    % Stocke les corrélations
    all_corrs_this_lc(i, :) = corrs_vector;
end

% Calcule la moyenne des corrélations *absolues* à travers tous les sujets
% On obtient un vecteur (1 x 6)
avg_abs_corrs = mean(abs(all_corrs_this_lc), 1);

% Trouve l'index *local* (parmi les 6 SCs) du SC avec la plus haute corrélation
% best_sc_local_idx sera un nombre de 1 à 6
[~, best_sc_local_idx] = max(avg_abs_corrs);

% Mappe cet index (1-6) à un index de PAIRE (1, 2, ou 3)
% Si best_sc_local_idx est 1 ou 2, pair_index = 1
% Si best_sc_local_idx est 3 ou 4, pair_index = 2
% Si best_sc_local_idx est 5 ou 6, pair_index = 3
pair_index = ceil(best_sc_local_idx / 2);

% Récupère la PAIRE d'indices correspondante (ex: [5, 6])
SC_associate_LC = target_sc_pairs{pair_index};

% [Optionnel] Logger
% fprintf('  -> Best correlated SC local idx %d, maps to pair index %d. Passing indices [%d, %d] to GLM.\n', ...
%          best_sc_local_idx, pair_index, SC_associate_LC(1), SC_associate_LC(2));

% ====================================================================
% FIN DE LA LOGIQUE MODIFIÉE
% ====================================================================
    % GLM
    job = [];
    job = nirs.modules.GLM_corr_three_SC; % Ton module custom
    job.type = 'OLS'; % Choisi le type de GLM
    job.AddShortSepRegressors_corr_three_SC = 1; % Décide de prendre ou non le signal des shorts channels (0 : non / 1 : oui)
    job.whichSC = SC_associate_LC; % <-- On passe l'index relatif du SC le plus corrélé
    
    
    sub_stats = job.run(Hb);
    
    for jx = 1:size(sub_stats,2)
        for cond = 1:length(conditions)
            ab = sub_stats(jx).table;
            
            % Trouve la ligne de résultat pour le LC courant (nxn)
            ab_idx = find(ab.source == tbl_probe.source(current_lc_global_idx) & ...
                          ab.detector == tbl_probe.detector(current_lc_global_idx) & ...
                          strcmp(ab.type, tbl_probe.type(current_lc_global_idx)) & ...
                          strcmp(ab.cond, conditions{cond}) );
            
            % Vérifie qu'on a bien trouvé une ligne (évite les erreurs)
            if ~isempty(ab_idx)
                % Il peut y avoir plusieurs types (hbo/hbr)
                for find_idx = 1:length(ab_idx)
                    flag = flag +1;
                    idx = ab_idx(find_idx); % Index de la ligne trouvée
                    table_close_SC.source(flag) = ab.source(idx);
                    table_close_SC.detector(flag) = ab.detector(idx);
                    table_close_SC.type(flag) = ab.type(idx);
                    table_close_SC.cond(flag) = ab.cond(idx);
                    table_close_SC.beta(flag) = ab.beta(idx);
                    table_close_SC.sub(flag) = jx;
                end
            else
                % Affiche un avertissement si aucune stat n'est trouvée
                % fprintf('Attention: Stats non trouvées pour LC idx %d, Sujet %d, Cond %s\n', ...
                %         current_lc_global_idx, jx, conditions{cond});
            end
        end
    end
end


% Nettoie la table des lignes vides pré-allouées
table_close_SC = table_close_SC(1:flag,:);


dm = nirs.util.getdesign_matrix(Hb(1));
imagesc(dm);
% save


save('H:/Autres ordinateurs/Mon ordinateur/Yann/Data Analyse/1. Analyse NIRS Toolbox/Analyse_TPP/Analyse CRNL/Script_R/Resultat_GLM/OLS_trial_block/stats_trial_block_OLS_MostCorrelated.mat','sub_stats')

%créé un fichier par sujet qu'il faudra nettoyer dans R
    writetable(table_close_SC,'H:/Autres ordinateurs/Mon ordinateur/Yann/Data Analyse/1. Analyse NIRS Toolbox/Analyse_TPP/Analyse CRNL/1. Script pour figure papier/V vs A R Analyses/OLS_trial_block_correlated_three/beta_table_GLM_MostCorrelated_three_SC_NEW.csv','Delimiter',',')