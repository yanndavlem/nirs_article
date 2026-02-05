%% Importer les donn�es + mettre les chemins pour les toolbox

clc;clear;close all
cd('H:/Autres ordinateurs/Mon ordinateur/Yann/Data Analyse/1. Analyse NIRS Toolbox/Analyse_TPP/Analyse CRNL/Scripts_Matlab/')
% ajoute au path MATLAB la NIRSToolbox (mettre l'accès à la place des trois petits
% points
dir_nirtoolbox = 'C:/fNIRS AnalyzIR/nirs-toolbox-master/';
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

%% Preprocess

% % Resample
% j = [];
% j = nirs.modules.Resample( j );
% j.Fs = 2; 
% raw = j.run(raw);

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

bad_chan = bad_channels_indices(raw_trim,0.6);

%pretraitement des données par la NIRS Toolbox
jobs = [];
jobs = nirs.modules.OpticalDensity(jobs);
OD = jobs.run(raw_trim); % transforme tes données en densité optiques

jobs = [];
jobs = nirs.modules.TDDR(); % ajoute aux "jobs" la correction TDDR
jobs = eeg.modules.BandPassFilter(jobs); % ajoute aux jobs un filtre bandpass
jobs.lowpass = 0.12; % choisis ton passe-bas
jobs.highpass = 0.01; % choisis ton passe-haut
jobs = nirs.modules.BeerLambertLaw(jobs); % ajoute aux jobs la transformation de BeerLambert
Hb = jobs.run(OD); % applique tous les jobs précédent aux densités optiques, tu obtiens des données d'hémoglobine oxy/desoxy



% GLM
job = [];
job = nirs.modules.GLM;
job.type = 'OLS'; % Choisi le type de GLM
job.AddShortSepRegressors = 1; % D�cide de prendre ou non le signal des shorts channels (0 : non / 1 : oui)


dm = nirs.util.getdesign_matrix(Hb(1));
imagesc(dm);
% save

sub_stats = job.run(Hb);
save('H:/Autres ordinateurs/Mon ordinateur/Yann/Data Analyse/1. Analyse NIRS Toolbox/Analyse_TPP/Analyse CRNL/Script_R/Resultat_GLM/OLS_trial_block/stats_trial_block_OLS.mat','sub_stats')

%cr�� un fichier par sujet qu'il faudra nettoyer dans R
for i = 1:size(sub_stats,2)
a = sub_stats(i).table;
writetable(a,['H:/Autres ordinateurs/Mon ordinateur/Yann/Data Analyse/1. Analyse NIRS Toolbox/Analyse_TPP/Analyse CRNL/1. Script pour figure papier/V vs A R Analyses/OLS_trial_block/suj_' num2str(i) '_all_sc_v2.csv'],'Delimiter',',')
end

output_path = 'H:/Autres ordinateurs/Mon ordinateur/Yann/Data Analyse/1. Analyse NIRS Toolbox/Analyse_TPP/Analyse CRNL/1. Script pour figure papier/SCI_All_analysis/bad_channels_list_all_sc.csv';
writetable(bad_chan, output_path)

