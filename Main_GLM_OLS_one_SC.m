clc;clear;close all
cd('H:/Autres ordinateurs/Mon ordinateur/Yann/Data Analyse/1. Analyse NIRS Toolbox/Analyse_TPP/Analyse CRNL/Scripts_Matlab/')
% ajoute au path MATLAB la NIRSToolbox (mettre l'accès à la place des trois petits
% points
dir_nirtoolbox = 'C:/fNIRS AnalyzIR/nirs-toolbox-master/';
addpath(genpath(dir_nirtoolbox));

tbl_roi = readtable('H:/Autres ordinateurs/Mon ordinateur/Yann/Data Analyse/1. Analyse NIRS Toolbox/Analyse_TPP/Analyse CRNL/Scripts_Matlab/correspondance_one_channel.csv');

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

%% Fait l'analyse pour ne prendre que la r�gression d'un short channel

lc_idx = find(Hb(1).probe.link.ShortSeperation == 0); % recupere les numeros de canaux des long channels
sc_idx = find(Hb(1).probe.link.ShortSeperation == 1); % recupere les numeros de canaux des shorts channels
lc_roi = tbl_roi(tbl_roi.SC == 0,:); %r�cup�re les LC sur les ROI
conditions = {'trig1','trig2','trig3'};
sz = [size(Hb,1)*length(lc_idx)*3,6];
varNames = {'source','detector','type','cond','beta','sub'};
varTypes = {'double','double','string','string','double','double'};
table_close_SC = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
flag = 0;
%Simplifier la lecture
tbl_probe = Hb(1).probe.link;
full_SC_associate_LC = [];

% Analayse one channel per lobe
% roi = ['temp_g','occi','temp_d'];

for nxn = 1:height(lc_roi)
    
    %Recupere les SC par rapport aux longs channels correspondant
    if strcmp(lc_roi.ROI(nxn),'temp_g')
        if strcmp(lc_roi.Type(nxn),'''hbo''')
            good_SC = 5;
        else
            good_SC = 6;
        end
    elseif strcmp(lc_roi.ROI(nxn),'occi')
        if strcmp(lc_roi.Type(nxn),'''hbo''')
            good_SC = 25;
        else
            good_SC = 26;
        end
        elseif strcmp(lc_roi.ROI(nxn),'temp_d')
        if strcmp(lc_roi.Type(nxn),'''hbo''')
            good_SC = 43;
        else
            good_SC = 44;
        end
    end

    %Attribue le SC associe au long channel
    SC_associate_LC = find(sc_idx == good_SC);
   %check les SC s�lectionn�s
   full_SC_associate_LC = [full_SC_associate_LC;nxn,lc_roi.ROI(nxn),good_SC];

    % GLM
    job = [];
    job = nirs.modules.GLM_close_SC;
    job.type = 'OLS'; % Choisi le type de GLM
    job.AddShortSepRegressors_close_SC = 1; % D�cide de prendre ou non le signal des shorts channels (0 : non / 1 : oui)
    job.whichSC = SC_associate_LC;
    
    
    sub_stats = job.run(Hb);
    for jx = 1:size(sub_stats,2)
        for cond = 1:length(conditions)
            ab = sub_stats(jx).table;
            
            ab_idx = find(ab.source == tbl_probe.source(lc_idx(nxn)) & ab.detector == tbl_probe.detector(lc_idx(nxn)) & strcmp(ab.type,tbl_probe.type(lc_idx(nxn))) & strcmp(ab.cond,conditions {cond}) );
            flag = flag +1;
            table_close_SC.source(flag) = ab.source(ab_idx);
            table_close_SC.detector(flag) = ab.detector(ab_idx);
            table_close_SC.type(flag) = ab.type(ab_idx);
            table_close_SC.cond(flag) = ab.cond(ab_idx);
            table_close_SC.beta(flag) = ab.beta(ab_idx);
            table_close_SC.sub(flag) = jx;
            
        end
    end
end

dm = nirs.util.getdesign_matrix(Hb(1));

imagesc(dm);
% save


save('H:/Autres ordinateurs/Mon ordinateur/Yann/Data Analyse/1. Analyse NIRS Toolbox/Analyse_TPP/Analyse CRNL/Script_R/Resultat_GLM/OLS_trial_block/stats_trial_block_OLS.mat','sub_stats')

%cr�� un fichier par sujet qu'il faudra nettoyer dans R
writetable(table_close_SC,'H:/Autres ordinateurs/Mon ordinateur/Yann/Data Analyse/1. Analyse NIRS Toolbox/Analyse_TPP/Analyse CRNL/1. Script pour figure papier/V vs A R Analyses/OLS_trial_block_one_sc/beta_table_GLM_one_SC_NEW.csv','Delimiter',',')

