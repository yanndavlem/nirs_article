classdef AddShortSeperationRegressors_mean_sc < nirs.modules.AbstractModule
    %% AddShortSeperationRegressors - Adds short seperation data as regressors to the GLM model
    %
    
    properties
         scICA;  % use single channel ICA instead of PCA for defining regressors
    end
    
    methods
        function obj = AddShortSeperationRegressors_mean_sc(prevJob )
            obj.name = 'AddShortSeperationRegressors_mean_sc';
            obj.scICA = false;
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
function data = runThis( obj, data )
    for i = 1:numel(data)
        % Vérifie si des canaux courts existent pour cette donnée
        if ~nirs.util.hasshortdistances(data(i))
            % S'il n'y en a pas, on passe à la donnée suivante
            continue;
        end
        
        % 1. Récupérer les indices de TOUS les canaux courts
        % lstss contiendra une liste d'indices, ex: [5, 6, 13, 14, 21, 22, ...]
        lstss = find(data(i).probe.link.ShortSeperation);
        
        % Si pour une raison quelconque, aucun canal court n'est trouvé, on continue
        if isempty(lstss)
            continue;
        end

        % 2. Séparer les indices en HbO (impairs) et HbR (pairs)
        % On sélectionne un indice sur deux, en commençant par le premier pour HbO
        idx_hbo = lstss(1:2:end);
        % On sélectionne un indice sur deux, en commençant par le deuxième pour HbR
        idx_hbr = lstss(2:2:end);
        
        % 3. Calculer la moyenne pour HbO et créer le régresseur
        if ~isempty(idx_hbo)
            % Extrait les données de tous les canaux courts HbO
            data_ss_hbo = data(i).data(:, idx_hbo);
            
            % Calcule la moyenne le long de la 2ème dimension (à travers les colonnes)
            % pour obtenir un unique vecteur de signal moyen
            mean_hbo = mean(data_ss_hbo, 2);
            
            % Crée la structure de régresseur pour la moyenne HbO
            st_hbo = nirs.design.StimulusVector;
            st_hbo.name = 'SS_Mean_HbO'; % Nom clair pour le régresseur
            st_hbo.time = data(i).time;
            st_hbo.vector = mean_hbo;
            st_hbo.regressor_no_interest = true; % Important pour le GLM
            
            % Ajoute le nouveau régresseur aux stimuli existants
            data(i).stimulus(st_hbo.name) = st_hbo;
        end
        
        % 4. Calculer la moyenne pour HbR et créer le régresseur
        if ~isempty(idx_hbr)
            % Extrait les données de tous les canaux courts HbR
            data_ss_hbr = data(i).data(:, idx_hbr);
            
            % Calcule la moyenne pour HbR
            mean_hbr = mean(data_ss_hbr, 2);

            % Crée la structure de régresseur pour la moyenne HbR
            st_hbr = nirs.design.StimulusVector;
            st_hbr.name = 'SS_Mean_HbR'; % Nom clair pour le régresseur
            st_hbr.time = data(i).time;
            st_hbr.vector = mean_hbr;
            st_hbr.regressor_no_interest = true;
            
            % Ajoute le nouveau régresseur
            data(i).stimulus(st_hbr.name) = st_hbr;
        end
    end
end
        
        
    end
    
end

