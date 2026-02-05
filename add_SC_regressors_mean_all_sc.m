function data_sc = add_SC_regressors_mean_all_sc(data)
% Cette fonction ajoute des régresseurs basés sur la moyenne des signaux de tous les canaux courts (SC).
% Elle identifie automatiquement les canaux SC pour chaque sujet, les sépare en HbO et HbR,
% calcule la moyenne de chaque groupe, et les ajoute comme deux régresseurs distincts.

for i = 1:length(data)
    % 1. Récupérer la table de description des canaux pour le sujet actuel
    link_table = data(i).probe.link;
    
    % 2. Trouver les indices de TOUS les canaux courts pour CE sujet
    all_sc_indices = find(link_table.ShortSeperation == 1);
    
    % 3. Séparer ces indices en listes pour HbO et HbR en vérifiant le type de chacun
    hbo_sc_indices = [];
    hbr_sc_indices = [];
    
    for j = 1:length(all_sc_indices)
        current_idx = all_sc_indices(j);
        if strcmp(link_table.type{current_idx}, 'hbo')
            hbo_sc_indices(end+1) = current_idx;
        elseif strcmp(link_table.type{current_idx}, 'hbr')
            hbr_sc_indices(end+1) = current_idx;
        end
    end
    
    % --- 4. Traitement pour les canaux HbO ---
    if ~isempty(hbo_sc_indices)
        % Extraire les données des canaux courts HbO identifiés
        data_hbo_sc = data(i).data(:, hbo_sc_indices);
        % Calculer la moyenne de ces signaux
        mean_hbo_sc = mean(data_hbo_sc, 2);
        
        % Créer et ajouter le régresseur pour HbO
        st_hbo = nirs.design.StimulusVector;
        st_hbo.name = 'SC_mean_hbo';
        st_hbo.time = data(i).time;
        st_hbo.vector = mean_hbo_sc;
        data(i).stimulus(st_hbo.name) = st_hbo;
    else
        warning('Sujet %d : Aucun canal court de type HbO trouvé.', i);
    end
    
    % --- 5. Traitement pour les canaux HbR ---
    if ~isempty(hbr_sc_indices)
        % Extraire les données des canaux courts HbR identifiés
        data_hbr_sc = data(i).data(:, hbr_sc_indices);
        % Calculer la moyenne de ces signaux
        mean_hbr_sc = mean(data_hbr_sc, 2);
        
        % Créer et ajouter le régresseur pour HbR
        st_hbr = nirs.design.StimulusVector;
        st_hbr.name = 'SC_mean_hbr';
        st_hbr.time = data(i).time;
        st_hbr.vector = mean_hbr_sc;
        data(i).stimulus(st_hbr.name) = st_hbr;
    else
        warning('Sujet %d : Aucun canal court de type HbR trouvé.', i);
    end
end

% Retourner la structure de données modifiée
data_sc = data;
end