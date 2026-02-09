function data_sc = add_SC_regressors_mean_three_sc(data)
% Cette fonction ajoute des régresseurs basés sur la moyenne de canaux courts (SC) spécifiques.
% Au lieu d'une décomposition orthogonale, elle calcule la moyenne des signaux
% des SC pour HbO et HbR séparément et les ajoute comme régresseurs.

% Définir les indices des colonnes des canaux courts à utiliser.
% Dans une configuration standard, 3 de ces canaux seront HbO et 3 seront HbR.
sc_indices_to_use = [5, 6, 25, 26, 43, 44];

for i = 1:length(data)
    % Récupérer la table de description des canaux pour le sujet actuel
    link_table = data(i).probe.link;

    % --- Traitement pour les canaux HbO ---
    
    % Identifier les indices qui correspondent à des canaux HbO parmi la liste fournie
    hbo_indices = [];
    for k = 1:length(sc_indices_to_use)
        idx = sc_indices_to_use(k);
        % Vérifier que l'index est valide et que le type est 'hbo'
        if idx <= height(link_table) && strcmp(link_table.type{idx}, 'hbo')
            hbo_indices(end+1) = idx;
        end
    end
    
    % S'assurer qu'au moins un canal HbO a été trouvé
    if ~isempty(hbo_indices)
        % Extraire les données des canaux courts HbO sélectionnés
        dd_hbo = data(i).data(:, hbo_indices);
        
        % Calculer la moyenne de ces canaux (dimension 2 = à travers les colonnes)
        mean_hbo_sc = mean(dd_hbo, 2);
        
        % Créer et ajouter le vecteur moyen comme un nouveau régresseur
        st_hbo = nirs.design.StimulusVector;
        st_hbo.name = 'SC_mean_hbo'; % Nom du régresseur pour HbO
        st_hbo.time = data(i).time;
        st_hbo.vector = mean_hbo_sc;
        data(i).stimulus(st_hbo.name) = st_hbo;
    else
        warning('Sujet %d : Aucun canal court de type HbO trouvé avec les indices fournis.', i);
    end

    % --- Traitement pour les canaux HbR ---
    
    % Identifier les indices qui correspondent à des canaux HbR parmi la liste fournie
    hbr_indices = [];
    for k = 1:length(sc_indices_to_use)
        idx = sc_indices_to_use(k);
        % Vérifier que l'index est valide et que le type est 'hbr'
        if idx <= height(link_table) && strcmp(link_table.type{idx}, 'hbr')
            hbr_indices(end+1) = idx;
        end
    end

    % S'assurer qu'au moins un canal HbR a été trouvé
    if ~isempty(hbr_indices)
        % Extraire les données des canaux courts HbR sélectionnés
        dd_hbr = data(i).data(:, hbr_indices);
        
        % Calculer la moyenne de ces canaux
        mean_hbr_sc = mean(dd_hbr, 2);
        
        % Créer et ajouter le vecteur moyen comme un nouveau régresseur
        st_hbr = nirs.design.StimulusVector;
        st_hbr.name = 'SC_mean_hbr'; % Nom du régresseur pour HbR
        st_hbr.time = data(i).time;
        st_hbr.vector = mean_hbr_sc;
        data(i).stimulus(st_hbr.name) = st_hbr;
    else
        warning('Sujet %d : Aucun canal court de type HbR trouvé avec les indices fournis.', i);
    end
end

% Retourner la structure de données modifiée
data_sc = data;
end