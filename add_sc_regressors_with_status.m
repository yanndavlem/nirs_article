function [data_out, sc_regressors_valid] = add_sc_regressors_with_status(data_in, subject_id_for_log_external)
    data_out = data_in;
    sc_regressors_valid = false;

    % Utiliser l'ID externe pour les logs si fourni, sinon prendre de data_in.description
    if nargin >= 2 && ~isempty(subject_id_for_log_external)
        current_id_for_log = subject_id_for_log_external;
    elseif ~isempty(data_in.description)
        current_id_for_log = data_in.description;
    else
        current_id_for_log = 'ID_Manquant_Dans_Fonction';
    end

    if ~nirs.util.hasshortdistances(data_out) || isempty(data_out.data)
        fprintf('Sujet %s: Pas de canaux courts définis ou pas de données.\n', current_id_for_log);
        return;
    end

    lstss = find(data_out.probe.link.ShortSeperation);
    if isempty(lstss)
        fprintf('Sujet %s: Aucun canal marqué comme ShortSeperation=true.\n', current_id_for_log);
        return;
    end
    
    dd_orig = data_out.data(:, lstss);

    if size(dd_orig,1) < 10 || size(dd_orig,2) == 0
        fprintf('Sujet %s: Pas assez de données temporelles ou pas de canaux courts pour les régresseurs SC.\n', current_id_for_log);
        return;
    end
    
    all_nan_cols = all(isnan(dd_orig), 1);
    if all(all_nan_cols)
        fprintf('Sujet %s: Toutes les données des canaux courts sont NaN. Impossible de générer des régresseurs SC.\n', current_id_for_log);
        return;
    end
    dd_orig(:, all_nan_cols) = 0;

    dd_fixed = dd_orig;
    for k_col = 1:size(dd_fixed, 2)
        if any(isnan(dd_fixed(:, k_col)))
            dd_fixed(:, k_col) = fillmissing(dd_fixed(:, k_col), 'linear', 'EndValues', 'nearest');
            if any(isnan(dd_fixed(:, k_col)))
                dd_fixed(isnan(dd_fixed(:, k_col)), k_col) = 0;
            end
        end
    end
    
    dd_centered = dd_fixed - mean(dd_fixed, 1);
    min_variance_explained_sc = 0.01;
    max_sc_regressors = 5;
    dd_for_pca = dd_centered;
    col_variance = var(dd_for_pca, 0, 1);
    dd_for_pca = dd_for_pca(:, col_variance > 1e-10);
    
    if isempty(dd_for_pca) || size(dd_for_pca,2) == 0
        fprintf('Sujet %s: Pas de variance dans les canaux courts après nettoyage. Pas de régresseurs SC.\n', current_id_for_log);
        return;
    end

    final_sc_regressors = []; % Initialiser
    try
        [~, score, ~, ~, explained] = pca(dd_for_pca, 'Centered', false);
        
        % Logique de sélection des composantes (inchangée pour l'instant, mais vérifiez si c'est ce que vous voulez)
        explained_cumulative = cumsum(explained);
        num_components_to_keep_variance = find(explained_cumulative >= 90, 1, 'first'); % Expliquer 90%
        if isempty(num_components_to_keep_variance)
            num_components_to_keep_variance = size(score,2); 
        end
        num_components_to_keep = min(num_components_to_keep_variance, max_sc_regressors); 
        num_components_to_keep = min(num_components_to_keep, size(score,2));
        
        if num_components_to_keep > 0
           final_sc_regressors = score(:, 1:num_components_to_keep);
        end
        
    catch ME_pca
        warning('Sujet %s: PCA sur les canaux courts a échoué: %s. Pas de régresseurs SC.', current_id_for_log, ME_pca.message);
        % final_sc_regressors reste vide
    end

    if ~isempty(final_sc_regressors)
        for j = 1:size(final_sc_regressors, 2)
            st = nirs.design.StimulusVector;
            st.regressor_no_interest = true;
            st.name = ['SS_PCA' num2str(j)];
            st.time = data_out.time;
            st.vector = final_sc_regressors(:, j);
            data_out.stimulus(st.name) = st;  
        end
        sc_regressors_valid = true;
        fprintf('Sujet %s: %d régresseurs SC (PCA) ajoutés.\n', current_id_for_log, size(final_sc_regressors, 2));
    else
        fprintf('Sujet %s: Aucun régresseur SC valide généré.\n', current_id_for_log);
    end
end