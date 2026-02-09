classdef GLM_one_SC < nirs.modules.AbstractGLM
    %% GLM_close_SC - GLM wrapper avec sélection d'un short channel spécifique
    %
    % Options:
    %     basis       - a Dictionary object containing temporal bases using stim name as key
    %     verbose     - flag to display progress
    %     trend_func  - a function that takes in a time vector and returns trend regressors
    %     type        - {OLS, NIRS-SPM, or [AR-IRLS]} 
    %     AddShortSepRegressors - flag to include short separation data
    %     whichSC     - index of the specific short channel to use as regressor
    
    properties
        type;
        AddShortSepRegressors = false;
        whichSC = 1; % Index du short channel spécifique à utiliser
        useGPU = false;
        precisionSingle = false;
        options;
    end
    
    methods
        function obj = GLM_close_SC(prevJob)
            if nargin > 0, obj.prevJob = prevJob; end
            
            obj.name = 'GLM model with specific short channel';
            obj.basis('default') = nirs.design.basis.Canonical();
            obj.type = 'OLS';
        end
        
        function obj = set.type(obj, type)
            validtypes = {'OLS', 'NIRS-SPM', 'AR-IRLS', 'MV-GLM', 'Nonlinear'};
            if(~ismember(type, validtypes))
                disp('type must be one of : ')
                disp(strvcat(validtypes));
                return;
            else
                obj.type = type;
            end
            
            % use the call functions to evoke any special messages or conditions
            switch(obj.type)
                case('OLS')
                    j = nirs.modules.OLS();
                case('AR-IRLS')
                    j = nirs.modules.AR_IRLS();
                case('NIRS-SPM')
                    j = nirs.modules.NIRS_SPM_GLM();
                case('MV-GLM')
                    j = nirs.modules.MultiVarGLM();
                    disp(['Note: Inputs expected to be optical density']);
                case('Nonlinear')
                    j = nirs.modules.nonlin_GLM();
                    obj.basis = j.basis;
                otherwise
                    error('type not recognized');
            end
            
            obj.citation = j.citation;
            obj.options = [];
            pj = obj.prevJob;
            flds = fields(j);
            lst = find(ismember(flds, fields(obj)));
            for i = 1:length(lst)
                obj.(flds{lst(i)}) = j.(flds{lst(i)});
            end
            lst = find(~ismember(flds, fields(obj)));
            for i = 1:length(lst)
                obj.options = setfield(obj.options, flds{lst(i)}, j.(flds{lst(i)}));
            end
            obj.prevJob = pj;
        end
        
        function S = runThis(obj, data)
            % Créer une copie des données pour préserver les originales
            data_with_sc = data;
            
            % Vérifier si on utilise les short channels comme régresseurs
            if (obj.AddShortSepRegressors)
                for i = 1:numel(data)
                    if (~nirs.util.hasshortdistances(data(i)))
                        continue;
                    end
                    
                    % Identifier tous les short channels
                    lstss = find(data(i).probe.link.ShortSeperation);
                    
                    % Si un short channel spécifique est demandé, conserver uniquement celui-ci
                    if obj.whichSC > 0 && obj.whichSC <= length(lstss)
                        % Pour chaque court circuit présent dans les données
                        if obj.whichSC <= length(lstss)
                            % Récupérer les données du short channel spécifié
                            dd = data(i).data(:, lstss(obj.whichSC));
                            
                            % Normaliser le signal
                            dd = dd - mean(dd);
                            
                            % Traiter les NaN éventuels
                            tmp = nirs.core.Data;
                            tmp.time = data(i).time;
                            tmp.data = dd;
                            j = nirs.modules.FixNaNs;
                            tmp = j.run(tmp);
                            dd = tmp.data;
                            
                            % Créer le régresseur avec ce signal
                            st = nirs.design.StimulusVector;
                            st.regressor_no_interest = true;
                            st.name = ['SS_SC' num2str(obj.whichSC)];
                            st.time = data(i).time;
                            st.vector = dd;
                            st.vector = st.vector - mean(st.vector);
                            st.vector = st.vector ./ sqrt(var(st.vector));
                            data_with_sc(i).stimulus(st.name) = st;
                            
                            % Créer un second régresseur avec la dérivée première
                            st = nirs.design.StimulusVector;
                            st.regressor_no_interest = true;
                            st.name = ['SS_SC' num2str(obj.whichSC) '_diff1'];
                            st.time = data(i).time;
                            st.vector = [diff(dd); 0];
                            st.vector = st.vector - mean(st.vector);
                            st.vector = st.vector ./ sqrt(var(st.vector));
                            data_with_sc(i).stimulus(st.name) = st;
                        end
                    end
                end
            end
            
            % Utiliser l'implémentation standard du GLM avec les données préparées
            switch(obj.type)
                case('OLS')
                    j = nirs.modules.OLS();
                case('AR-IRLS')
                    j = nirs.modules.AR_IRLS();
                case('NIRS-SPM')
                    j = nirs.modules.NIRS_SPM_GLM();
                case('MV-GLM')
                    j = nirs.modules.MultiVarGLM();
                case('Nonlinear')
                    j = nirs.modules.nonlin_GLM();
                otherwise
                    error('type not recognized');
            end
            
            j.basis = obj.basis;
            j.verbose = obj.verbose;
            j.trend_func = obj.trend_func;
            j.goforit = obj.goforit;
            try; j.useGPU = obj.useGPU; end
            try; j.precisionSingle = obj.precisionSingle; end
            
            if(~isempty(obj.options))
                flds = fields(obj.options);
                for i = 1:length(flds)
                    j.(flds{i}) = obj.options.(flds{i});
                end
            end
            
            % Exécuter le GLM standard avec nos régresseurs personnalisés
            S = j.run(data_with_sc);
            
            % Nettoyer les résultats pour ne garder que les conditions d'intérêt
            for idx = 1:length(S)
                [~, Stim] = nirs.design.createDesignMatrix(data(idx).stimulus, data(idx).time, obj.basis);
                StimNew = unique(nirs.getStimNames(S(idx)));
                
                lst = []; lst2 = [];
                for j = 1:length(Stim)
                    ss = Stim{j};
                    if(~isempty(strfind(ss, ':')))
                        ss(strfind(ss, ':'):end) = [];
                    end
                    if(ss(1) == 'x' && ~ismember(ss, data(idx).stimulus.keys))
                        ss(1) = [];
                    end
                    st = data(idx).stimulus(ss);
                    if(~ismember('regressor_no_interest', fields(st)))
                        lst = [lst j];
                    else
                        if(~st.regressor_no_interest)
                            lst = [lst j];
                        else
                            lst2 = [lst2 j];
                        end
                    end
                end
                StimRm = {Stim{lst2}};
                Stim = {Stim{lst}};
                
                SS = {};
                for i = 1:length(StimNew)
                    for j = 1:length(Stim)
                        if(~isempty(ismember(StimNew{i}, Stim{j}, 'rows')) && ...
                                isempty(strfind(StimNew{i}, 'SS_SC')) && ...
                                ~ismember(StimNew{i}, StimRm))
                            SS{end+1} = StimNew{i};
                        end
                    end
                end
                SS = unique(SS);
                for ii = 1:length(SS)
                    SS{ii} = ['^' SS{ii} '$'];
                end
                
                % Garder seulement les stimuli d'intérêt dans les résultats
                j = nirs.modules.KeepStims;
                j.listOfStims = SS;
                j.regex = true;
                S(idx) = j.run(S(idx));
            end
        end
    end
end