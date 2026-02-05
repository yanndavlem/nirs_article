    
function data_sc = add_SC_regressors_one_sc(data,whichSC)
% Cette fonction ajoute les signaux bruts de 6 canaux courts spécifiques 
% en tant que régresseurs de nuisance individuels.
for i = 1:length(data)
    
    dd=data(i).data(:,whichSC);
    
    for sc = 1:length(whichSC)
        st=nirs.design.StimulusVector;
        st.regressor_no_interest=false;
        st.name=['SC_' num2str(sc) '_PCA'];
        st.time=data(i).time;
        st.vector=dd(:,sc);
        data(i).stimulus(st.name)=st;    
    end        
end

data_sc = data;



