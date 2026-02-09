function data_sc = add_SC_regressors_no_PCA(data)

for i = 1:length(data)
    
    sc_idx=find(data(i).probe.link.ShortSeperation);
    dd=data(i).data(:,sc_idx);              
    dd=dd-ones(size(dd,1),1)*mean(dd,1);
    
    for sc = 1:length(sc_idx)
        st=nirs.design.StimulusVector;
        st.regressor_no_interest=false;
        st.name=['SC_' num2str(sc) '_PCA'];
        st.time=data(i).time;
        st.vector=dd(:,sc);
        data(i).stimulus(st.name)=st;    
    end        
end

data_sc = data;





