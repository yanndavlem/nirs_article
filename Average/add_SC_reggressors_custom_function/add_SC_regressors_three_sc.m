function data_sc = add_SC_regressors_three_sc(data)

for i = 1:length(data)
    
%     length_idx=find(data(i).probe.link.ShortSeperation);
    sc_idx = [5, 6, 25, 26, 43, 44];
    dd=data(i).data(:,sc_idx);
    dd=dd-ones(size(dd,1),1)*mean(dd,1);
    dd=orth(dd);
    
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