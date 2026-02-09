function data_sc = add_SC_regressors_PCA_close_SC(data,whichSC)

for i = 1:length(data)
    
    dd=data(i).data(:,whichSC);
    
    st=nirs.design.StimulusVector;
    st.regressor_no_interest=false;
    st.name=['SC_' num2str(whichSC) '_PCA'];
    st.time=data(i).time;
    st.vector=dd;
    data(i).stimulus(st.name)=st;
    
end

data_sc = data;





