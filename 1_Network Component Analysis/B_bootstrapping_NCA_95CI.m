function B_bootstrapping_NCA_95CI(NCA_struct)

% get only the P matrices and arrange them in a new structure
for ix=1:numel(fieldnames(NCA_struct))/5
    my_field1 = strcat('P_final',num2str(ix));
    P{ix}=NCA_struct.(my_field1);
end


for i=1:size(P{1},2) %go through all timepoints / conditions in each P
    disp(i)
        for j=1:size(P,2)   % take out from each P matrix the corresponding timepoint
            P_time(:,j)=P{j}(:,i);  
        end   
        
        for k=1:size(P_time,1)  % for each transcription factor, calculate the mean and 95%-CI across all computed randomised calculation
            x=P_time(k,:);
            P_mean(k,i)=mean(x);
            SEM = std(x)/sqrt(length(x));               % Standard Error
            ts = tinv([0.025  0.975],length(x)-1);      % T-Score
            CI = mean(x) + ts*SEM;                      % Confidence Intervals
            P_lo(k,i)=CI(1);
            P_up(k,i)=CI(2);
        end
end
    
save('NCA_result_log10.mat', 'P_mean', 'P_lo', 'P_up')
    

