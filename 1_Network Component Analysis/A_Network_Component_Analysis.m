%% code is adapted according to Liao et al. 2003


function A_Network_Component_Analysis

%constrain iterations of algorithm
options = optimoptions('lsqlin','MaxIter',2000);        %constrain iterations of algorithm
warning off all;




%% load gene expression data E and connectivity matrix A from supplementary tables

[FileName,PathName,FilterIndex]=uigetfile('*.xlsx','Select Supplementray Tables');
[status,sheets,~]= xlsfinfo([PathName FileName]);
[~,~,raw1] = xlsread([PathName FileName],9);
[~,~,raw2] = xlsread([PathName FileName],10);
A = cell2mat(raw1(3:end, 2:end));
E = cell2mat(raw2(5:end, 2:end));




%% generate boundaries for connections

lb_a = A;
ub_a = A;
for k = 1:size(A,1)   
    
    for l = 1:size(A,2)
        if A(k,l) == 1                          % bounderies for positiv interaction: realmin / 1
            lb_a(k,l) = realmin;    
            ub_a(k,l) = 1;
        end
        
        if A(k,l) == -1                         % bounderies for negative interaction: -1 / -realmin
            lb_a(k,l) = -1;
            ub_a(k,l) = -realmin;
        end
        
        if A(k,l) == 0                          % positions without interactions have no bounderies
            lb_a(k,l) = 0;
            ub_a(k,l) = 0;
        end
    end
end

lb_a = lb_a'; ub_a = ub_a';                     % transpose boundery matrices




%% calculate transcription factor activity

for ix = 1:100                                  % 100 randomized iterations 

    resnorms_AP = [];                           

    Ar = A;                                     % new randomized connectivity matrix with arbitrary values between -1 and 1
    R = rand(size(A));
    Ar = A.*R; 

    for i=1:size(E,2)                           % initial calculation of transcription factor activities P using the linear least square solver lsqlin and the randomized connection matrix Ar
        P_(:,i)=lsqlin(Ar,E(:,i),[],[],[],[],[],[],[], options);
    end

    E_2=E';      
    A_2=A';

    for j=1:10000                               % iterative calculation of A and P update
        disp(j)
        pause(2);
        P_2=P_';                   
   
    
    
        for i=1:size(E_2,2)                     %A update with boundries
            g=find(A_2(:,i));                   %check which positions of A should be calculated / in which positions we have a connection
            [A_3(g,i), resnorm]=lsqlin(P_2(:,g),E_2(:,i),[],[],[],[],lb_a(g,i),ub_a(g,i),[], options);      %calculation using lsqlin and boundaries for A
            resnorm_A(i)=resnorm;               %the squared 2-norm of the residual for each optimisation of A
        end
    
        resnorms_A(j)=sum(resnorm_A);    %sum of residuals of A for each optimisation   
    
        A_=A_3';

    
        for i=1:size(E,2)                           %P update
            [P_(:,i),resnorm, residual]=lsqlin(A_,E(:,i),[],[],[],[],[],[],[], options);
            resnorm_P(i)=resnorm;                   %squared 2-norm of the residual for each optimisation
            residual_P(i,:)=residual;               %the residual which gives difference between E and E_calc
        end
        
        
        resnorms_P(j)=sum(resnorm_P);                   %sum of residuals of P for each optimisation
        resnorms_AP(j)=resnorms_A(j)*resnorms_P(j);     %combination of the squared 2-norm of the residual sums of A and P
        disp(resnorms_AP)

        Ps{j}=P_;           %save each update P and A after each round
        As{j}=A_;  
    
    
        if j>2                                      %terminate the calculation if the combined combination of the squared 2-norm of the residual sums of A and P increases or changes less then 1 perent between two iterations 
             if 1-(resnorms_AP(j)/resnorms_AP(j-1))<0.01 | resnorms_P(j)>resnorms_P(j-1)
                 break
             end
        end

    end
    
    f=j;

    warning on all

    P_final=Ps{f};                                                  %get final P after termination of calculation
    A_final=As{f};                                                  %get final A after termination of calculation
    E_final=A_final*P_final;                                        %get final E_calc calculated from A and P
    residuals_final=residual_P;                                     %get final difference between E_calc and E

    residual_P_abs=abs(residual_P);
    E_abs=abs(E_final);
    rel_error=residual_P_abs'./E_abs*100;                           %get final relative difference between E_calc and E for each gene
    rel_error_gene=mean(rel_error');                                %average deviation of each gene
    out=quantile(rel_error_gene, 0.95);                             
    rel_error_gene(find(rel_error_gene>out))=[];                    %exclude highest five percent of deviation
    rel_mean_deviation_all=mean(rel_error_gene);                    %final average deviation of E_calc and E

    my_field1 = strcat('P_final',num2str(ix));
    my_field2 = strcat('A_final',num2str(ix));
    my_field4 = strcat('resnorms_AP',num2str(ix));
    my_field5 = strcat('residuals_final',num2str(ix));
    my_field6 = strcat('rel_mean_deviation_all',num2str(ix));

    NCA_struct.(my_field1)=P_final;                                 % save
    NCA_struct.(my_field2)=A_final;
    NCA_struct.(my_field4)=resnorms_AP;
    NCA_struct.(my_field5)=residuals_final;
    NCA_struct.(my_field6)=rel_mean_deviation_all;

end

save('NCA_struct_log10.mat', 'NCA_struct')

B_bootstrapping_NCA_95CI(NCA_struct)


