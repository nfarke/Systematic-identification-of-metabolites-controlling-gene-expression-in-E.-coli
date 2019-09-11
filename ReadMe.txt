Matlab code for the calculation of the transcription factor activities (folder 1) and correlation (folder 2).

Software : Matlab 2018b



1_Network Component Analysis

- Run A_Network_Component_Analysis.m 
- Load SupplementaryTables.xlsx

*********************Information****************************
A_Network_Component_Analysis.m 	- calculates the transcription factor activites in 100 runs with randomized starting points
B_bootstrapping_NCA_95CI.m 	- calculates the average transcription factor activity and the 95% - confidence intervals (CI)

Input: 
Connectivity matrix and log10 -   gene expression data imported from the supplementary tables.

Output: 
NCA_struct_log10.mat 		- log10 - transcription factor activities, connectivity strength, squared 2-norm of the residual, the residual, and the final average deviation from the real gene expression data (accuracy) for each of the 100 iterations
NCA_result_log10.mat 		- average log10 - transcription factor activities ('P_mean') of the 100 randomized calculations and the 95%-CI ('P_up' & 'P_lo')

Note:
With randomized starting points for each of the 100 calculations, final average transcription factor activities can slightly differ from the final results from the Supplementary Tables.




2_Kinetic Correlation


- Run Kinetic_Correlation.m 
- Load SupplementaryTables.xlsx

*********************Information****************************
Kinetic_Correlation.m 		- calculates the correlation coefficient of a activatory or inhibitory kinetic between the transcription factor activity and the metabolite level

Input: 
Relative metabolite abundances and and log10 - transcription factor activities imported from the supplementary tables.

Output:
corr_result.mat			- correlation coefficients for the kinetic fit 
					RSQ0 - without a shift between transcription factor activity and metabolit level
					RSQ1 - with a shift of -1 between transcription factor activity and metabolite level
					RSQ  - best correlation coefficient between transcription factor activity and metabolite level
					LAG  - the corresponding time shift for correlation coefficient 