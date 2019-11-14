pairStimulation.m
This script generates model responses to pairwise stimulations.
	Input: 	M0 vs M1 validation exp automated_mRNA only.xlsx
     		M0 vs M2 validation exp automated_mRNA only.xlsx
     		Stimuli chart.xlsx
     		modelODE.m
			modelParams.m
			PairstimuliScreening_3inputs.m
			y0.mat

	Output: simulation results/
			M0 vs M1 validation exp automated_mRNA only.xlsx_Screening_raw.txt
			M0 vs M1 validation exp automated_mRNA only.xlsx_Screening.txt
			M0 vs M2 validation exp automated_mRNA only.xlsx_Screening_raw.txt
			M0 vs M2 validation exp automated_mRNA only.xlsx_Screening.txt
			Screening_output.txt
			Screening_percentMatch.txt
			
screeningProfile.R
This script analyzes and plots model responses to pairwise stimulations and single input stimulations.
The output includes Fig 4(ABC), 5A, and supp Fig 1(ABC).
	Input: 	M0vsM1_matched.csv
			M0vsM1_matched.csv
			./simulation results/Screening_output.txt
			./simulation results/Screening_percentMatch.txt
	Output: simulation results/
			Screening_hclust1_stimuli-0.7nontemp12.csv - this is the stimulus combination clustering
			Screening_hclust2_nodes-0.7nontemp12.csv
			Signaling_Modules_for_Fig5A0.7nontemp12.csv - this is the clustering/node modules used in Fig 5A
			
			plots/
			Screening_module_dim2_IFNg0.7_12.png - Fig 5A
			Screening_PCA2D_dim0.7nontemp12.png - Fig 4B
			Screening_PCA2D_f_contri_dim0.7nontemp12.png - Fig 4C
			Screening_PCA2D_f_contri_single0.7nontemp12.png - supp Fig 1C
			Screening_PCA2D_single0.7nontemp12.png - supp Fig 1B
			Screening_single0.7nontemp12.png - supp Fig 1A
			Screening0.7nontemp12.png - Fig 4A

Coloring Node Modules for Fig5A.pdf and Coloring Node Module for Fig5A.png
	Model nodes were colored by their modules in Fig 5A. Only modules shown in Fig 5A were colored (MAPKs, MFkB, STAT1, STAT3, PI3K, and STAT6).
	Generated from macrophage_test_original6.3p.cys after importing Signaling_Modules_for_Fig5A0.7nontemp12.csv as a node attribute.

simulation results/Stimulus_Combination_clustering_manual_assignment_of_hclust1
	Based on Screening_hclust1_stimuli-0.7nontemp12.csv. Manually labeling each stimulus combination cluster to match the y axis labels of Fig 4A.