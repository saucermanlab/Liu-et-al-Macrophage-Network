* Please make sure the working directory is set to the correct path (folder Liu-et-al-Macrophage-Network) before running the codes.

./Model Simulation/M1_M2_Validation.m
This script generates model responses to LPS+IFNg and IL4 stimulation.
The output includes Fig 2A and 2B.
	Input: ./Input/
			mac validation sum M1.xlsx
			mac validation sum M2.xlsx
    			LPS+IFNg validation RNASeq.xlsx
    			IL4 validation RNASeq.xlsx
    			mac validation sum M1 PM.xlsx
    			mac validation sum M2 PM.xlsx
    			modelODE.m
			modelParams.m
			y0.mat
		./Model Simulation/QuantValidation_3inputs.m

	Output: ./Model Simulation/simulation results/
			LPS+IFNg validation RNASeq.txt
			LPS+IFNg validation RNASeq.txt
			IL4 validation RNASeq_raw.txt
			IL4 validation RNASeq.txt
			mac validation sum M1in0.7_raw.txt
			mac validation sum M1in0.7.txt
			mac validation sum M2in0.7_raw.txt
			mac validation sum M2in0.7.txt
			macmodelvalidation_M2in0.7.txt
			macmodelvalidation_M1in0.7.txt
			
		./Model Simulation/plots/
			Arg1_mrna0.7.tif - Fig 2B
			CCL17_mrna0.7.tif
			IKBa_mrna0.7.tif
			IL1_mrna0.7.tif - Fig 2B
			IL6_mrna0.7.tif
			IL10_mrna0.7.tif
			IL12_mrna0.7.tif
			iNOS_mrna0.7.tif - Fig 2B
			IRF1_mrna0.7.tif
			macmodelvalidation M1 in0.7.tif - Fig 2A
			macmodelvalidation M2 in0.7.tif - Fig 2A
			MMP3_mrna0.7.tif
			MMP7_mrna0.7.tif
			MMP9_mrna0.7.tif
			MMP12_mrna0.7.tif
			PPARg_mrna0.7.tif
			SOCS1_mrna0.7.tif - Fig 2B
			SOCS3_mrna0.7.tif
			TNFa_mrna0.7.tif
			VEGF_mrna0.7.tif

./Model Simulation/M1_M2_Validation.R
This script analyzes and plots model responses to LPS+IFNg and IL4 stimulation.
The output includes Fig 2C.
	Input: 	./Model Simulation/simulation results/
			LPS+IFNg validation RNASeq.txt
			IL4 validation RNASeq.txt

	Output: ./Model Simulation/plots/
			Validation_both_dim_Seq_cscaled0.7142_alt_col.png - Fig 2C


./Model Simulation/IFNg_IL4_Validation.m
This script generates model responses to IFNg and IL4 combined stimulation.
	Input: 	./Input/
			IFNg validation GSE84520.xlsx
			IFNg+IL4 validation GSE84520.xlsx
    			IL4 validation Illum-GSE84520.xlsx
    			modelODE.m
			modelParams.m
			y0.mat
		./Model Simulation/QuantValidation_3inputs.m

	Output: ./Model Simulation/simulation results/
			IFNg validation GSE84520_act.txt
			IFNg validation GSE84520_raw.txt
			IFNg validation GSE84520_validation.txt
			IFNg+IL4 validation GSE84520_act.txt
			IFNg+IL4 validation GSE84520_raw.txt
			IFNg+IL4 validation GSE84520_validation.txt
			IL4 validation Illum-GSE84520_act.txt
			IL4 validation Illum-GSE84520_raw.txt
			IL4 validation Illum-GSE84520_validation.txt

./Model Simulation/IFNg_IL4_Validation.R
This script analyzes and plots model responses to IFNg and IL4 combined stimulation.
The output includes Fig 5B and 5C.
	Input: 	./Model Simulation/simulation results/
			IFNg validation GSE84520_validation.txt
			IFNg+IL4 validation GSE84520_validation.txt
			IL4 validation Illum-GSE84520_validation.txt
	Output: ./Model Simulation/plots/
			Validation_84520_142Arg1_mrna.png - Fig 5C
			Validation_84520_142IL4Ra_mrna.png - Fig 5C
			Validation_84520_142iNOS_mrna.png - Fig 5C
			Validation_84520_142SOCS1_mrna.png - Fig 5C
			Validation_84520_142TNFa_mrna.png
			Validation_84520_cscaled_142_reverse.png - Fig 5B

./Pairwise Simulation/pairStimulation.m
This script generates model responses to pairwise stimulations.
	Input: 	./Input/
			LPS+IFNg validation RNASeq.xlsx
    			IL4 validation RNASeq.xlsx
     			Stimuli chart.xlsx
     			modelODE.m
			modelParams.m
			y0.mat
		./Pairwise Simulation/PairstimuliScreening_3inputs.m

	Output: ./Pairwise Simulation/simulation results/
			LPS+IFNg validation exp automated_mRNA only.xlsx_Screening_raw.txt
			LPS+IFNg validation exp automated_mRNA only.xlsx_Screening.txt
			IL4 validation exp automated_mRNA only.xlsx_Screening_raw.txt
			IL4 validation exp automated_mRNA only.xlsx_Screening.txt
			Screening_output.txt
			Screening_percentMatch.txt
			
./Pairwise Simulation/screeningProfile.R
This script analyzes and plots model responses to pairwise stimulations and single input stimulations.
The output includes Fig 4(ABC), 5A, and supp Fig 1(ABC).
	Input: 	./Input/
			M0vsM1_matched.csv
			M0vsM1_matched.csv
		./Pairwise Simulation/simulation results/Screening_output.txt
		./Pairwise Simulation/simulation results/Screening_percentMatch.txt
	Output: ./Pairwise Simulation/simulation results/
			Screening_hclust1_stimuli-0.7nontemp12.csv - this is the stimulus combination clustering
			Screening_hclust2_nodes-0.7nontemp12.csv
			Signaling_Modules_for_Fig5A0.7nontemp12.csv - this is the clustering/node modules used in Fig 5A
			
		./Pairwise Simulation/plots/
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

./Sensitivity Screening/runSensScreening.m
This script runs sensitivity analysis (node KO) with all single stimuli and all stimulus combinations.
	Input: 	./Input/
			Stimuli chart.xlsx
     			modelODE.m
			modelParams.m
			y0.mat
		./Sensitivity Screening/sensAnalysisScr.m

	Output: ./Sensitivity Screening/simulation results/
			macmodelSens_IFNb+_0.7.txt			macmodelSens_IFNb+IL1_0.7.txt			macmodelSens_IFNb+IL4_0.7.txt			бн			macmodelSens_TNFa+IL6_0.7.txt			macmodelSens_TNFa+IL10_0.7.txt			macmodelSens_TNFa+IL12_0.7.txt


./Sensitivity_M1M2_4h/runSens.m
This script plots sensitivity analysis (node KO) results with simulated M1 and M2 stimulations.
	Input: ./Input/
			macmodelSensLPS+IFNg.txt
			macmodelSensIL4.txt
     			modelODE.m
			modelParams.m
			y0.mat
		./Sensitivity_M1M2_4h/
			sensAnalysis.m
			xticklabel_rotate.m

	Output: ./Sensitivity Screening/
			Sensitive Matrix LPS+IFNg.tiff
			Sensitive Matrix LPS+IFNg Dim.tiff - Fig 3B
			Sensitive Matrix IL4.tiff
			Sensitive Matrix IL4 Dim.tiff - Fig 3C
			Sens_Accum_M1-M2.tiff - Fig 3A
