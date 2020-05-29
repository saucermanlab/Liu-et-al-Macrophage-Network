M1_M2_Validation.m
This script generates model responses to LPS+IFNg and IL4 stimulation.
The output includes Fig 2A and 2B.
	Input: 	mac validation sum M1.xlsx
			mac validation sum M2.xlsx
    		M0 vs M1 validation exp automated_mRNA only.xlsx
    		M0 vs M2 validation exp automated_mRNA only.xlsx
    		mac validation sum M1 PM.xlsx
    		mac validation sum M2 PM.xlsx
    		modelODE.m
			modelParams.m
			QuantValidation_3inputs.m
			y0.mat

	Output: simulation results/
			M0 vs M1 validation exp automated_mRNA onlyin0.7_raw.txt
			M0 vs M1 validation exp automated_mRNA onlyin0.7.txt
			M0 vs M2 validation exp automated_mRNA onlyin0.7_raw.txt
			M0 vs M2 validation exp automated_mRNA onlyin0.7.txt
			mac validation sum M1in0.7_raw.txt
			mac validation sum M1in0.7.txt
			mac validation sum M2in0.7_raw.txt
			mac validation sum M2in0.7.txt
			macmodelvalidation_M2in0.7.txt
			macmodelvalidation_M1in0.7.txt
			
			plots/
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

M1_M2_Validation.R
This script analyzes and plots model responses to LPS+IFNg and IL4 stimulation.
The output includes Fig 2C.
	Input: 	simulation results/
			M0 vs M1 validation exp automated_mRNA onlyin0.7.txt
			M0 vs M2 validation exp automated_mRNA onlyin0.7.txt
	Output: plots/
			Validation_both_dim_Seq_cscaled0.7142_alt_col.png - Fig 2C


IFNg_IL4_Validation.m
This script generates model responses to IFNg and IL4 combined stimulation.
	Input: 	IFNg validation GSE84520.xlsx
			IFNg+IL4 validation GSE84520.xlsx
    		M0 vs M2 validation Illum-GSE84520.xlsx
    		modelODE.m
			modelParams.m
			QuantValidation_3inputs.m
			y0.mat
	Output: simulation results/
			IFNg validation GSE84520_act.txt
			IFNg validation GSE84520_raw.txt
			IFNg validation GSE84520_validation.txt
			IFNg+IL4 validation GSE84520_act.txt
			IFNg+IL4 validation GSE84520_raw.txt
			IFNg+IL4 validation GSE84520_validation.txt
			M0 vs M2 validation Illum-GSE84520_act.txt
			M0 vs M2 validation Illum-GSE84520_raw.txt
			M0 vs M2 validation Illum-GSE84520_validation.txt

IFNg_IL4_Validation.R
This script analyzes and plots model responses to IFNg and IL4 combined stimulation.
The output includes Fig 5B and 5C.
	Input: 	simulation results/
			IFNg validation GSE84520_validation.txt
			IFNg+IL4 validation GSE84520_validation.txt
			M0 vs M2 validation Illum-GSE84520_validation.txt
	Output: plots/
			Validation_84520_142Arg1_mrna.png - Fig 5C
			Validation_84520_142IL4Ra_mrna.png - Fig 5C
			Validation_84520_142iNOS_mrna.png - Fig 5C
			Validation_84520_142SOCS1_mrna.png - Fig 5C
			Validation_84520_142TNFa_mrna.png
			Validation_84520_cscaled_142_reverse.png - Fig 5B