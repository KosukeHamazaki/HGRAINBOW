# HGRAINBOW
#### Supplementary data and scripts used in the article [Kosuke Hamazaki and Hiroyoshi Iwata, “RAINBOW: Haplotype-based genome-wide association study using a novel SNP-set method”, PLOS Computational Biology, 16(2): e1007663, 2020.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007663)
#### Author : Kosuke Hamazaki (hamazaki@ut-biomet.org)
#### Date : 2019/03/28
Here, I explain the structure of this repository.

----------

- `HGRAINBOW`
	- `data.zip` : This folder contains datasets analyzed in this study. Please download and then decompress.
		- `extra` : This folder contain datasets other than genotype and phenotype.
			- `Table_S1_origin.csv` : The original data for `Table S1`.
			- `additive_relationship_matrix.RData` : The additive genetic relationship matrix used in this study.
			-  `haplotype_block_list.csv` : The list of haplotype block information used in the SNP-set GWAS methods.
			-  `map2.csv` : The physical map corresponding to the haplotype block list.
		-  `genotype` : The marker genotype used in this study.
			-  `L3024_core_extract_L414_ind1A_ind1B_MAF_0.025_geno.tsv` : The marker genotype of 414 accessions used in this study.
			-  `L3024_core_extract_L414_ind1A_ind1B_MAF_0.025_haplo1.tsv` : The marker haplotype of 414 accessions used in this study.
			-  `L3024_core_extract_L414_ind1A_ind1B_MAF_0.025_haplo2.tsv` : The marker haplotype of 414 accessions used in this study.
		-  `phenotype` : The phenotypic data simulated in this study.
			-  `number_of_causals=2_direction_of_effect=minus_trial_no=1.csv` : The phenotypic data for the repulsion scenario.
			-  `number_of_causals=2_direction_of_effect=plus_trial_no=1.csv` : The phenotypic data for the coupling scenario.
			-  `seeds_number_of_causals=2_direction_of_effect=minus_trial_no=1.csv` : The random seeds used for simulating phenotypic data for the repulsion scenario.
			-  `seeds_number_of_causals=2_direction_of_effect=plus_trial_no=1.csv` : The random seeds used for simulating phenotypic data for the coupling scenario.
	-  `scripts` : This folder contains scripts with the `R` language used in this study.
		-  `0.0_Rice_HGRAINBOW_subpop_list_to_generate_haplotype_data.R` : Extract subpopulation list from 3,000 accessions.
		-  `0.1_Rice_HGRAINBOW_haplotype_data.txt`  : Extract haplotype data of 414 accessions (`plink` and `vcftools`) adn estimate haplotype blocks by `plink` (not `R`!).
		-  `0.2_Rice_HGRAINBOW_modyfing_haplotype_block_list.R` : Modify haplotype block data estimated by `plink` into the format as `data/extra/haplotype_block_list.csv`.
		-  `0.3_Rice_HGRAINBOW_Simulation_of_phenotypic_values_and_some_preparation.R` : Simulate phenotypic values for both scenarios, coupling and repulsion.
		-  `0.4_Rice_HGRAINBOW_subpop_list_for_Table_S1.R` : Scripts for generating `Table_S1_origin.csv`.
		-  `1.1_Rice_HGRAINBOW_score_SKAT_geneset.R` : The function to run `SKAT` as haplotype-based GWAS method.
		-  `1.2_Rice_HGRAINBOW_haplotype_group_fixed_GWAS.R` : The function to run `HGF` as haplotype-based GWAS method.
		-  `1.3_Rice_HGRAINBOW_SS_gwas.R` : The function to calculate summary statistics from GWAS results.
		-  `2.1_Rice_HGRAINBOW_Haplotype_based_GWAS_for_RAINBOW_package_paper.R` : Perform each GWAS method and save the results.
		-  `3.1_Rice_HGRAINBOW_Summary_of_HGRAINBOW_Haplotype_based_GWAS_for_RAINBOW.R` : Summary GWAS results.
		-  `3.2_Rice_HGRAINBOW_Manhattan_plot_from_the_results.R` : Draw Manhattan plots from GWAS results.
		-  `3.3_Rice_HGRAINBOW_Summary_plot_of_HGRAINBOW_Haplotype_based_GWAS_for_RAINBOW.R` : Draw figures in the article from GWAS results.
		-  `3.4_Rice_HGRAINBOW_Manhattan_plot_for_overwhelming_results.R` : Draw Manhattan plots for overwhelming results (`Figure 4`).
