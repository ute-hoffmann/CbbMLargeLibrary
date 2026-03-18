# The fitness landscape of a Form II rubisco in a photosynthetic bacterium guides engineering of oxygen tolerance

Here, code and info for the analyses presented in the manuscript "The fitness landscape of a Form II rubisco in a photosynthetic bacterium guides engineering of oxygen tolerance" are collected.

[RubisCO_LargeLib_plasmid.dna](RubisCO_LargeLib_plasmid.dna) is the plasmid used for integrating expression cassette for all CbbM variants genomically.

## Overview data to Figure and Tables

Data to some of the Figures and Tables shown in Hoffmann & Knave et al. are available in this repository. The following table gives an overview:

| Figure  | Description | Code  | Data  |
| ------------- | ------------- | ------------- |------------- |
| Fig. 1C  | Correlation between weighted mean fitness under different growth conditions  | [.Rmd file with analyses](code/EDA.Rmd)  | [File 1 as pdf](EDA_output/pdf/wfitness_CLN2_CLO2_withoutSignif_wo-labels.pdf), [File 2 as pdf](EDA_output/pdf/wfitness_CLO2_LD_withoutSignif_withoutLabeling.pdf)  |
| Fig. 1D  | Density plots of saturational/combinatorial library at different growth conditions  | [.Rmd file with analyses](code/EDA.Rmd)  | [File as pdf](EDA_output/pdf/density_compareCombiSatur.pdf)  |
| Fig. 2A  | Heat map of DMS effect  | [.Rmd file with analyses](code/EDA.Rmd), [python visualization file 1](EDA_output/heatMap/matrix_plot_CLN2.py), [file 2](EDA_output/heatMap/matrix_plot_CLO2.py), [file 3](EDA_output/heatMap/matrix_plot_LD.py)  | [File 1](EDA_output/heatMap/subst_matrix_CLN2.pdf), [file 2](EDA_output/heatMap/subst_matrix_CLO2.pdf), [file 3](EDA_output/heatMap/subst_matrix_LD.pdf)  |
| Fig. 2B  | Bar plot with fitness per number of exchanges  | [.Rmd file with analyses](code/EDA.Rmd)  | [File as pdf](EDA_output/pdf/bar_fitness_data_numberMuts_onlyExpected.pdf)  |
| Fig. 2C  | Percent of Higher-order variants with positive/neutral/negative effect  | [.Rmd file with analyses](code/EDA.Rmd) and [epistasis.py](code/epistasis.py)  | [File as pdf](EDA_output/epistasis/Positive_negativeExchanges_higherOrdervariants.pdf)  |
| Fig. 2D  | Growth curves extracted from pool for epistatic effect  | [.Rmd file with analyses](code/EDA.Rmd)  | Plot labelled "Epistatic effects and 'strain E'" in [.Rmd file with analyses](code/EDA.Rmd) |
| Fig. 3A  | Scatter plot normalized weighted mean fitness to EVmutation score  | [.Rmd file with analyses](code/EDA.Rmd)  | [File as pdf](EDA_output/pdf/scatter_norm_all_EVcoup.pdf)  |
| Fig. 3C  | Scatter plot normalized weighted mean fitness to additive score  | [.Rmd file with analyses](code/EDA.Rmd)  | [File as pdf](EDA_output/pdf/scatter_additive_combi.pdf)  |
| Supp. Fig. 2  | Growth curves highlighting K200R/K214R and base variant   | [.Rmd file with analyses](code/EDA.Rmd)  | [File as pdf](EDA_output/pdf/all_sgRNAtargets_timeLinePlot.pdf)  |
| Supp Fig. 3A  | Scatter plot normalized weighted mean fitness to relative solvent accessibility  | [.Rmd file with analyses](code/EDA.Rmd)  | [File as pdf](EDA_output/pdf/relSAS_norm.pdf)  |
| Supp. Fig. 3B  | Scatter plot normalized weighted mean fitness to conservation  | [.Rmd file with analyses](code/EDA.Rmd)  | [File as pdf](EDA_output/pdf/conservation_norm.pdf)  |
| Supp. Fig. 4A  | Volcano plot proteomics data generation 4/generation 0  | [.Rmd file for MS data](MS_firstCultivation/plot_MS_data.Rmd)  | [Fig. as pdf](MS_firstCultivation/Figures/Volcano_highCO2_gen0-4_vers2.pdf)  |
| Supp. Fig. 4B | Volcano plot proteomics data generation 8/generation 0  | [.Rmd file for MS data](MS_firstCultivation/plot_MS_data.Rmd)  | [Fig. as pdf](Volcano_highCO2_gen0-8_vers2.pdf)  |
| Supp. Fig. 7A  | Correlation EVcouplings to DeepSequence  | [.Rmd file with analyses](code/EDA.Rmd)  | [File as pdf](EDA_output/pdf/scatter_DeepSeq_EVcoup.pdf)   |
| Supp. Fig. 7B  | Correlation EVcouplings to MSA Transformer  | [.Rmd file with analyses](code/EDA.Rmd)  | [File as pdf](EDA_output/pdf/scatter_MSA_Transform_EVcoup.pdf)  |
| Supp. Fig. 8A  | Scatter plot normalized weighted mean fitness to EVmutation score  | [.Rmd file with analyses](code/EDA.Rmd)  | [File as pdf](EDA_output/pdf/scatter_norm_all_EVcoup.pdf)  |
| Supp. Fig. 8B  | Scatter plot normalized weighted mean fitness to DeepSequence  | [.Rmd file with analyses](code/EDA.Rmd)  | [File as pdf](EDA_output/pdf/scatter_DeepSeq_norm.pdf)  |
| Supp. Fig. 8C  | Scatter plot normalized weighted mean fitness to MSA Transformer  | [.Rmd file with analyses](code/EDA.Rmd)  | [File as pdf](EDA_output/pdf/scatter_MSA_Transform_norm.pdf)  |
| Supp. Fig. 9A  | Scatter plot normalized weighted mean fitness to EVmutation score without V323  | [.Rmd file with analyses](code/EDA.Rmd)  | [File as pdf](EDA_output/pdf/scatter_norm_all_EVcoup_noV323.pdf)  |
| Supp. Fig. 9B  | Scatter plot normalized weighted mean fitness to DeepSequence without V323  | [.Rmd file with analyses](code/EDA.Rmd)  | [File as pdf](EDA_output/pdf/scatter_DeepSeq_norm_noV323.pdf)  |
| Supp. Fig. 9C  | Scatter plot normalized weighted mean fitness to MSA Transformer without V323  | [.Rmd file with analyses](code/EDA.Rmd)  | [File as pdf](EDA_output/pdf/scatter_MSA_Transform_norm_noV323.pdf)  |
| Supp. Fig. 10A  | Scatter plot normalized weighted mean fitness to additive score  | [.Rmd file with analyses](code/EDA.Rmd)  | [File as pdf](EDA_output/pdf/scatter_additive_combi.pdf)  |
| Supp. Fig. 10B  | Scatter plot normalized weighted mean fitness to ProteinNPT predictions, with V323  | [.Rmd file with analyses](code/EDA.Rmd)  | [File as pdf](EDA_output/pdf/scatter_proteinNPT_norm.pdf)  |
| Supp. Fig. 10C  | Scatter plot normalized weighted mean fitness to ProteinNPT predictions, without V323  | [.Rmd file with analyses](code/EDA.Rmd)  | [File as pdf](EDA_output/pdf/scatter_proteinNPT_norm_noV323.pdf)  |

| Table  | Description | Code  | Data  |
| ------------- | ------------- | ------------- |------------- |
| Supp. Table 1  | Positions used in the DMS part of the library  | see below (Library design)  |  |
| Supp. Table 2  | Overview of library composition  | [.Rmd file with analyses](code/EDA.Rmd)  |   |
| Supp. Table 3  | Overview of variants with higher fitness value than base variant in different conditions  | [.Rmd file with analyses](code/EDA.Rmd)  |   |
| Supp. Table 4  | Proteomics data after 4 and 8 generations  | [.Rmd file for MS data](MS_firstCultivation/plot_MS_data.Rmd)  | [File as csv](MS_firstCultivation/tables/CCM_table_forPaper.csv), data for all proteins in MS comparisons in [same directory](MS_firstCultivation/tables)  |
| Supp. Table 5  | Variants of CbbMbase with the highest normalized fitness score  | [.Rmd file with analyses](code/EDA.Rmd)  | [File as csv](EDA_output/goodVariants.csv)  |

## Library design, zero-shot predictions, ProteinNPT and MoCHI

Library design was performed as published in [doi.org/10.1021/acssynbio.5c00065](https://doi.org/10.1021/acssynbio.5c00065) using code available in the [EVmut_inSilico repository](https://github.com/ute-hoffmann/EVmut_inSilico).

Files are collected in [library_design](library_design).

The [EVcouplings web server](https://v2.evcouplings.org/) was run on a version of CbbM with two amino acid exchanges (as detailed in the corresponding manuscript) and including the N-terminal Strep tag that was included for experiments. A bitscore of 0.3 was chosen. [This model](library_design/I_ab368e89ebe9421898840b5c014db5b3_TARGET_b0.3.model) was obtained and used for further analyses and predictions (DeepSequence, EVmutation).

### DeepSequence

#### Setting up conda environment

DeepSequence predictions were run on a Google Cloud Product n1-highmem-2 virtual machine with a NVIDIA T4 GPU. To set up the conda environment used for training, a [.yml file](https://github.com/fhalab/MLDE/blob/main/deep_sequence.yml) from [MLDE](https://github.com/fhalab/MLDE) was used.

```
conda env create -f deep_sequence.yml
```

#### Code adjustment DeepSequence

We adjusted the code available for [DeepSequence on GitHub](https://github.com/debbiemarkslab/DeepSequence) to ensure a correct assignment of residue numbers by exchanging [the following line](https://github.com/debbiemarkslab/DeepSequence/blob/fd5c272b273385c0b87fce9dec410a9d74463393/DeepSequence/helper.py#L249)

```
mut_seq[self.uniprot_focus_col_to_focus_idx[pos]]
```

by

```
mut_seq[self.focus_cols.index(self.uniprot_focus_col_to_focus_idx[pos])]
```

#### *in silico* evolution

Code from [this repository](https://github.com/chloechsu/combining-evolutionary-and-assay-labelled-data) was used to train the DeepSequence VAE ([train_vae.sh](https://github.com/chloechsu/combining-evolutionary-and-assay-labelled-data/blob/main/scripts/train_vae.sh)) on [this alignment](library_design/TARGET_b0.3_I_.a2m). All code used is collected in [library_design/code_DeepSeq](library_design/code_DeepSeq).

To run predictions, different bash wrapper scripts were used. For prediction of single sites, [1_single-site_mutagenesis_GallionellaRubisco_I.sh](library_design/code_DeepSeq/1_single-site_mutagenesis_GallionellaRubisco_I.sh) was used. Then, for *in silico* evolution, [2-1_complete_inSilico_Evolution_I.sh](library_design/code_DeepSeq/2-1_complete_inSilico_Evolution_I.sh) was used, followed by python scripts [2-2_extract_best_in_trajectory_I.py](library_design/code_DeepSeq/2-2_extract_best_in_trajectory_I.py), [2-3_extract_point-mutations_number-mutations_I.py](library_design/code_DeepSeq/2-3_extract_point-mutations_number-mutations_I.py) and [2-4_count-mutations_I.py](library_design/code_DeepSeq/2-4_count-mutations_I.py). After deciding on a subset of positions of interest, wrapper script [3_inSilico_Evolution_I_selectedPositions.sh](library_design/code_DeepSeq/3_inSilico_Evolution_I_selectedPositions.sh) was run.

Files needed as input for the different scripts are collected in [library_design/input](library_design/input).

### Predict fitness scores for all variants

Scripts for predicting fitness scores for combinatorial library members are collected in the folder [insilico_predictions](insilico_predictions).

#### EVcouplings

in EVcouplings_predict, use evmut_insilico.yml to create conda environment

```
conda activate evmut_insilico
python scripts/predict_allMembers_largeLibrary.py
```

#### DeepSequence

DeepSequence was run on a Google Cloud Product n1-highmem-2 virtual machine with a NVIDIA T4 GPU. The VAE was trained as described above. To run predictions, use [bash wrapper script](library_design/code_DeepSeq/4_Gallionella-I_calc-delta-E-specific-variants.sh) for python code [DeepSequence_predict_allMembers_largeLibrary_GallionellaI.py](library_design/code_DeepSeq/DeepSequence_predict_allMembers_largeLibrary_GallionellaI.py):

```
conda activate deep_sequence
bash 4_Gallionella-I_calc-delta-E-specific-variants.sh
```

#### MSA Transfomer and ProteinNPT

For both predictions, code from the [ProteinNPT repository](https://github.com/OATML-Markslab/ProteinNPT) was used, specifically the [commit from Dec 23rd, 2024](https://github.com/OATML-Markslab/ProteinNPT/commit/abe963b6afa15e76e8e6b37f01049047a27447c7). For MSA Transformer, the ProteinNPT pipeline runs an ensemble of five MSATransformer instances with different initation seeds.

### ProteinNPT predictions

The [commit from Dec 23rd, 2024](https://github.com/OATML-Markslab/ProteinNPT/commit/abe963b6afa15e76e8e6b37f01049047a27447c7) was used for all work with ProteinNPT. Used an adaptation of [pipeline.sh](https://github.com/OATML-Markslab/ProteinNPT/blob/master/scripts/pipeline.sh), specifically [pipeline_predictNewVariants.sh](ProteinNPT/pipeline_predictNewVariants.sh). A file was created including all experimental data from the two sequencing runs and a set of variants that seemed interesting based on high fitness values in experimental screen. Code for creating this file is collected in [ProteinNPT/2025-02-04_create_testSet](ProteinNPT/2025-02-04_create_testSet). Output (e.g. [stderr](ProteinNPT/2025-02-10_predictions_new_variants/2025-02-08_pipelineNewVariants_stderr.txt) and [stout](ProteinNPT/2025-02-10_predictions_new_variants/2025-02-08_pipelineNewVariants_stdout.txt)) of the pipeline is collected in [ProteinNPT/2025-02-10_predictions_new_variants](ProteinNPT/2025-02-10_predictions_new_variants).

[list_variants.txt](ProteinNPT/2025-02-04_create_testSet/list_variants.txt) includes single amino acid exchanges and combinatorial variants, which were combined in every possible manner using script [1_create_file_with_variants.py](ProteinNPT/2025-02-04_create_testSet/1_create_file_with_variants.py), yielding [list_combinations.txt](ProteinNPT/2025-02-04_create_testSet/list_combinations.txt). [Script 2_create_mutant_sequences_fromNewVariants.py](ProteinNPT/2025-02-04_create_testSet/2_create_mutant_sequences_fromNewVariants.py) was used to create [full-length sequences of these variants](ProteinNPT/2025-02-04_create_testSet/list_new_sequences.txt), which were combined with the list of experimental data for the measured library using (3_create_file_with_newVariants_compatible_with_trainingFile.py)[ProteinNPT/2025-02-04_create_testSet/3_create_file_with_newVariants_compatible_with_trainingFile.py]. Importantly, duplicated sequences were removed, yielding table [2025-02-09_predict_new_variants.csv](ProteinNPT/2025-02-04_create_testSet/2025-02-09_predict_new_variants.csv), which was used as input for [pipeline_predictNewVariants.sh](ProteinNPT/2025-02-10_predictions_new_variants/pipeline_predictNewVariants.sh).

### MoCHI predictions

All files related to MoCHI predictions are collected in [EDA_output/MoCHI](EDA_output/MoCHI). Output from EDA.Rmd was used as input for a python script [Create_training_set.py](EDA_output/MoCHI/Create_training_set.py), preparing tables used as input for MoCHI. MoCHI was run in a [mamba environment](EDA_output/MoCHI/mochi.yaml) and the following commands were used to predict fitness values for higher-order variants:

```
run_mochi.py --model_design model_design_CLN2.txt --output_directory CLN2 --predict
run_mochi.py --model_design model_design_CLO2.txt --output_directory CLO2 --predict
run_mochi.py --model_design model_design_LD.txt --output_directory LD --predict
```

Predictions of these three MoCHI runs were compiled for further R analysis using python scripts [Create_tableMutantCode.py](EDA_output/MoCHI/Create_tableMutantCode.py) and [CompileMoCHIpredictions_oneTable.py](EDA_output/MoCHI/CompileMoCHIpredictions_oneTable.py).

## Bioinformatic analysis of PacBio data

Files and input/output are collected in [PacBio_code](PacBio_code).

- [PacBio_barcode_processing.ipynb](PacBio_code/PacBio_barcode_processing.ipynb): scripts used for analysis, partially described below
- [ref_cbbM_onlyInsert.fa](PacBio_code/ref_cbbM_onlyInsert.fa) is the 1,675 bp long template spanning the promoter region, coding sequence of *cbbM* and a small region up- and downstream of the N20 barcode. This file was used as reference for mapping using minimap2 v2.28 (r1209). After quality-filtering, the resulting sam file was converted to fasta using the Jupyter notebook mentioned above.
- Cutadapt v3.5 was run was run twice with sequences exactly up- and downstream of the *cbbM* coding sequence and the N20 barcodes:

```
cutadapt -a ACTGAAACATCTTAATCATGCTAAGGAGGTTTTCTA...TCTAGAATCGCCGAAAGTAATTCAACTCCATTAA --discard-untrimmed -o CDS_trim.fasta aligned_seq.fasta > cutadapt_CDSreport.txt

cutadapt -a TCTAGAATCGCCGAAAGTAATTCAACTCCATTAA...TCTAGATGCTTACTAGTTACCGCGGCCAGGCAT --discard-untrimmed -o barcode_trim.fasta aligned_seq.fasta 2> cutadapt_barcode_report.txt
```

- The above mentioned jupyter notebook was used to assign barcodes to coding sequences and translate the latter to protein sequences. In the same step, the set of extracted barcodes was separated into three subsets according to the number of their occurrences. Furthermore, coding sequences translated into protein products shorter than 400 amino acids were discarded. 
- For barcodes occurring more than once, the assigned coding sequences were compared and the sequence occurring most frequently was used as consensus sequence. In case there were ties, these barcodes were omitted from future analysis.
- In a final step, all barcodes were combined into one file, amino acid exchanges compared to the CbbM(I) sequence were determined and the combination of barcodes and assigned mutations were formatted to be compatible with downstream applications.
- The nf-core-crispripipeline was run once with all detected barcodes and in a following steps, only barcodes which occurred in the *Synechocystis* library were used for follow-up analyses. Run Nextflow once with [unfiltered_barcodes.fasta](PacBio_code/output/unfiltered_barcodes.fasta), then select  all_counts.tsv (mapping of Illumina data from first large library cultivation) For further analyses, mutant variants which were not part of the originally designed library were only kept if they were represented by at least two barcodes

- [createTable_singleCombinatorial.py](PacBio_code/createTable_singleCombinatorial.py): script to create file with overview of mutant variant name and if saturational/combinatorial/...

## Illumina NGS bioinformatic analysis and exploratory data analysis

### Illumina sequencing

All files related to sequencing are collected in [Sequencing](Sequencing).

#### Sequencing on 24th May 2024

On a NextSeq2000, a P2 flow cell was used to obtain results for the first cultivation of the large Rubisco library. For demultiplexing, BCL Convert v2.4.0 was used. Downloaded using the command "bs download project --name large_Rubisco_firstCultiv" (compare [Basespace CLI documentation](https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview)). File [sort_files.sh](input/sort_files.sh) was used to combine all fastq files into one folder.

#### Sequencencing on 08th November 2024

Cultivation of large library with 5% CO2, 70% N2, 20% O2 in continuous light or light-dark cycles. A P2 flow cell was used on a NextSeq2000 system. BCL Convert v2.4.0 was used for demultiplexing.

### nf-core-crispripipeline analysis 

#### Input

Pipeline commit used for analysis: [d9e9b1812b2e093a49f290cd0b04fc407ff1d5a2](https://github.com/MPUSP/nf-core-crispriscreen/commit/d9e9b1812b2e093a49f290cd0b04fc407ff1d5a2) from 31st Oct 2024. Nextflow version 24.04.4, build 5917 (compare Nextflow report in Results folder)

```
conda activate env_nf
nextflow run ../../pipelines/nf-core-crispriscreen/ -c custom.config -profile singularity --input "input/samplesheet_1stCultiv_2ndCultiv.csv" --fasta "input/selected_barcodes_Nextflow.fa" --outdir "nf-core-crispripipeline_Results" --five_prime_adapter GTCTAGAatcgccgaaagtaattcaactccattaa...TCTAGATGCTTACTAGTTACCGCGGCCA --error_rate 0.2 --filter_mapq=1 --run_mageck false --gene_fitness false
```

Results are collected in directory [nf-core-crispripipeline_Results](nf-core-crispripipeline_Results).

#### Weighting

CRISPRi correlation and efficiency as performed by the pipeline only makes limited sense, so different weighting strategies were tried. Compare description in [compare_weightingStrategies.Rmd](code/compare_weightingStrategies.Rmd) and [compare_weightingStrategies.html](code/compare_weightingStrategies.html). Code and output are collected in [WeightingStrategies_fitnessCalc](WeightingStrategies_fitnessCalc). Final decision on [calculate_fitness_Weighted_lfcSE_Pearson_sameWeight.R](WeightingStrategies_fitnessCalc/calculate_fitness_Weighted_lfcSE_Pearson_sameWeight.R). 

Followed by analysis using code based on [fitness summary file of nf-core-crispriscreen pipeline](https://github.com/MPUSP/nf-core-crispriscreen/blob/master/bin/fitness_summary.Rmd).

### Further analyses

[This .Rmd file](code/EDA.Rmd) collects relevant analyses underlying Figures in the manuscript. Epistasis was determined using script [epistasis.py](code/epistasis.py).

## Citations & GitHub repositories

- [EVmutation repository](https://github.com/debbiemarkslab/EVmutation) and corresponding manuscripts [https://doi.org/10.1038/nbt.3769](https://doi.org/10.1038/nbt.3769) as well as [https://doi.org/10.1093/bioinformatics/bty862](https://doi.org/10.1093/bioinformatics/bty862)
- [https://doi.org/10.1038/s41592-018-0138-4](https://doi.org/10.1038/s41592-018-0138-4) and the [repository](https://github.com/debbiemarkslab/DeepSequence)
- [https://doi.org/10.1016/j.cels.2021.07.008](https://doi.org/10.1016/j.cels.2021.07.008) and the [corresponding repository](https://github.com/fhalab/MLDE)
- [https://doi.org/10.1038/s41587-021-01146-5](https://doi.org/10.1038/s41587-021-01146-5) and the [corresponding repository](https://github.com/chloechsu/combining-evolutionary-and-assay-labelled-data)
- [doi.org/10.1021/acssynbio.5c00065](https://doi.org/10.1021/acssynbio.5c00065) and [EVmut_inSilico repository](https://github.com/ute-hoffmann/EVmut_inSilico)
- [MSA Transformer](https://proceedings.mlr.press/v139/rao21a.html?utm_source=miragenews&utm_medium=miragenews&utm_campaign=news)
- [ProteinNPT repository](https://github.com/OATML-Markslab/ProteinNPT), compare the [manuscript](https://openreview.net/forum?id=AwzbQVuDBk)
- [MoCHI paper](https://doi.org/10.1186/s13059-024-03444-y) and [repository](https://github.com/lehner-lab/MoCHI)
