**PREPRINT**

This repository and website documents all analyses, summary, tables and figures associated with the following PREPRINT: Integration of GWAS and TWAS to elucidate the genetic architecture of natural variation for leaf cuticular conductance in maize (https://www.biorxiv.org/content/10.1101/2021.10.26.465975v1).

**Abstract**

The cuticle, a hydrophobic layer of cutin and waxes synthesized by plant epidermal cells, is the major barrier to water loss when stomata are closed. Dissecting the genetic architecture of natural variation for maize leaf cuticular conductance (gc) is important for identifying genes relevant to improving crop productivity in drought-prone environments. To this end, we performed an integrated genome- and transcriptome-wide association study (GWAS/TWAS) to identify candidate genes putatively regulating variation in leaf gc. Of the 22 plausible candidate genes identified, five were predicted to be involved in cuticle precursor biosynthesis and export, two in cell wall modification, nine in intracellular membrane trafficking, and seven in the regulation of cuticle development. A gene encoding an INCREASED SALT TOLERANCE1-LIKE1 (ISTL1) protein putatively involved in intracellular protein and membrane trafficking was identified in GWAS and TWAS as the strongest candidate causal gene. A set of maize nested near-isogenic lines that harbor the ISTL1 genomic region from eight donor parents were evaluated for gc, confirming the association between gc and ISTL1 in a haplotype-based association analysis. The findings of this study provide novel insights into the role of regulatory variants in the development of the maize leaf cuticle, and will ultimately assist breeders to develop drought-tolerant maize for target environments.

**Data availability**

Much of the supporting data and output from the analyses documented here are too large for GitHub.

**Scripts 1 to 3 are related to the expression data processing, BLUP calculation and outlier removal**

1.format_GE_filtering.R performs filtering for libary quality, rlog transformation on read counts files and remove genes that has 0 values for more than 50% of the population

2.BLUP_combGeneExp.R calculates BLUPs for each gene to combine gene expression from two environments, and PEER residuals are calculated for TWAS

3.OutlierRM_Transcription.R identifies outliers for PEER residuals and changes them to NA's

**Script 4 performs TWAS analysis**

4.TWAS_rrBLUP.R performs TWAS using the rrBLUP R package and selected candidate genes as descried in the manuscript

**Scripts 7.0 to 7.2 are related to GWAS, peak identification and candidate gene selection**

7.0.GAPIT_GWAS.R performs GWAS for gc and generates Manhattan plot for GWAS results

7.1.define_GWAS_peaks.R identifies peaks using top 0.01% association SNPs in GWAS

7.2.candidate_selection_GWAS.R identifies candidate genes in +/-200kb windows surrounding GWAS peaks

**Script 9 generates the gene-SNP pair for the Fisher's combined test (FCT)**

9.SNP_gene_pair.R generates the gene-SNP pair for the Fisher's combined test (FCT)

**Script 10 is related to the FCT**

10.SNP_assign_Fisher.R performs the FCT and selected the top 1% associated genes for gc

**Scripts 11 and 12 are related to random forest regression for cuticular conductance using cuticular waxes and cuticle features**

11.process_wax_raw_data.R removes outliers for cuticular waxes and cuticle features and calculates plot means for these features

12.RF_ce_pred_cforest.R performs random forest regression for cuticular conductance using cuticular waxes and cuticle features

**Scripts 13.1 to 13.5 are related to haplotype analysis for the ISTL1 genomic region**

13.1_format_vcf_to_haploview_format.py formats genotype in vcf files to phased format that can be imported into Haploview

13.2_split_by_chr.py splits genotype files based on markers' chromosome location

13.3_run_haploview.sh is a bash file that run Haploview in the linux system

13.4_build_block_dict.py generates haplotype files using Haploview blocks and genotype files

13.5.haplotype_gwas_asreml.R performs haplotype-based association analysis for the ISTL1 region in the diversity panel

**Script 14 is related to haplotype-based analysis for the ISTL1 region among nNILs**

14.hapotype_effect_nNILs_combEnv.R estimates halotype effects for top associated haploblocks identifed in the diversity panel for nNILs
