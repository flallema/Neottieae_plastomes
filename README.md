# Neottieae_plastomes

Scripts used for quality control of sequencing data associated with the following paper: Thirteen new plastid genomes from mixotrophic and autotrophic species provide insights into heterotrophy evolution in Neottieae orchids. Lallemand et al. 2019.

1. QualityControl_DB.cpp . This script is computed in C++. The goal of this script is control the quality of the sequence file in input. The amount of ATGCN and the number of redundant sequences are calculated and the boxplot graphes are then created.

2. SearchRedondances.cpp : This script is computed in C++. The goal of this script is the dentification of redondant sequences from rawdata by couple Read1-Read2 or in one single read.

3. Trimming.cpp : This script is computed in C++. The goal of this script is the trimming (as described in the readme) of the sequences.

4. param_trim_q30_l30_min30_N0_L151.txt . This file lists all the trimming parameters that we used. This file is used by the script Trimming.cpp

5. README_EPGV_DataTransfer_Illumina_Sequencing_v1.2.3_EN.pdf : This document explains the analyzing pipeline of our partner (EPGV team).
