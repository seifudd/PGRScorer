#!/bin/bash
#$ -cwd

###############################################################################################################
#
# Script to calculate polygenic/genomic risk score (using PLINK 2.0)
# Download PLINK 2.0 if not available on a computing cluster already = https://www.cog-genomics.org/plink/2.0/
#
# 1. Input: File (space delimited): "SNPid (matching with genotype data)" "Effect Allele" "Beta/Odds Ratio"
# Description: This file was provided from the UK Biobank contains the following information:
#              "chr" "position" "rsid" "allele1" "allele2" "effect_allele" "beta"
# Description: This file is the combined beta (across 3 data sets) from the metaGRS CAD score developed 
#              using the UK Biobank data [PMID: 30309464]
# 
# 2. Input: Files BED, BIM and FAM file [PLINK binary format]: Genotype data
# Description: https://www.cog-genomics.org/plink/2.0/input#bed
#
# 3. Ouput: Score file:
# Description: https://www.cog-genomics.org/plink/2.0/formats#sscore
#
###############################################################################################################

# Working directory, change accordingly
workingdir="/data/NHLBI_BCB/Arai_Lab/00-metaGRS"

# IF using hg19 coordinates for SNPid, metaGRS Beta/Odds Ratio file name
beta_file_name="metaGRS_hg19_20180205_v2.txt"

# IF using hg38 coordinates for SNPid, metaGRS Beta/Odds Ratio file name
# beta_file_name="metaGRS_hg38_06102020.txt"

# BED, BIM and FAM file
# for e.g. if calculating the score in AGES and the PLINK binary genotype files are labelled as: AGES.bed, AGES.bim and AGES.fam
genotype_file_name="AGES"

# Calculate polygenic/genomic score in PLINK 2.0
# For a full list of options = https://www.cog-genomics.org/plink/2.0/score

plink --bfile "$workingdir/$genotype_file_name" \
      --score "$workingdir/$beta_file_name" \
      3 \
      6 \
      7 \
      header-read \
      no-mean-imputation \
      ignore-dup-ids \
     --out "$workingdir/${genotype_file_name}.grs.score"

# Output will be named for e.g. AGES.grs.score
# Description of output = https://www.cog-genomics.org/plink/2.0/formats#sscore

