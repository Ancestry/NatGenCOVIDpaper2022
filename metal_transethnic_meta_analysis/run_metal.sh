#!/usr/bin/env bash

###define genotype source (output of prepared by prepare_gwas_output_for_trans_meta_analysis.R), phenotype file path prefix
output_directory=/path/to/gwas_summary_stats/
pheno_file_path=/path/to/phenotypes

#list of phenotypes for which a separate trans-ethnic meta analysis would be run
PHENOS=(
COVIDHOSPITALIZED_2COVIDHOSP_1COVIDNOHOSP
COVIDHOSPITALIZEDPOP_2COVIDHOSP_1POPULATION
COVIDWITHPOSITIVEHOUSEMATE_2COVIDPOS_1COVIDNEG
COVIDWITHPOSITIVEHOUSEMATEPOP_2POPULATION_1COVIDNEGWITHEXP
MILDASYMPTOMATICCOVID_2SYMPTOMATIC_1MILDASYMPTOMATIC
COVIDSWABTEST_2COVIDPOS_1COVIDNEG
COVIDSWABTESTPOP_2COVIDPOS_1POPULATION
)

#############Run METAL on males and females separately for the EUR population##############
for pheno in $PHENOS
do
  for ethnicity in eur
  do
    for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
    do
      f_metal_in=$output_directory'/'$pheno'/gwas_by_eth/'$ethnicity'_chr'$chr'_females_input_sumstats_forMETAL.txt'
      m_metal_in=$output_directory'/'$pheno'/gwas_by_eth/'$ethnicity'_chr'$chr'_males_input_sumstats_forMETAL.txt'
      metal_out=$output_directory'/'$pheno'/te_meta/METAL-IVW_'$ethnicity'_chr'$chr'_sexMeta'
    
      cat > $metal_out'.config' <<EOL
AVERAGEFREQ ON
MINMAXFREQ ON
MARKER   MarkerName
ALLELE   Allele1 Allele2
FREQ     Freq1
EFFECT   Effect
PVAL     P-value
SCHEME STDERR
STDERR StdErr
PROCESS ${f_metal_in}
PROCESS ${m_metal_in}
OUTFILE ${metal_out} .tbl
ANALYZE HETEROGENEITY
... 
EOL

    /path/to/software/metal $metal_out'.config' > $metal_out'.METAL.log'
    done
  done
done
##############################Do Trans-Ethnic Meta across EUR sex meta analysis, LAT, and AA cohort-specific GWAS########################################

for pheno in $PHENOS
do
  for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
    do
  
    metal_out_transEth=$output_directory'/'$pheno'/te_meta/METAL-IVW_EUR_LAT_AA_chr'$chr'_TEMeta'
  
    eur_metal_in=$output_directory'/'$pheno'/te_meta/METAL-IVW_eur_chr'$chr'_sexMeta1.tbl'
    lat_metal_in=$output_directory'/'$pheno'/gwas_by_eth/lat_chr'$chr'_input_sumstats_forMETAL.txt'
    aa_metal_in=$output_directory'/'$pheno'/gwas_by_eth/aa_chr'$chr'_input_sumstats_forMETAL.txt'
    
    cat > $metal_out_transEth'.config' <<EOL
AVERAGEFREQ ON
MINMAXFREQ ON
MARKER   MarkerName
ALLELE   Allele1 Allele2
FREQ     Freq1
EFFECT   Effect
PVAL     P-value
SCHEME STDERR
STDERR StdErr
PROCESS ${eur_metal_in}
PROCESS ${lat_metal_in}
PROCESS ${aa_metal_in}
#PROCESS ${eas_metal_in}
OUTFILE ${metal_out_transEth} .tbl
ANALYZE HETEROGENEITY
... 
EOL
  /path/to/software/metal $metal_out_transEth'.config' > $metal_out_transEth'.METAL.log'

  done
done