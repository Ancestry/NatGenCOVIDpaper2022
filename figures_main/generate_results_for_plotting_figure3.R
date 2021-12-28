#!/usr/bin/env Rscript

#USE: Rscript executable_figure_3B.R [path/to/te_meta_summary_stats_by_pheno/] [path/to/output_directory]
#te_meta_summary_stats_by_pheno stands for summary statistics of trans-ethnic meta analyses for each phenotype

#Calculate summary statistics of protective effect associations for COVID 19 phenotypes
require(dplyr)
require(data.table)
require(ggplot2)
require(R.utils)

args = commandArgs(trailingOnly=TRUE)

phenos <- c('COVIDHOSPITALIZED_2COVIDHOSP_1COVIDNOHOSP','COVIDHOSPITALIZEDPOP_2COVIDHOSP_1POPULATION','COVID_SEVERITY_CONTINUOUS_PC',
            'COVIDSWABTEST_2COVIDPOS_1COVIDNEG','COVIDSWABTESTPOP_2COVIDPOS_1POPULATION','COVIDWITHPOSITIVEHOUSEMATE_2COVIDPOS_1COVIDNEG',
            'COVIDWITHPOSITIVEHOUSEMATEPOP_2POPULATION_1COVIDNEGWITHEXP','MILDASYMPTOMATICCOVID_2SYMPTOMATIC_1MILDASYMPTOMATIC')

te_meta_path <- args[1]
output_path <- args[2]

freq_breaks <- c(0.01,0.05,0.1,0.2,0.3,0.4,0.5)

window_size <- 500000

filter_top_signal_in_window <- function(dat_chri, pval_thr){
  #dat_chri is already filtered for a single chromosome
  top_signals <- NULL
  #dat_chri <- df
  #pval_thr <- 0.0001
  dat_chri <-  dat_chri %>% arrange(BP) %>% mutate(window = cut_width(BP, window_size, boundary = 1)) %>% group_by(window) %>% arrange(P_value) %>% slice(1) #split the chromosome into windows of specified size, then take the top SNP in each window
  dat_chri <- data.frame(dat_chri)
  dat_chri <- dat_chri %>% arrange(BP) 
  dat_chri <- dat_chri %>% filter(P_value<=pval_thr) %>% mutate(prior_diff = BP - lag(BP), next_diff= lead(BP) - BP) #remove SNPs below the P-value threshold, then see if any adjacent windows are too close (difference between top SNPs < specified window size)
  rmove <- dat_chri %>% filter(P_value>lag(P_value) | P_value>lead(P_value)) %>% filter(prior_diff<window_size | next_diff<window_size) #get a list of SNPs that are too close
  dat_chri <- dat_chri %>% filter(!MarkerName %in% rmove$MarkerName) #remove the SNPs that are too close; retain the SNP with the more significant Pval
  excluded_vars <- c("prior_diff", "next_diff", "window")#clean up the data frame
  dat_chri  <- dat_chri  %>% select(-excluded_vars)
  #colnames(dat_chri)[which(names(dat_chri) == "P")] <- P
  top_signals <- rbind(top_signals, dat_chri)
  #}
  return(top_signals)
}

for(pval_thr in c(1e-2,1e-3,1e-4,1e-5)){
  print(pval_thr)
  #calculate the stats for each pheno
  protective_by_maf_pheno <- lapply(phenos, function(pp){
    print(pp)
    Sys.time()
    df_te_sum_stat <- do.call('rbind',lapply(c(1:22,'X'), function(chri){
      #Filter associations present in all 3 cohorts
      df <- fread(paste0(te_meta_path,pp,'/te_meta/METAL-IVW_EUR_LAT_AA_chr',chri,'_TEMeta1.tbl')) %>%
        select(c('MarkerName','Freq1','Effect','Direction','P-value')) %>%
        mutate(BP = as.numeric(sapply(strsplit(MarkerName,':'),'[[',2))) %>%
        filter(!grepl('?',Direction,fixed = T)) %>%
        rename('P_value'='P-value')
      
      #Get top SNP in 500kb window passing p-value threshold
      df_top_in_window <- filter_top_signal_in_window(dat_chri = df, pval_thr = pval_thr)
      #Get effect size for minor allele
      df_top_in_window[df_top_in_window$Freq1>0.5,]$Effect <- -1*df_top_in_window[df_top_in_window$Freq1>0.5,]$Effect
      #Set Freq1 to MAF
      df_top_in_window[df_top_in_window$Freq1>0.5,]$Freq1 = 1-df_top_in_window[df_top_in_window$Freq1>0.5,]$Freq1
      df_top_in_window
    }))
    Sys.time()
    
    #classify MAF into frequency bins
    df_te_sum_stat$freq_level <- cut(df_te_sum_stat$Freq1, breaks = freq_breaks, include.lowest = T)
    #Identify protective associations
    df_te_sum_stat$Protective <- df_te_sum_stat$Effect < 0 
    #summarize counts of protective effect associations in each frequency bin 
    tab_protective_maf <- table(df_te_sum_stat$freq_level, df_te_sum_stat$Protective)
    df_protective_maf <- data.frame(maf_bin = rownames(tab_protective_maf),
                                    protective_YES_count = tab_protective_maf[,'TRUE'],
                                    protective_NO_count = tab_protective_maf[,'FALSE'], stringsAsFactors = F,
                                    row.names = NULL)
    df_protective_maf$protective_YES_proportion <- round(df_protective_maf$protective_YES_count/(df_protective_maf$protective_YES_count+df_protective_maf$protective_NO_count),4)
    df_protective_maf$protective_NO_proportion <- round(df_protective_maf$protective_NO_count/(df_protective_maf$protective_YES_count+df_protective_maf$protective_NO_count),4)
    df_protective_maf
  })
  names(protective_by_maf_pheno) <- phenos
  
  #calculate the stats for the rest of the 7 phenotypes for each pheno
  protective_by_maf_pheno_rest <- lapply(phenos, function(pp){
    #pp <- phenos[1]
    #print(pp)
    df <- protective_by_maf_pheno[[pp]]
    df$protective_YES_count <- df$protective_NO_count <- 0
    for(pp_r in setdiff(phenos,pp)){
      df$protective_YES_count <- df$protective_YES_count+protective_by_maf_pheno[[pp_r]]$protective_YES_count
      df$protective_NO_count <- df$protective_NO_count+protective_by_maf_pheno[[pp_r]]$protective_NO_count
    }
    df$protective_YES_proportion <- round(df$protective_YES_count/(df$protective_YES_count+df$protective_NO_count),4)
    df$protective_NO_proportion <- round(df$protective_NO_count/(df$protective_YES_count+df$protective_NO_count),4)
    df
  })
  names(protective_by_maf_pheno_rest) <- phenos
  
  #save stats to file
  saveRDS(list(protective_effects_stats_by_maf_pheno = protective_by_maf_pheno, protective_effects_stats_by_maf_pheno_rest = protective_by_maf_pheno_rest),
          file = paste0(output_path,'protective_effects_enrichment/protective_effects_enrichmentprotective_effects_stats_1_vs_rest_pvalueThreshold_',pval_thr,'.rds'))
}
