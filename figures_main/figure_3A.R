#!/usr/bin/env Rscript

#USE: Rscript executable_figure_3A.R [path/to/te_meta_summary_stats_by_pheno/] [path/to/output_directory]
#te_meta_summary_stats_by_pheno stands for summary statistics of trans-ethnic meta analyses for each phenotype

args = commandArgs(trailingOnly=TRUE)

require(dplyr)
require(data.table)
require(ggplot2)
#require(R.utils)

phenos <- c('COVIDHOSPITALIZED_2COVIDHOSP_1COVIDNOHOSP','COVIDHOSPITALIZEDPOP_2COVIDHOSP_1POPULATION','COVID_SEVERITY_CONTINUOUS_PC',
            'COVIDSWABTEST_2COVIDPOS_1COVIDNEG','COVIDSWABTESTPOP_2COVIDPOS_1POPULATION','COVIDWITHPOSITIVEHOUSEMATE_2COVIDPOS_1COVIDNEG',
            'COVIDWITHPOSITIVEHOUSEMATEPOP_2POPULATION_1COVIDNEGWITHEXP','MILDASYMPTOMATICCOVID_2SYMPTOMATIC_1MILDASYMPTOMATIC')
te_meta_path <- args[1]
output_path <- args[2]
freq_breaks <- c(0.01,0.05,0.1,0.2,0.3,0.4,0.5)

window_size <- 500000

filter_top_signal_in_window <- function(dat_chri, pval_thr){
  top_signals <- NULL
  dat_chri <-  dat_chri %>% arrange(BP) %>% mutate(window = cut_width(BP, window_size, boundary = 1)) %>% group_by(window) %>% arrange(P_value) %>% slice(1) #split the chromosome into windows of specified size, then take the top SNP in each window
  dat_chri <- data.frame(dat_chri)
  dat_chri <- dat_chri %>% arrange(BP) 
  dat_chri <- dat_chri %>% filter(P_value<=pval_thr) %>% mutate(prior_diff = BP - lag(BP), next_diff= lead(BP) - BP) #remove SNPs below the P-value threshold, then see if any adjacent windows are too close (difference between top SNPs < specified window size)
  rmove <- dat_chri %>% filter(P_value>lag(P_value) | P_value>lead(P_value)) %>% filter(prior_diff<window_size | next_diff<window_size) #get a list of SNPs that are too close
  dat_chri <- dat_chri %>% filter(!MarkerName %in% rmove$MarkerName) #remove the SNPs that are too close; retain the SNP with the more significant Pval
  excluded_vars <- c("prior_diff", "next_diff", "window")#clean up the data frame
  dat_chri  <- dat_chri  %>% select(-excluded_vars)
  top_signals <- rbind(top_signals, dat_chri)
  return(top_signals)
}

for(pval_thr in c(1e-5)){
  print(pval_thr)
  protective_by_maf_pheno <- lapply(phenos, function(pp){
    print(pp)
    Sys.time()
    df_te_sum_stat <- do.call('rbind',lapply(c(1:22,'X'), function(chri){
      
      #Filter associations present in all 3 cohorts
      df <- fread(paste0(te_meta_path,pp,'/te_meta/METAL-IVW_EUR_LAT_AA_chr',chri,'_TEMeta1.tbl')) %>%
        #select(c('MarkerName','Freq1','Effect','Direction','P-value')) %>%
        mutate(BP = as.numeric(sapply(strsplit(MarkerName,':'),'[[',2))) %>%
        filter(!grepl('?',Direction,fixed = T)) %>%
        rename('P_value'='P-value') %>%
        mutate('Phenotype'=pp)
      
      #Get top SNP in 500kb window passing p-value threshold
      df_top_in_window <- filter_top_signal_in_window(dat_chri = df, pval_thr = pval_thr)
      #Get effect size for minor allele
      df_top_in_window[df_top_in_window$Freq1>0.5,]$Effect <- -1*df_top_in_window[df_top_in_window$Freq1>0.5,]$Effect
      #Set Freq1 to MAF
      df_top_in_window[df_top_in_window$Freq1>0.5,]$Freq1 = 1-df_top_in_window[df_top_in_window$Freq1>0.5,]$Freq1
      df_top_in_window
    }))
  })
}


#combine everything together
snps_for_plot <- do.call(rbind, protective_by_maf_pheno)


#fix the names
snps_for_plot$Phenotype[snps_for_plot$Phenotype=="COVIDHOSPITALIZED_2COVIDHOSP_1COVIDNOHOSP"] <- "Hospitalized/Not_Hospitalized"
snps_for_plot$Phenotype[snps_for_plot$Phenotype=="COVIDHOSPITALIZEDPOP_2COVIDHOSP_1POPULATION"] <- "Hospitalized/Unscreened"
snps_for_plot$Phenotype[snps_for_plot$Phenotype=="COVIDSWABTEST_2COVIDPOS_1COVIDNEG"] <- "Positive/Negative"
snps_for_plot$Phenotype[snps_for_plot$Phenotype=="COVIDSWABTESTPOP_2COVIDPOS_1POPULATION"] <- "Positive/Unscreened"
snps_for_plot$Phenotype[snps_for_plot$Phenotype=="COVIDWITHPOSITIVEHOUSEMATE_2COVIDPOS_1COVIDNEG"] <- "Exposed_Positive/Exposed_Negative"
snps_for_plot$Phenotype[snps_for_plot$Phenotype=="COVIDWITHPOSITIVEHOUSEMATEPOP_2POPULATION_1COVIDNEGWITHEXP"] <- "Unscreened/Exposed_Negative"
snps_for_plot$Phenotype[snps_for_plot$Phenotype=="MILDASYMPTOMATICCOVID_2SYMPTOMATIC_1MILDASYMPTOMATIC"] <- "Symptomatic/Paucisymptomatic"
snps_for_plot$Phenotype[snps_for_plot$Phenotype=="COVID_SEVERITY_CONTINUOUS_PC"] <- "Continuous_Severity_Score"

###write out the file
write.table(snps_for_plot, paste0(output_path, "TopSNPs_WithDirections.csv"), sep=",", quote=F, row.names = F)

####Make a plot of the % protective suggestive minor alleles by Phenotype.
#look at direction of effect
prot_sum <- snps_for_plot %>% mutate(direction_of_effect=ifelse(sign(Effect)==1, "Risk", "Protective")) %>%
  group_by(Phenotype, direction_of_effect) %>%
  summarise(n=n()) %>%
  spread(direction_of_effect, n, fill =0) %>%
  mutate(total=Risk+Protective, percent_protective=Protective/total) %>%
  data.frame()

prot_sum$Phenotype <- factor(prot_sum$Phenotype, ordered = TRUE, 
                             levels = c("Hospitalized/Not_Hospitalized", "Hospitalized/Unscreened", "Positive/Negative", "Positive/Unscreened",
                                        "Continuous_Severity_Score", "Exposed_Positive/Exposed_Negative", "Unscreened/Exposed_Negative", "Symptomatic/Paucisymptomatic"))



overall_mean=sum(prot_sum$Protective)/sum(prot_sum$total)

p <- ggplot(prot_sum, aes(x=Phenotype, y=percent_protective, size = total)) +
  theme(strip.text.x = element_text(size = 10)) +
  xlab("") +
  ylab("% Protective Suggestive Lead SNPs")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, face = "italic"))+
  #theme(legend.position = "none") +
  scale_y_continuous(labels = function(percent_protective) paste0(percent_protective*100, "%")) +
  geom_point(aes(fill = percent_protective), shape=21, color = "black", stroke = 1) +
  scale_fill_gradient2(low = "#CA3542", mid = "white", high = "#27647B", midpoint = 0.5) +
  guides(fill=FALSE, size=guide_legend(title=str_wrap("Total # Associated Loci", 10)))+
  geom_hline(yintercept = overall_mean, linetype='dotted')

ggsave(paste0(output_path, "figure3A_percent_protective.eps"), plot=p, device = "eps")


