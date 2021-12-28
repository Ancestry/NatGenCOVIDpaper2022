#!/usr/bin/env Rscript

#USE: Rscript executable_figure_3B.R [path/to/te_meta_summary_stats_by_pheno/] [path/to/output_directory]
#te_meta_summary_stats_by_pheno stands for summary statistics of trans-ethnic meta analyses for each phenotype

args = commandArgs(trailingOnly=TRUE)

library(stringr)
library(dplyr)
library(ggplot2)

te_meta_path <- args[1]
output_path <- args[2]

setwd(paste0(te_meta_path, "protective_effects_enrichment/"))

files <- list.files(path = paste0(te_meta_path, "protective_effects_enrichment/"))

out <- NULL
for(file in files){
  ob <- readRDS(file)
  pval <- str_match(file, "pvalueThreshold_\\s*(.*?)\\s*.rds")[,2]
  
  int_out <- NULL
  for(i in 1:8){
    #get the right dataframes
    pheno_df1 <- data.frame(ob[[1]][i])
    pheno_df2 <- data.frame(ob[[2]][i])
    colnames(pheno_df1) <- c("maf_bin", "protective_YES_count", "protective_NO_count", "protective_YES_proportion", "protective_NO_proportion")
    colnames(pheno_df2) <- c("maf_bin", "protective_YES_count", "protective_NO_count", "protective_YES_proportion", "protective_NO_proportion")
    
    print(pheno_df1)
    
    maf_bins <- pheno_df1$maf_bin
    pheno_name <- names(ob[[1]][i])
    
    snps_array <- array(c(pheno_df1$protective_YES_count[1],
                          pheno_df2$protective_YES_count[1],
                          pheno_df1$protective_NO_count[1],
                          pheno_df2$protective_NO_count[1],
                          pheno_df1$protective_YES_count[2],
                          pheno_df2$protective_YES_count[2],
                          pheno_df1$protective_NO_count[2],
                          pheno_df2$protective_NO_count[2],
                          pheno_df1$protective_YES_count[3],
                          pheno_df2$protective_YES_count[3],
                          pheno_df1$protective_NO_count[3],
                          pheno_df2$protective_NO_count[3],
                          pheno_df1$protective_YES_count[4],
                          pheno_df2$protective_YES_count[4],
                          pheno_df1$protective_NO_count[4],
                          pheno_df2$protective_NO_count[4],
                          pheno_df1$protective_YES_count[5],
                          pheno_df2$protective_YES_count[5],
                          pheno_df1$protective_NO_count[5],
                          pheno_df2$protective_NO_count[5],
                          pheno_df1$protective_YES_count[6],
                          pheno_df2$protective_YES_count[6],
                          pheno_df1$protective_NO_count[6],
                          pheno_df2$protective_NO_count[6]),
                        dim=c(2,2,6),
                        dimnames = list(
                          Phenotype = c(pheno_name, "Else"),
                          Protective = c("Yes", "No"),
                          MAF = maf_bins))
    
    cmh <- mantelhaen.test(snps_array)
    cmh_out <- c(cmh$estimate, cmh$p.value, pheno_name)
    int_out <- rbind(int_out, cmh_out)
  }
  
  add_pval <- cbind(int_out , pval)
  out <- rbind(out, add_pval)
}

#manipulate the results
res <- data.frame(out, stringsAsFactors = F)
colnames(res) <- c("or", "p", "pheno", "pval")
rownames(res) <- NULL
res$or <- as.numeric(res$or)
res$p <- as.numeric(res$p)
res$pval <- as.numeric(res$pval)


res <- res %>% mutate(Enrichment=ifelse(or<1, -1/or, or)) %>%
  mutate(significant=ifelse(p<1E-100, "***", 
                            ifelse(p<1E-10, "**",
                                   ifelse(p<0.0015625, "*", ""))))

#fix the names
res$pheno[res$pheno=="COVIDHOSPITALIZED_2COVIDHOSP_1COVIDNOHOSP"] <- "Hospitalized/Not_Hospitalized"
res$pheno[res$pheno=="COVIDHOSPITALIZEDPOP_2COVIDHOSP_1POPULATION"] <- "Hospitalized/Unscreened"
res$pheno[res$pheno=="COVIDSWABTEST_2COVIDPOS_1COVIDNEG"] <- "Positive/Negative"
res$pheno[res$pheno=="COVIDSWABTESTPOP_2COVIDPOS_1POPULATION"] <- "Positive/Unscreened"
res$pheno[res$pheno=="COVIDWITHPOSITIVEHOUSEMATE_2COVIDPOS_1COVIDNEG"] <- "Exposed_Positive/Exposed_Negative"
res$pheno[res$pheno=="COVIDWITHPOSITIVEHOUSEMATEPOP_2POPULATION_1COVIDNEGWITHEXP"] <- "Unscreened/Exposed_Negative"
res$pheno[res$pheno=="MILDASYMPTOMATICCOVID_2SYMPTOMATIC_1MILDASYMPTOMATIC"] <- "Symptomatic/Paucisymptomatic"
res$pheno[res$pheno=="COVID_SEVERITY_CONTINUOUS_PC"] <- "Continuous_Severity_Score"




res$Phenotype <- factor(res$pheno, ordered = TRUE, 
                        levels = rev(c("Hospitalized/Not_Hospitalized", "Hospitalized/Unscreened", "Positive/Negative", "Positive/Unscreened",
                                       "Continuous_Severity_Score", "Exposed_Positive/Exposed_Negative", "Unscreened/Exposed_Negative", "Symptomatic/Paucisymptomatic")))

xax <- expression(paste("Trans-ancestry ", italic("P"), "-value Threshold"))

p <- ggplot(res, aes(factor(pval), Phenotype, fill= Enrichment)) + 
  geom_tile() +
  scale_fill_gradient2(low = "#CA3542", mid = "white", high = "#27647B")+
  geom_text(aes(label=significant))+
  theme_classic() +
  labs(title = "", x = xax, y = "", fill = "Enrichment OR") +
  theme(axis.text.y = element_text(face = "italic"))

ggsave(output_path, "figure3B_cmh_enrichment_figure.eps", p, device="eps")
