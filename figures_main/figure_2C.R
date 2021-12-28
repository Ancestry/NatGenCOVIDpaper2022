#!/usr/bin/env Rscript

# get effect sizes for replication of SNPs in the meta-analysis results. Also create forest plot of results.

#USE: Rscript figure_2C.R [/path/to/HGI_REGN_SNPs_to_replicate.txt] [/path/to/gwas_te_meta_summary_stats_by_pheno/] [/path/to/output_directory]
#gwas_te_meta_summary_stats_by_pheno stands for directory containing summary statistics of trans-ethnic meta analyses as well as individual cohort GWAS for each phenotype

library(tidyverse)
library(ggforestplot)

args = commandArgs(trailingOnly=TRUE)

hgi_regn_snps_list_path <- args[1]
gwas_te_meta_path <- args[2]
output_path <- args[3]

# Use zgrep on the command line to pull the snps out, rather than reading everything in:
cmd <- paste0("zgrep -F -f ", hgi_regn_snps_list_path," ", gwas_te_meta_path, "/te_meta/*Trans* > snp_list_sumstats.txt")
system(cmd)
# (then use emacs query-replace-regexp to edit away extra text)
# and to fetch the original sumstats that contain MAFs:
cmd <- paste0("grep -F -f ", hgi_regn_snps_list_path, " ", gwas_te_meta_path, "/te_meta/*.allChr.*tbl > snp_list_AFs.txt")
system(cmd)

sumstats <- read.csv("snp_list_sumstats.txt", header = T, stringsAsFactors = F)
freqstats <- read.csv("snp_list_AFs.txt", header = T, stringsAsFactors = F, sep="\t")

# check whether each SNP has the same Allele1
for(snp in unique(sumstats$MarkerName)){
  print(sumstats[sumstats$MarkerName==snp,c("Phenotype", "MarkerName", "Allele1", "Effect")])
}
# Yes. That's good!

# now look whether each snp maintains the same minor allele freq
for(snp in unique(freqstats$MarkerName)){
  print(freqstats[freqstats$MarkerName==snp,c("Phenotype","MarkerName", "Allele1", "Freq1")])
  print(freqstats[freqstats$MarkerName==snp,"Freq1"] > 0.5)
}
# every snp has the same minor allele across phenotypes

# make a little table of which SNPs have Allele1's that are minor alleles
alleles <- freqstats[freqstats$Phenotype=="COVIDHOSPITALIZED_2COVIDHOSP_1COVIDNOHOSP", c("MarkerName", "Allele1", "Freq1")]
alleles$Allele1_is_minor <- ifelse(alleles$Freq1 < 0.5, "TRUE", "FALSE")

# join up the minor allele column with sumstats from original file read
sumstats_alleles <- merge(sumstats, alleles[,-3], by=c("MarkerName", "Allele1"))

# and now make an effect column that has the direction flipped if the effect allele is not the minor allele
# (so we are always thinking of the effect in terms of the minor allele)
sumstats_alleles$Effect_of_minor <- ifelse(sumstats_alleles$Allele1_is_minor, sumstats_alleles$Effect, - sumstats_alleles$Effect)

# rename markers
sumstats_alleles$MarkerName[sumstats_alleles$MarkerName=="9:136149229:T:C"] <- "ABO (rs505922)"
sumstats_alleles$MarkerName[sumstats_alleles$MarkerName=="6:31121426:G:A"] <- "CCHCR1 (rs143334143)"
sumstats_alleles$MarkerName[sumstats_alleles$MarkerName=="19:4723670:C:A"] <- "DPP9 (rs2277732)"
sumstats_alleles$MarkerName[sumstats_alleles$MarkerName=="21:34615210:T:C"] <- "IFNAR2 (rs13050728)"
sumstats_alleles$MarkerName[sumstats_alleles$MarkerName=="12:103014757:C:A"] <- "IGF1 (rs10860891)"
sumstats_alleles$MarkerName[sumstats_alleles$MarkerName=="17:47868921:G:A"] <- "KAT7 (rs9903642)"
sumstats_alleles$MarkerName[sumstats_alleles$MarkerName=="3:45889921:A:T"] <- "SLC6A20/LZTFL1 (rs35081325)"
sumstats_alleles$MarkerName[sumstats_alleles$MarkerName=="12:113392182:A:G"] <- "OAS3 (rs7310667)"
sumstats_alleles$MarkerName[sumstats_alleles$MarkerName=="3:45835417:G:A"] <- "SLC6A20/LZTFL1 (rs73062389)"
sumstats_alleles$MarkerName[sumstats_alleles$MarkerName=="3:45838300:G:A"] <- "SLC6A20/LZTFL1 (rs2531743)"
sumstats_alleles$MarkerName[sumstats_alleles$MarkerName=="7:54647894:A:C"] <- "STM2A (rs622568)"
sumstats_alleles$MarkerName[sumstats_alleles$MarkerName=="21:42845642:C:A"] <- "TMPRSS2 (rs2298661)"

# Rename phenotypes, and fix "hosptialized"
sumstats_alleles$Phenotype[sumstats_alleles$Phenotype=="Positive_Negative"] <- "Positive/Negative"
sumstats_alleles$Phenotype[sumstats_alleles$Phenotype=="Positive_Unscreened"] <- "Positive/Unscreened"
sumstats_alleles$Phenotype[sumstats_alleles$Phenotype=="ExposedPositive_ExposedNegative"] <- "Exposed_Positive/Exposed_Negative"
sumstats_alleles$Phenotype[sumstats_alleles$Phenotype=="Unscreened_ExposedNegative"] <- "Unscreened/Exposed_Negative"
sumstats_alleles$Phenotype[sumstats_alleles$Phenotype=="Hospitalized_NotHosptialized"] <- "Hospitalized/Not_Hospitalized"
sumstats_alleles$Phenotype[sumstats_alleles$Phenotype=="Hospitalized_Unscreened"] <- "Hospitalized/Unscreened"
sumstats_alleles$Phenotype[sumstats_alleles$Phenotype=="Symptomatic_Paucisymptomatic"] <- "Symptomatic/Paucisymptomatic"
sumstats_alleles$Phenotype[sumstats_alleles$Phenotype=="Continuous_Severity_Score"] <- "Continuous_Severity_Score"

# write the replication effect sizes out to a table for the paper supplement      
write.table(sumstats_alleles, file="SuppTab4_replication_effect_sizes.tsv", quote = F, sep = "\t", row.names = F)

# for plotting, drop out SNPs that don't have any results with p<0.05
# first identify them
exclude_snps <- c()
for(snp in unique(sumstats_alleles$MarkerName)){
  exclude <- sumstats_alleles[sumstats_alleles$MarkerName==snp,"Pvalue"] > 0.05
  if(all(exclude)){exclude_snps <- c(exclude_snps, snp)}
}

# now drop them from results we will plot
to_plot <- sumstats_alleles[ ! sumstats_alleles$MarkerName %in% exclude_snps, ]

####### plot ORs 

# group phenotypes by susceptibility and severity, and then by novelty and increasing sample size
to_plot$Phenotype <- factor(to_plot$Phenotype, levels=rev(c("Positive/Unscreened","Positive/Negative", 
                                                            "Unscreened/Exposed_Negative", "Exposed_Positive/Exposed_Negative", 
                                                            "Hospitalized/Unscreened", "Hospitalized/Not_Hospitalized",
                                                            "Symptomatic/Paucisymptomatic", "Continuous_Severity_Score")))


# pick same colors as for power plot (figure 2A)
colors <- rev(c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00"))

OR_plot <- ggforestplot::forestplot(
  df = to_plot,
  name = MarkerName,
  estimate = Effect_of_minor,
  se = StdErr,
  pvalue = Pvalue,
  colour = Phenotype,
  xlab = "Odds Ratio & 95% CI",
  logodds = TRUE
  # title="SNPs from previous studies replicated in at least one phenotype"
) + theme(axis.text.y = element_text(face="italic", size=12), legend.text=element_text(size=12)) +
  scale_color_manual(values=colors) + geom_stripes(odd = "#00000000", even = "#FFFFFF00")
ggsave(output_path, "figure2C_forest_plot.pdf", OR_plot, device="pdf")


### SAME THING WITH BETA SCALE
OR_plot <- ggforestplot::forestplot(
  df = to_plot,
  name = MarkerName,
  estimate = Effect_of_minor,
  se = StdErr,
  pvalue = Pvalue,
  colour = Phenotype,
  xlab = "Beta & 95% CI",
  logodds = FALSE
  # title="SNPs from previous studies replicated in at least one phenotype"
) + theme(axis.text.y = element_text(face="italic", size=12), legend.text=element_text(size=12)) +
  scale_color_manual(values=colors) 
ggsave(output_path, "figure2C_forest_plot_beta_scale.pdf", OR_plot, device="pdf")