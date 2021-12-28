#!/usr/bin/env Rscript

#USE: Rscript figure_2B.R [path/to/HGI_REGN_SNPs_to_replicate.csv] [path/to/gwas_te_meta_summary_stats_by_pheno/] [path/to/output_directory]
#gwas_te_meta_summary_stats_by_pheno stands for directory containing summary statistics of trans-ethnic meta analyses as well as individual cohort GWAS for each phenotype

library(data.table)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

gwas_te_meta_path <- args[2]
output_path <- args[3]

#read in the hgi and regeneron SNPs for replication analyses
regn <- read.csv(args[1])
chr_names <- as.character(unique(regn$CHR))

pheno_list <- c("COVIDHOSPITALIZED_2COVIDHOSP_1COVIDNOHOSP" ,
                "COVIDHOSPITALIZEDPOP_2COVIDHOSP_1POPULATION",
                "COVIDSWABTEST_2COVIDPOS_1COVIDNEG",
                "COVIDSWABTESTPOP_2COVIDPOS_1POPULATION",
                "MILDASYMPTOMATICCOVID_2SYMPTOMATIC_1MILDASYMPTOMATIC",
                "COVIDWITHPOSITIVEHOUSEMATE_2COVIDPOS_1COVIDNEG",
                "COVIDWITHPOSITIVEHOUSEMATEPOP_2POPULATION_1COVIDNEGWITHEXP",
                "COVID_SEVERITY_CONTINUOUS_PC")


replicate_regn_leadsnps <- NULL
for(pheno in pheno_list){
  phenotype_name=pheno
  chr_names <- as.character(unique(regn$CHR))
  
  chr_temp <-NULL
  for(i in chr_names){
    
    #pull in AA SNPs 
    aa_gwas_path <- paste0(gwas_te_meta_path, phenotype_name, "/gwas_by_eth/aa_chrXX_input_sumstats_forMETAL_updated.txt")
    
    aa_gwas <- NULL
    infile_name <- str_replace(aa_gwas_path, "XX", i)
    temp <- fread(infile_name, header=T, stringsAsFactors = F) %>% 
      filter(MarkerName %in% regn$SNP) %>%
      rename(P=`P-value`) %>% 
      mutate(Phenotype=phenotype_name, ethnicity="AA") %>% 
      select(CHR, MarkerName, BP, Allele1, Allele2, Effect, P, StdErr, Freq1, Phenotype, ethnicity)
    aa_gwas <- rbind(aa_gwas, temp)
    temp <- NULL
    
    
    #pull in LAT SNPs
    lat_gwas_path <- paste0(gwas_te_meta_path, phenotype_name, "/gwas_by_eth/lat_chrXX_input_sumstats_forMETAL_updated.txt")
    
    lat_gwas <- NULL
    infile_name <- str_replace(lat_gwas_path, "XX", i)
    temp <- fread(infile_name, header=T, stringsAsFactors = F) %>% 
      filter(MarkerName %in% regn$SNP) %>%
      rename(P=`P-value`) %>% 
      mutate(Phenotype=phenotype_name, ethnicity="LAT") %>% 
      select(CHR, MarkerName, BP, Allele1, Allele2, Effect, P, StdErr, Freq1, Phenotype, ethnicity)
    lat_gwas <- rbind(lat_gwas, temp)
    temp <- NULL
    
    
    #pull in EUR SNPs
    eur_gwas_path <- paste0(gwas_te_meta_path, phenotype_name, "/te_meta/METAL-IVW_eur_chrXX_sexMeta1.tbl")
    
    eur_gwas <- NULL
    infile_name <- str_replace(eur_gwas_path, "XX", i)
    temp <- fread(infile_name, header=T, stringsAsFactors = F) %>% 
      filter(MarkerName %in% regn$SNP) %>%
      rename(P=`P-value`) %>% 
      mutate(Phenotype=phenotype_name, ethnicity="EUR")
    temp <- temp %>% mutate(SNP=MarkerName) %>% 
      separate(SNP, c("CHR", "BP", "delete1", "delete2")) %>% 
      mutate(Allele1=toupper(Allele1), Allele2=toupper(Allele2)) %>%
      select(CHR, MarkerName, BP, Allele1, Allele2, Effect, P, StdErr, Freq1, Phenotype, ethnicity)
    eur_gwas <- rbind(eur_gwas, temp)
    temp <- NULL
    
    #join everything together
    sig <- rbind(aa_gwas, lat_gwas, eur_gwas)
    
    #re-join to get effect estimates, P-values, and allele frequencies for each study at all SNPs in sig
    aa_gwas_reduced <- aa_gwas %>% filter(MarkerName %in% sig$MarkerName) %>% 
      select(MarkerName, AA_A1=Allele1, AA_beta=Effect, AA_SEbeta=StdErr, AA_Freq=Freq1, AA_P=P)
    
    lat_gwas_reduced <- lat_gwas %>% filter(MarkerName %in% sig$MarkerName) %>% 
      select(MarkerName, LAT_A1=Allele1, LAT_beta=Effect, LAT_SEbeta=StdErr, LAT_Freq=Freq1, LAT_P=P)
    
    eur_gwas_reduced <- eur_gwas %>% filter(MarkerName %in% sig$MarkerName) %>% 
      select(MarkerName, EUR_A1=Allele1, EUR_beta=Effect, EUR_SEbeta=StdErr, EUR_Freq=Freq1, EUR_P=P)
    
    all_reduced <- full_join(aa_gwas_reduced, lat_gwas_reduced, by="MarkerName")
    all_reduced <- full_join(all_reduced, eur_gwas_reduced)
    
    #align everything to the same allele
    all_reduced <- all_reduced %>% mutate(align_allele=ifelse(!is.na(EUR_A1), EUR_A1,
                                                              ifelse(!is.na(LAT_A1), LAT_A1, AA_A1))) %>%
      mutate(AA_beta=ifelse(AA_A1==align_allele, AA_beta, -AA_beta)) %>%
      mutate(AA_Freq=ifelse(AA_A1==align_allele, AA_Freq, 1-AA_Freq)) %>%
      mutate(LAT_beta=ifelse(LAT_A1==align_allele, LAT_beta, -LAT_beta)) %>%
      mutate(LAT_Freq=ifelse(LAT_A1==align_allele, LAT_Freq, 1-LAT_Freq)) %>%
      mutate(EUR_beta=ifelse(EUR_A1==align_allele, EUR_beta, -EUR_beta)) %>%
      mutate(EUR_Freq=ifelse(EUR_A1==align_allele, EUR_Freq, 1-EUR_Freq))
    
    all_reduced <- all_reduced %>% select(MarkerName, align_allele, EUR_P, LAT_P, AA_P, EUR_beta, LAT_beta, AA_beta, EUR_SEbeta, LAT_SEbeta, AA_SEbeta, EUR_Freq, LAT_Freq, AA_Freq)
    
    #join with the sig data frame and match effect direction
    sig <- left_join(sig, all_reduced, by="MarkerName")
    sig <- sig %>% mutate(Effect=ifelse(Allele1==align_allele, Effect, -Effect)) %>%
      mutate(Freq1=ifelse(Allele1==align_allele, Freq1, 1-Freq1))
    
    
    #read in the meta analysis and join that P-value and I2
    meta_path <- paste0(gwas_te_meta_path, phenotype_name, "/te_meta/METAL-IVW_EUR_LAT_AA_chrXX_TEMeta1.tbl")
    
    meta <- NULL
    infile_name <- str_replace(meta_path, "XX", i)
    temp <- data.frame(fread(infile_name, header=T, stringsAsFactors = F)) %>% 
      select(MarkerName, meta_a1=Allele1, meta_P=P.value, meta_beta=Effect, meta_SEbeta=StdErr, meta_dir=Direction, meta_Isq=HetISq, meta_hetPval=HetPVal) %>%
      mutate(meta_a1=toupper(meta_a1))
    meta  <- rbind(meta, temp)
    
    sig <- left_join(sig, meta, by="MarkerName") %>%
      mutate(meta_beta=ifelse(meta_a1==align_allele, meta_beta, -meta_beta)) %>%
      mutate(Allele2=ifelse(Allele1==align_allele, Allele2, Allele1)) %>%
      select(-Allele1, -Allele2, -meta_a1) %>%
      rename(Allele1=align_allele)
    
    chr_temp <- rbind(chr_temp, sig)
    
  }
  replicate_regn_leadsnps<- rbind(replicate_regn_leadsnps, chr_temp)
}  


regn_rep <- replicate_regn_leadsnps %>% 
  group_by(MarkerName, Phenotype) %>%
  select(-Effect, -P, -StdErr, -ethnicity) %>%
  mutate(CHR=as.integer(CHR), BP=as.integer(BP)) %>%
  rename(SNP=MarkerName)%>%
  data.frame()


regn_rep <- left_join(regn_rep, regn, by=c("CHR", "BP", "SNP"))
regn_rep <- regn_rep %>% mutate(Reported_Beta=ifelse(Allele1==Effect_Allele, Reported_Beta, -Reported_Beta))
regn_rep <- regn_rep %>% filter(Use==1) %>%
  mutate(replicate=ifelse(meta_P<0.05 & sign(meta_beta)==sign(Reported_Beta), 1, 0)) %>%
  mutate(gene_snp=paste0(Gene,"-", SNP)) %>%
  group_by(gene_snp, Phenotype) %>%
  slice(1) %>%
  mutate(CHR=as.integer(CHR), BP=as.integer(BP)) %>%
  data.frame()


#make a replication matrix of locus x phenotype
rep_mat <- regn_rep %>%
  mutate(logP=ifelse(replicate==1, -log10(meta_P), log10(1))) %>%
  select(gene_snp, Phenotype, logP) %>%
  spread(gene_snp, logP)


#update the phenotype names
rep_mat$Phenotype[rep_mat$Phenotype=="COVIDHOSPITALIZED_2COVIDHOSP_1COVIDNOHOSP"] <- "Hospitalized/Not_Hospitalized"
rep_mat$Phenotype[rep_mat$Phenotype=="COVIDHOSPITALIZEDPOP_2COVIDHOSP_1POPULATION"] <- "Hospitalized/Unscreened"
rep_mat$Phenotype[rep_mat$Phenotype=="COVIDSWABTEST_2COVIDPOS_1COVIDNEG"] <- "Positive/Negative"
rep_mat$Phenotype[rep_mat$Phenotype=="COVIDSWABTESTPOP_2COVIDPOS_1POPULATION"] <- "Positive/Unscreened"
rep_mat$Phenotype[rep_mat$Phenotype=="COVIDWITHPOSITIVEHOUSEMATE_2COVIDPOS_1COVIDNEG"] <- "Exposed_Positive/Exposed_Negative"
rep_mat$Phenotype[rep_mat$Phenotype=="COVIDWITHPOSITIVEHOUSEMATEPOP_2POPULATION_1COVIDNEGWITHEXP"] <- "Unscreened/Exposed_Negative"
rep_mat$Phenotype[rep_mat$Phenotype=="MILDASYMPTOMATICCOVID_2SYMPTOMATIC_1MILDASYMPTOMATIC"] <- "Symptomatic/Paucisymptomatic"
rep_mat$Phenotype[rep_mat$Phenotype=="COVID_SEVERITY_CONTINUOUS_PC"] <- "Continuous_Severity_Score"


#update the genotype names
rep_mat <- rep_mat %>% rename("ABO (rs505922)"="ABO-9:136149229:T:C") %>%
  #rename("ACSF3 (rs4782327)"="ACSF3-16:89184135:G:C") %>%
  rename("CCHCR1 (rs143334143)"="CCHCR1-6:31121426:G:A") %>%
  #rename("CCNG1 (rs79833209)"="CCNG1-5:162727453:C:T") %>%
  rename("DPP9 (rs2277732)"="DPP9-19:4723670:C:A") %>%
  #rename("FPR1 (rs12461764)"="FPR1-19:52242750:G:T") %>%
  rename("IFNAR2 (rs13050728)"="IFNAR2-21:34615210:T:C") %>%
  rename("IGF1 (rs10860891)"="IGF1-12:103014757:C:A") %>%
  rename("KAT7 (rs9903642)"="KAT7-17:47868921:G:A") %>%
  rename("SLC6A20/LZTFL1 (rs35081325)"="LZTFL1-3:45889921:A:T") %>%
  rename("OAS3 (rs7310667)"="OAS3-12:113392182:A:G") %>%
  rename("SLC6A20/LZTFL1 (rs73062389)"="SLC6A20-3:45835417:G:A") %>%
  rename("SLC6A20/LZTFL1 (rs2531743)"="SLC6A20-3:45838300:G:A") %>%
  rename("STM2A (rs622568)"="STM2A/VSTM2A-7:54647894:A:C") %>%
  rename("TMPRSS2 (rs2298661)"="TMPRSS2-21:42845642:C:A")


row.names(rep_mat) <- rep_mat[,1]
rep_mat <- rep_mat[,-1]
rep_mat <- data.matrix(rep_mat)

newnames <- lapply(
  rownames(rep_mat),
  function(x) bquote(italic(.(x))))


library(pheatmap)
p <- pheatmap(rep_mat, display_numbers = T, color = colorRampPalette(c('white', 'red'))(7), 
              fontsize_number = 10, scale="none",
              labels_row = as.expression(newnames))

ggsave(output_path, "figure2B_ReplicationHeatmap.eps", p, device="eps")
