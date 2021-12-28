#!/usr/bin/env Rscript
#Script to re-format GWAS summary stats output so that they can be read by METAL for trans-ethnic meta analysis
#use format Rscript [path/to/sumstats.logistic] [path/to/population_allele.frq] [path/to/outputfilename.metal]


args = commandArgs(trailingOnly=TRUE)

#load libraries
library("dplyr")

#load in the summarry stats file
infile_path = args[1]
eur <- read.table(infile_path, header=T, stringsAsFactors = F, comment.char="") %>% rename(CHR=X.CHROM)

#get A1 and A2.  A1 is the effect allele
eur <- eur %>% mutate(A2=ifelse(ALT1==A1, REF, ALT1))

#eur <- eur %>% select(CHR, SNP=ID, BP=POS, A1, A2, OR, P, N=OBS_CT, STAT=Z_STAT, SE=LOG.OR._SE)
eur <- eur %>% select(CHR, SNP=ID, BP=POS, A1, A2, beta=BETA, P, N=OBS_CT, STAT=Z_STAT, SE) %>%
  mutate(OR = exp(beta))

#Keep only the things we need and output new file
eur <- eur %>% select(CHR, SNP, BP, A1, A2, OR, P, N, STAT, beta, SE)

#add population-specific allele frequencies in
afs_infile <- args[2]
afs <- read.table(afs_infile, header=T, stringsAsFactors=F, comment.char="") %>% rename(CHR=X.CHROM, SNP=ID, A1=ALT, MAF=ALT_FREQS)
eur <-left_join(eur, afs, by="SNP")
eur <- eur %>% mutate(ALLELE_FREQ=ifelse(A1.x==A1.y, MAF, 1-MAF)) %>% select(CHR.x, SNP, BP, A1.x, A2, OR, beta, P, N, SE, ALLELE_FREQ)
colnames(eur) <- c('CHR','MarkerName','BP','Allele1','Allele2','OR','Effect','P-value','N','StdErr','Freq1')
#get rid of 'Inf' SNPs
eur <- eur %>% filter(!is.infinite(Effect))

outfile_path <- args[3]#usually defined with the suffix _input_sumstats_forMETAL.txt 
write.table(eur, outfile_path, row.names=F, quote=F, sep="\t")

#Run METAL trans-ethnic meta analysis using the script run_metal.sh 