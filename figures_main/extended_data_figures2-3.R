#!/usr/bin/env Rscript

library(data.table)
library(tidyverse)
library(ggplot2)

##### generate manhattan plots and QQplots

pheno_list <- c("Continuous_Severity_Score",
"ExposedPositive_ExposedNegative",
"Hospitalized_NotHosptialized",
"Hospitalized_Unscreened",
"Positive_Negative",
"Positive_Unscreened",
"Symptomatic_Paucisymptomatic",
"Unscreened_ExposedNegative")

data_path <- "./"
file_suffix <- ".TransAncestryMeta.METAL.SumStats.csv.gz"
out_path_meta <- "./visualization/"

for(pheno in pheno_list){
  
phenotype_name=pheno
meta_path <- meta_path <- paste0(data_path, phenotype_name,file_suffix)
out_path_meta <- paste0(out_path_meta, phenotype_name,"_TransAncestryMeta_")

# read file
imp_meta <- fread(meta_path, header = T, stringsAsFactors = F)

#split out SNP name in g1 to get chrom and position
g1 <- imp_meta %>% mutate(SNP=MarkerName) %>% separate(MarkerName, c("CHR", "BP", "A1_temp", "A2_temp"), ":")
rm(imp_meta)
g1$CHR <- ifelse(g1$CHR=="X", 23, g1$CHR)
g1$CHR <- as.numeric(g1$CHR)
g1$BP <- as.numeric(g1$BP)
g1$P <- as.numeric(g1$P)

#join the data together for plotting
joint <- g1 %>% select(CHR, SNP, BP, A1=A1_temp, A2=A2_temp, beta=Effect, P=`Pvalue`, Direction, HetPVal)

# Filter out any SNPs that are missing a result from one or more study in the meta
joint <- joint %>% filter(!stringr::str_detect(Direction, '\\?'))


#do manhattan plotting
joint = joint[order(joint$CHR,joint$BP),]
joint$plot = c(1:nrow(joint))
midpoints=
  unlist(lapply(c(1:23),
                function(i){
                  return(mean(joint$plot[which(joint$CHR==i)]))
                }))

joint_filtered <- joint %>% filter(CHR>0 & CHR<24) %>% 
  filter(P <= 0.05)

#make the actual plot
p1=ggplot(joint_filtered, aes(x=plot, y=-1*log10(P)))+
  geom_point(aes(color=factor(CHR%%2)), size=0.5)+
  scale_colour_manual(values = c("#48C9B0", "#138D75"))+
  geom_hline(yintercept=-log10(5E-08), colour="#990000", linetype="dashed")+
  geom_hline(yintercept=-log10(1E-05), colour="#051094", linetype="dashed")+
  scale_x_continuous(breaks=midpoints,labels=c(1:23))+
  ylab(expression(paste("-log"[10], plain("P-value"))))+
  xlab("Chromosome")+
  guides(color=FALSE)+
  #ylim(-log10(0.05), -log10(pmin(min(joint_filtered$P, na.rm=T, 5E-08)))+0.1)+
  geom_text(aes(label=ifelse(P<5E-08,as.character(SNP),'')),hjust=0,vjust=0, size=2)+
  ggtitle(paste0(phenotype_name, " Trans-Ancestry Meta Analysis"))+
  theme_classic()+
  theme(axis.text=element_text(size=6))

#save the output  
out_manh_path_meta <- paste0(out_path_meta, "manhattanPlot.png")
ggsave(plot=p1, filename=out_manh_path_meta, width=8, height=3.5, units="in", device="png")



#####make a qqplot for the imputed GWAS and output lambda

#compute GIF
inflation <- function(ps) {
  ps <- ps[!is.na(ps)]
  z = qnorm(ps / 2)
  lambda = round(median(z^2) / 0.454, 3)
  lambda
}

#write a plotting function
gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}

#loop over the plotting function and add on the GIF annotation
ps <- joint$P[!is.na(joint$P)]
q1 <- gg_qqplot(ps) +
  ggtitle(paste0(phenotype_name, " Trans-Ancestry Meta Analysis"))+
  theme_bw(base_size = 10) +
  annotate(
    geom = "text",
    x = -Inf,
    y = Inf,
    hjust = -0.15,
    vjust = 1 + 0.15 * 3,
    label = sprintf("λ = %.2f", inflation(ps)),
    size = 8
  ) +
  theme(
    axis.ticks = element_line(size = 0.5),
    panel.grid = element_blank()
  )

#save the output
out_qq_path_meta <- paste0(out_path_meta, "qqPlot.png")
ggsave(plot=q1, filename=out_qq_path_meta, width=7, height=4, units="in", device="png")

}
