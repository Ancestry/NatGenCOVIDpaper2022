library(dplyr)
library(stringr)
library(tidyr)
library(LDlinkR) #gets 1KG LD info https://cran.r-project.org/web/packages/LDlinkR/vignettes/LDlinkR.html
library(RColorBrewer)



c3_pvals <- read.table("Desktop/pvals_chr3only.csv", header=T, sep=",", stringsAsFactors = F)
c3_pvals <- c3_pvals %>% arrange(Position)
snps <- unique(c3_pvals$rsID)

#you will need to get your own LDlink token at https://ldlink.nci.nih.gov/?tab=apiaccess if you wish to repeat this analysis
LDlink_token=""

#get some LD quantities
eur_r2 <- LDmatrix(snps, pop = "CEU", r2d = "r2", token = LDlink_token, file = FALSE) %>% 
  arrange(match(RS_number, c3_pvals$rsID))
eur_dp <- LDmatrix(snps, pop = "CEU", r2d = "d", token = LDlink_token, file = FALSE) %>% 
  arrange(match(RS_number, c3_pvals$rsID))


rownames(eur_r2) <- unique(c3_pvals$Label.Name)
positions <- unique(c3_pvals$Position)


#Plot LDmap with LDheatmap
library(LDheatmap) #https://cran.r-project.org/web/packages/LDheatmap/vignettes/LDheatmap.pdf

#make r^2 heatmap
rtest <- data.matrix(eur_r2)
rtest <- rtest[1:5,2:6]
rtest_heatmap <- LDheatmap(rtest, flip = TRUE, genetic.distances = positions, 
                           SNP.name=colnames(rtest),
                           color=rev(brewer.pal(n = 9, name = "Reds")),
                           add.map = T,
                           title = "",
                           geneMapLabelY = 2,
                           geneMapLocation = 0.1)
#add genes to the plot
rtest_plusgenes <- LDheatmap.addGenes(rtest_heatmap, chr="chr3", genome="hg19", splice_variants=F)



#make D' heatmap
dtest <- data.matrix(eur_dp)
dtest <- dtest[1:5,2:6]
dtest_heatmap <- LDheatmap(dtest, flip = TRUE, genetic.distances = positions, 
                           LDmeasure = "D'",
                           SNP.name=colnames(rtest),
                           add.map = T,
                           title = "",
                           geneMapLabelY = 2,
                           geneMapLocation = 0.1,
                           color=rev(brewer.pal(n = 9, name = "Blues")))


#making the ggplot of P-values
library(ggplot2)

#order the phenotypes
c3_pvals$Phenotype <- factor(c3_pvals$Phenotype, levels = c(levels = c("Positive/Unscreened",
                                                                       "Positive/Negative",
                                                                       "Unscreened/Exposed_Negative",
                                                                       "Exposed_Positive/Exposed_Negative",
                                                                       "Hospitalized/Unscreened", 
                                                                       "Hospitalized/Not_Hospitalized",
                                                                       "Symptomatic/Paucisymptomatic",
                                                                       "Continuous_Severity_Score")))



#order the label names
c3_pvals$Label.Name <- factor(c3_pvals$Label.Name, levels=c("rs73062389 - HGI C2 EUR Susceptibility (Oct 2020)",
                                                            "rs2271616 - HGI C2 Meta Susceptibility (Jan 2021)",
                                                            "rs2531743 - Horowitz et al. Susceptibility (Dec 2020)",
                                                            "rs10490770 - HGI B2 Meta Severity (Jan 2021)", 
                                                            "rs35081325 - HGI B2 EUR Severity (Oct 2020)"))




#make the -log10(P-value) plot
colors <- brewer.pal(n = 8, name = "Paired")


chr3_pval_plot <- ggplot(c3_pvals, aes(x = Label.Name, y = -log10(Pval), fill = Phenotype)) + 
  geom_line(aes(group = Phenotype, color=Phenotype))+
  geom_point(shape = 21, aes(size=MAF)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12))+
  geom_hline(yintercept = -log10(5E-08), linetype='dotted', col = 'red')+
  geom_hline(yintercept = -log10(0.05), linetype='dotted', col = 'black') +
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)+
  xlab("")+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

chr3_pval_plot

ggsave("supp_fig6_nogridlines.eps", 
       #width = 8, height = 3, units="in",
       chr3_pval_plot, device="eps")


