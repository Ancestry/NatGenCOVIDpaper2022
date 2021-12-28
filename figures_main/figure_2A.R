#!/usr/bin/env Rscript

#USE: Rscript figure_2A.R

library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

pwr <- read.csv("figure_2A_power_analysis_by_phenotype.csv")

#change to long format
pwr_long <- pwr %>% gather(key="MAF", value="Power", MAF.0.01:MAF.0.50) %>%
  mutate(MAF=str_remove(MAF, "MAF.")) %>%
  mutate(MAF=as.numeric(MAF))

pwr_long$Phenotype <- factor(pwr_long$Phenotype, levels = c("Positive/Unscreened",
                                                            "Positive/Negative",
                                                            "Unscreened/Exposed_Negative",
                                                            "Exposed_Positive/Exposed_Negative",
                                                            "Hospitalized/Unscreened", 
                                                            "Hospitalized/Not_Hospitalized",
                                                            "Symptomatic/Paucisymptomatic"
                                                            ))

pwr_long$Alpha <- factor(pwr_long$Alpha, levels=c("0.05", "5e-08"))
pwr_long <- pwr_long %>% rename(`Control Type`=Control_type)


#plot
rep_power <- pwr_long %>% filter(Alpha=="0.05")

p <- ggplot(rep_power, aes(x = MAF, y = Power, color=Phenotype)) + 
  scale_color_brewer(palette = "Paired") +
  geom_point() +
  geom_line(aes(linetype = `Control Type`)) + 
  #facet_wrap(~ factor(Alpha)) +
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")
        )

#display.brewer.all(colorblindFriendly = TRUE)
#brewer.pal(n = 8, name = "Paired")

ggsave("figure2A_power_plot.eps", 
       width = 7, height = 5, units="in",
       p, device="eps")
