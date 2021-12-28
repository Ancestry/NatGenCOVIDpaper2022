#!/usr/bin/env Rscript

library(ggplot2)
library(RColorBrewer)


#' Generate sample size by response plot.
#'
#' @param data Prepared data.
#' @param questions Questions.
#'
#' @export
#'
#' @importFrom ggplot2 ggplot
#' 
#' @return Sample size by response plot.
phenotype_sample_size_by_response_plot = function(data, questions) {
  surveys_byq <- NULL
  p <- NULL
  for(q in questions){
    survey_by_question = data %>% filter(question.x == q)
    summ_table = table(survey_by_question$response.x, survey_by_question$clean_pheno, exclude = NULL)
    surveys_byq[[q]] = as.data.frame.matrix(summ_table)
    categorical_gradient = RColorBrewer::brewer.pal(3,"Set2")
    
    p[[q]] = ggplot2::ggplot(survey_by_question) + 
      geom_bar(aes_string(x = "response.x", y = "..count..", fill = "clean_cat_pheno"),
               position = "dodge")+
      guides(fill = guide_legend(title = "Cleaned Phenotype")) +
      scale_fill_manual(values = categorical_gradient) +
      theme_bw(base_size = 14) + ylab("Count") +
      ggtitle(paste(q)) +
      theme(axis.text.x = element_text(angle = 90, size = 10), plot.title = element_text(size = 10))
  }
  return(p)
}


#' Generate prevalence plot.
#'
#' @param data Prepared data.
#'
#' @export
#'
#' @importFrom ggplot2 ggplot
#' 
#' @return Prevalence plot.
prevalence_plot = function(data) {
  data = data %>%
    mutate(age_binned=cut_width(age.x, 10, boundary = 10)) %>%
    filter(sex_genetic %in% c("FEMALE", "MALE")) %>%
    filter(sex_genetic %in% c("FEMALE", "MALE") & !is.na(sex_genetic)) %>%
    filter(eth_cohort != "not_assigned")
  
  categorical_gradient = rev(brewer.pal(9,"Spectral"))
  p = ggplot2::ggplot(data) + geom_bar(aes_string(x = "clean_cat_pheno", y = "..prop..", fill = "age_binned", group = "age_binned"), position = "dodge") +
    guides(fill=guide_legend(title = "age group")) +
    facet_grid(eth_cohort ~ sex_genetic) +
    scale_fill_manual(values = categorical_gradient)
  return(p)
}
