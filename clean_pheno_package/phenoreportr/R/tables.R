#!/usr/bin/env Rscript

library(dplyr)
library(stargazer)


#' Generate prevalence table from prepared data.
#'
#' @param data Prepared data.
#'
#' @export
#'
#' @importFrom dplyr group_by filter mutate select summarise
#' @importFrom stargazer stargazer
#' 
#' @return Prevalence table.
prevalence_table = function(data) {
  prev_cleaned = data %>% 
    dplyr::filter(!is.na(eth_cohort)) %>%
    dplyr::group_by(eth_cohort, clean_pheno) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(prevalence = round(n / sum(n), digits = 3)) %>%
    dplyr::mutate(cc_ratio = lag(n, n = 2L)/n) %>% 
    dplyr::filter(clean_pheno == 2) %>%
    dplyr::select(-clean_pheno, n_cases = n) %>%
    dplyr::mutate(n_controls = n_cases*cc_ratio, n_total=n_cases+n_controls) %>%
    dplyr::select(eth_cohort, n_cases, n_controls, n_total, "control:case_ratio" = cc_ratio, prevalence)
  stargazer::stargazer(prev_cleaned, summary = FALSE, type = "text", rownames = FALSE, digits = 2, title = "Cleaned Phenotype")
}


#' Generate response sample size by question from prepared data.
#'
#' @param data Prepared data.
#' @param questions Questions.
#'
#' @export
#'
#' @importFrom dplyr group_by filter mutate summarise
#' @importFrom stargazer stargazer
#' 
#' @return Response sample size by question table.
response_sample_size_by_question_table = function(data, questions) {
  for(i in questions){
    prev_response = data %>%
      dplyr::filter(question.x == i) %>%
      dplyr::filter(!is.na(eth_cohort)) %>%
      dplyr::group_by(eth_cohort, response.x) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::mutate(response_prevalence = round(n / sum(n), digits = 3))
    
    prev_response_cast <- dcast(prev_response, response.x ~ eth_cohort, value.var = "response_prevalence")
    n_response_cast <- dcast(prev_response, response.x ~ eth_cohort, value.var = "n")
    
    stargazer(prev_response_cast, summary = FALSE, type = "text", rownames = FALSE, digits = 2)
    stargazer(n_response_cast, summary = FALSE, type = "text", rownames = FALSE, digits = 2)
  }
}


#' Generate responses combined table from prepared data.
#'
#' @param data Prepared data.
#'
#' @export
#'
#' @importFrom dplyr group_by filter mutate summarise
#' @importFrom stargazer stargazer
#' 
#' @return Responses combined table.
response_combined_table = function(data) {
  prev_response_all <- data %>%
    dplyr::filter(!is.na(eth_cohort)) %>%
    dplyr::group_by(eth_cohort, response.y) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(response_prevalence = round(n / sum(n), digits = 3))
  
  prev_response_all_cast <- dcast(prev_response_all, response.y ~ eth_cohort, value.var = "response_prevalence")
  n_response_all_cast <- dcast(prev_response_all, response.y ~ eth_cohort, value.var = "n")
  
  stargazer(prev_response_all_cast, summary = FALSE, type = "text", rownames = FALSE, digits = 2)
  stargazer(n_response_all_cast, summary = FALSE, type = "text", rownames = FALSE, digits = 2)
}






