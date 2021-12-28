#!/usr/bin/env Rscript

library(dplyr)
library(data.table)


#' Load ethnicity data.
#'
#' @param ethnicity_path   Path to ethnicity data directory.
#' @param default_cohorts  Default ethnicity cohorts.
#'
#' @export
#'
#' @importFrom data.table fread
#' @importFrom dplyr mutate select
#' 
#' @return Ethnicity data.
load_ethnicity = function(ethnicity_path, default_cohorts) {
  
  default_cohorts = unlist(strsplit(default_cohorts, ","))
  ethnicity_data = data.table(FID = character(), IID = character(), eth_cohort = character())
  
  for (cohort in default_cohorts) {
    cohort_str = unlist(strsplit(cohort, ":"))
    cohort_filename = cohort_str[1]
    cohort_column = cohort_str[2]
    cohort_filepath = file.path(ethnicity_path, cohort_filename, fsep = .Platform$file.sep)
    cohort_data = data.table::fread(cohort_filepath) %>% dplyr::mutate(eth_cohort = cohort_column) 
    ethnicity_data = rbind(cohort_data)
  }
  colnames(ethnicity_data) = c("FID", "IID", "eth_cohort")
  ethnicity_data = ethnicity_data %>% dplyr::select("IID", "eth_cohort")
  return(ethnicity_data)
}


#' Load phenotype data.
#'
#' @param pheno_path  Path to phenotypes data directory.
#' @param trait_name  Name of the trait.
#' @param ext         File extension (default .tsv).
#' 
#' @export
#' 
#' @importFrom data.table fread
#' 
#' @return Phenotype data.
load_pheno = function(pheno_path, trait_name, ext = ".tsv") {
  pheno_filepath = file.path(pheno_path, paste0(trait_name, ext), fsep = .Platform$file.sep)
  pheno_data = data.table::fread(pheno_filepath)
  return(pheno_data)
}


#' Extract questions from phenotype data.
#'
#' @param pheno_data  Phenotype data.
#' 
#' @export
#' 
#' @return Questions.
extract_questions = function(pheno_data) {
  questions = unlist(strsplit(as.character(pheno_data[1, "question"]), "&", fixed=TRUE))
  return(questions)
}


#' Load survey data and subset to questions of interest.
#'
#' @param survey_path  Path to survey file.
#' @param questions    Questions of interest.
#' 
#' @export
#' 
#' @importFrom data.table fread
#' @importFrom dplyr filter
#' 
#' @return Survey data.
load_survey = function(survey_path, questions) {
  survey_data = data.table::fread(survey_path) %>% dplyr::filter(question %in% questions)
  return(survey_data)
}


#' Load survey data, phenotype data, ethnicity data.
#'
#' @param survey_path      Path to survey file.
#' @param ethnicity_path   Path to ethnicity data directory.
#' @param pheno_path       Path to phenotypes data directory.
#' @param trait_name       Name of the trait.
#' @param default_cohorts  Default ethnicity cohorts.    
#' 
#' @export
#' 
#' @return Prepared data.
load_data = function(survey_path, ethnicity_path, pheno_path, trait_name, default_cohorts) {
  pheno_data = load_pheno(pheno_path = pheno_path, trait_name = trait_name)
  questions = extract_questions(pheno_data = pheno_data)
  ethnicity_data = load_ethnicity(ethnicity_path = ethnicity_path, default_cohorts = default_cohorts)
  survey_data = load_survey(survey_path = survey_path, questions = questions)
  
  names(survey_data)[names(survey_data) == "genotype_id"] <- "IID"
  
  joined_survey = dplyr::left_join(survey_data, pheno_data, by = "IID")
  joined_survey = dplyr::left_join(joined_survey, ethnicity_data, by = "IID")
  names(joined_survey)[names(joined_survey) == trait_name] <- "clean_pheno"
  joined_survey = joined_survey %>% dplyr::mutate(clean_cat_pheno = addNA(as.factor(ifelse(is.na(clean_pheno) | clean_pheno==-9, "not_assigned", clean_pheno))))
  joined_survey$eth_cohort[is.na(joined_survey$eth_cohort)] <- "not_assigned"
  return(list(data = joined_survey, questions = questions))
}