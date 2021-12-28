import(magrittr)

clean_phenotype = function(trait, survey) {
  m = modules::use(here::here("R", "phenotypes", "covid19", paste0(trait, ".R")))
  m$clean_phenotype(survey)
}

#' Generate clean phenotype files
#'
#' @param traits Name of the trait(s).
#' @param survey_file Path to survey file.
#' @param pheno_path Path to phenotypes data directory.
#'
#' @export
clean_phenotypes = function(traits, survey_file, pheno_path) {
  survey = data.table::fread(survey_file, na.strings=c("","NA"))
  
  for (trait in traits) {
    data = clean_phenotype(trait, survey)
    utils::write.table(data, file.path(pheno_path, paste0(trait, ".tsv")), row.names=F, quote=F, sep="\t")
    utils::write.table(data %>% dplyr::filter(related==0 & split=="gwas" & !is.na(trait)), 
                       file.path(pheno_path, paste0(trait, ".pheno")), row.names=F, quote=F, sep="\t")
  }
}
