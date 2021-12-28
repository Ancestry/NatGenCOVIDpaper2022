#!/usr/bin/env Rscript

"
phenoreportr command-line interface

Generate clean phenotypes.
Generate phenotype report.

Usage:
    phenoreport.R process phenotypes --survey=<file> --pheno=<path> --trait=<name>...
    phenoreport.R process report --survey=<file> --pheno=<path> --ethnicity=<path> --trait=<name>... 
                  [--cohorts=<default>] [--format=<type>] [--dir=<path>] 
    phenoreport.R list phenotypes
    phenoreport.R -h | --help
    phenoreport.R -v | --version

Options:
    -e, --ethnicity=<path>   Path to ethnicity cohorts.
    -f, --format=<type>      The R Markdown output format to convert to [default: html_document]
    -d, --dir=<path>         The output directory for the rendered output_file.
    -p, --pheno=<path>       Phenotype path.
    -s, --survey=<file>      Path to survey file.
    -t, --trait=<name>       Name of the trait.
    -c, --cohorts=<default>  Default cohorts [default: Europe.idvs:EUR,African_American.idvs:AA,Latino.idvs:LAT,East_Asia.idvs:EAS]
    -h, --help               Show this help message.
    -v, --version            Show package version.
" -> doc

library(docopt)
library(here)
library(phenoreportr)
library(modules)


phenoreportr_version <- as.character(packageVersion("phenoreportr"))
cmdargs <- docopt(doc, version = phenoreportr_version)


#' phenoreportr CLI
#'
#' phenoreportr command-line interface processor.
#'
#' @param cmdargs List of command-line arguments.
#'
#' @examples
cli <- function(cmdargs) {
  
  if (cmdargs$process & cmdargs$report) {
    
    default_cohorts = unlist(strsplit(cmdargs$cohorts, split=","))
    
    if (is.null(cmdargs$dir)) {
      dir = getwd()
    }
    else {
      dir = cmdargs$dir
    }
    
    for (trait_name in cmdargs$trait) {
      rmarkdown::render(here::here("R", "report.Rmd"),
                        params = list(
                          survey_file = cmdargs$survey,
                          pheno_path = cmdargs$pheno,
                          ethnicity_path = cmdargs$ethnicity,
                          trait_name = trait_name,
                          default_cohorts = cmdargs$cohorts
                        ),
                        output_format = cmdargs$format,
                        output_file = trait_name,
                        output_dir = dir)
    }
  }
  else if (cmdargs$process & cmdargs$phenotypes) {
    m = modules::use(here::here("R", "phenotypes", "clean_phenotype.R"))
    m$clean_phenotypes(traits=cmdargs$trait, survey_file=cmdargs$survey, pheno_path=cmdargs$pheno)
  }
  else if (cmdargs$list & cmdargs$phenotypes) {
    filenames = list.files(here::here("R", "phenotypes", "covid19"), pattern="\\.R$")
    basenames = sub(pattern="\\.R$", replacement="", x=filenames)
    print(basenames)
  }
}

cli(cmdargs)


# To find CLI executable:
# system.file("exec", "phenoreport.R", package = "phenoreportr")
