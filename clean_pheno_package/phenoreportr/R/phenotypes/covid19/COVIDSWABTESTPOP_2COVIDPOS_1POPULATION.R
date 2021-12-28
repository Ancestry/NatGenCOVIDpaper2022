import(magrittr)

clean_phenotype = function(survey) {
  data = survey %>%
    dplyr::filter(question=="619.survey.covid19_study.have_you_been_swab_tested_for_covid19_commonly_referred_to_as_coronaviru") %>%
    dplyr::mutate(response=gsub(" ", "", response)) %>%
    dplyr::mutate(fid=0) %>%
    dplyr::mutate(COVIDSWABTESTPOP_2COVIDPOS_1POPULATION=ifelse(response=="Yes,andwaspositive", 2 ,1)) %>%
    dplyr::select(fid, iid, question, response, related, split, platform, COVIDSWABTESTPOP_2COVIDPOS_1POPULATION) %>%
    dplyr::rename(FID=fid, IID=iid)
  
  return(data)
}
