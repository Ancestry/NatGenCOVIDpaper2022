import(magrittr)

clean_phenotype = function(survey) {
  ###Tested Positive Vs. Tested Negative
  covid <- survey %>% 
    dplyr::filter(question=="619.survey.covid19_study.have_you_been_swab_tested_for_covid19_commonly_referred_to_as_coronaviru") %>%
    dplyr::mutate(response=gsub(" ", "", response)) %>%
    dplyr::mutate(fid=0) %>%
    dplyr::mutate(COVIDSWABTEST_2COVIDPOS_1COVIDNEG=ifelse(response=="Yes,andwaspositive", 2,
                                                    ifelse(response=="Yes,andwasnegative",1, -9))) %>%
    dplyr::select(fid, iid, question, response, related, split, platform, COVIDSWABTEST_2COVIDPOS_1COVIDNEG) %>%
    dplyr::rename(FID=fid, IID=iid)
  
  return(covid)
}
