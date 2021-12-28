import(magrittr)

clean_phenotype = function(survey) {
  ###Tested Positive and Household Member Positive Vs. Tested Negative and Household Member Positive
  covid = survey %>% 
    dplyr::filter(question=="619.survey.covid19_study.have_you_been_swab_tested_for_covid19_commonly_referred_to_as_coronaviru") %>%
    dplyr::mutate(response1=gsub(" ", "", response)) %>%
    dplyr::select(iid, question1=question, response1, related, split, platform)
  
  house = survey %>%
    dplyr::filter(question=="636.survey.covid19_study.has_someone_in_your_household_tested_positive_for_covid19") %>%
    dplyr::mutate(response2=gsub(" ", "", response)) %>%
    dplyr::select(iid, question2=question, response2)
  
  house_covid = dplyr::left_join(covid, house, by="iid")
  house_covid = house_covid %>% 
    dplyr::mutate(fid=0) %>%
    dplyr::mutate(question=paste0(question1, "&", question2)) %>%
    dplyr::mutate(response=paste0(response1, "&", response2)) %>%
    dplyr::mutate(COVIDWITHPOSITIVEHOUSEMATEPOP_2POPULATION_1COVIDNEGWITHEXP=ifelse(response1=="Yes,andwasnegative" & response2=="Yes,atleastoneperson", 1, 2)) %>%
    dplyr::select(fid, iid, question, response, related, split, platform, COVIDWITHPOSITIVEHOUSEMATEPOP_2POPULATION_1COVIDNEGWITHEXP) %>%
    dplyr::rename(FID=fid, IID=iid)
  
  return(house_covid)
}
