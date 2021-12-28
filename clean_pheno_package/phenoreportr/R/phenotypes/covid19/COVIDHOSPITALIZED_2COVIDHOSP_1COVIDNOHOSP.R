import(magrittr)

clean_phenotype = function(survey) {
  ####Tested Positive and Hospitalized Vs Tested Positive and Not Hospitalized
  covid = survey %>% 
    dplyr::filter(question=="619.survey.covid19_study.have_you_been_swab_tested_for_covid19_commonly_referred_to_as_coronaviru") %>%
    dplyr::mutate(response1=gsub(" ", "", response)) %>%
    dplyr::select(iid, question1=question, response1, related, split, platform)
  
  hospitalized = survey %>%
    dplyr::filter(question=="623.survey.covid19_study.were_you_hospitalized_due_to_these_symptoms") %>%
    dplyr::mutate(response2=gsub(" ", "", response)) %>%
    dplyr::select(iid, question2=question, response2)
  
  hosp_covid = dplyr::left_join(covid, hospitalized, by="iid")
  hosp_covid = hosp_covid %>% 
    dplyr::mutate(fid=0) %>%
    dplyr::mutate(question=paste0(question1, "&", question2)) %>%
    dplyr::mutate(response=paste0(response1, "&", response2)) %>%
    dplyr::mutate(COVIDHOSPITALIZED_2COVIDHOSP_1COVIDNOHOSP=ifelse(response1=="Yes,andwaspositive" & response2=="Yes", 2,
                                                            ifelse(response1=="Yes,andwaspositive" & response2=="No", 1, -9))) %>%
    dplyr::select(fid, iid, question, response, related, split, platform, COVIDHOSPITALIZED_2COVIDHOSP_1COVIDNOHOSP) %>%
    dplyr::rename(FID=fid, IID=iid)
  
  return(hosp_covid)
}

