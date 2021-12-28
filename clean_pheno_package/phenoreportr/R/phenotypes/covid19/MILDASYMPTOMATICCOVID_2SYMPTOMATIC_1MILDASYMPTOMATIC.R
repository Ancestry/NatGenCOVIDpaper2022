import(magrittr)
import(stringr)

clean_phenotype = function(survey) {
  ###Tested Positive and mildly symptomatic or asymptomatic versus tested postive and symptomatic
  covid = survey %>% 
    dplyr::filter(question=="619.survey.covid19_study.have_you_been_swab_tested_for_covid19_commonly_referred_to_as_coronaviru") %>%
    dplyr::mutate(question=substr(question, 1, 60)) %>%
    dplyr::mutate(response1=gsub(" ", "", response)) %>%
    dplyr::mutate(response1=substr(response1, 1, 20)) %>%
    dplyr::select(iid, question1=question, response1, related, split, platform)
  
  asym = survey %>%
    dplyr::filter(str_detect(question, "620.")) %>%
    dplyr::mutate(question=substr(question, 1, 60)) %>%
    dplyr::mutate(response2=gsub(" ", "", response)) %>%
    dplyr::mutate(response2=substr(response2, 1, 20)) %>%
    dplyr::select(iid, question2=question, response2)
  
  sym = survey %>% 
    dplyr::filter(question=="621.survey.covid19_study.between_the_beginning_of_february_now_have_you_had_any_of_the_following_") %>%
    dplyr::mutate(symptom_group=ifelse(response %in% c("None", "Very mild", "Mild"), "None/Mild",
                                ifelse(response %in% c("Moderate", "Severe", "Very severe"), "Moderate/Severe", NA))) %>%
    dplyr::mutate(concat_response=paste0(symptom_group, "-", matrix_row))
  
  #they scored all 15 symptoms and they reported that all 15 symptoms were either none, very mild, or mild
  sym = sym %>% 
    dplyr::group_by(iid, symptom_group) %>%
    dplyr::summarise(count = n() ) %>%
    dplyr::mutate(mild_only=ifelse(count==15, "mild_only", "some_moderate_severe")) %>%
    dplyr::filter(symptom_group=="None/Mild") %>% 
    dplyr::mutate(question3="621...have_you_had_any_of_the_following")%>%
    dplyr::select(iid, question3, mild_only) %>%
    data.frame()
  
  #join the data frames together and make the final phenotype
  asym_covid = dplyr::left_join(covid, asym, by="iid")
  asym_covid = dplyr::left_join(asym_covid, sym, by="iid")
  
  asym_covid = asym_covid %>% 
    dplyr::mutate(fid=0) %>%
    dplyr::mutate(question=paste0(question1, "&", question2, "&", question3)) %>%
    dplyr::mutate(response=paste0(response1, "&", response2, "&", mild_only)) %>%
    dplyr::mutate(MILDASYMPTOMATICCOVID_2SYMPTOMATIC_1MILDASYMPTOMATIC=ifelse(response=="Yes,andwaspositive&Yes&NA" | response=="Yes,andwaspositive&Yes&some_moderate_severe", 2,
                                                                       ifelse(response=="Yes,andwaspositive&No&NA" | response=="Yes,andwaspositive&&mild_only" | response=="Yes,andwaspositive&Notsure&mild_only" | response=="Yes,andwaspositive&Yes&mild_only", 1, -9))) %>%
    dplyr::select(fid, iid, question, response, related, split, platform, MILDASYMPTOMATICCOVID_2SYMPTOMATIC_1MILDASYMPTOMATIC) %>%
    dplyr::rename(FID=fid, IID=iid)
  
  return(asym_covid)
}
