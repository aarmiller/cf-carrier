

##### This script contains the basic R script that was used to identify and 
##### build the CF Carrier data and matched controls. The identical approach 
##### was used to build the CF comparison population along with matched controls
##### NOTE: This script contains a combination of actual and psuedo code. T

# After thius script has run the identified enrolee ids will need to be extracted from the 
# raw Truven data

#### load required packages
library(tidyverse)
library(bit64)

#### Load in Required Data ####

# You will need to load 2 data sets

# carrier_enrolids <- <This should contain all enrolids>


##### Get distinct enrolids in each group ####
carrier_enrolids <- distinct(carrier_enrolids,enrolid)
cf_enrolids_any <- distinct(cf_enrolids_any,enrolid)


#### Identify families with cf or cftr mutations ####
cftr_family <- bind_rows(enrolids,cf_enrolids_any) %>%
  mutate(efamid=enrolid %/% 100) %>%
  select(efamid) %>%
  distinct() %>%
  mutate(efamid=as.integer(efamid))

cftr_familial <- inner_join(cftr_family,ccae_enroll,by="efamid")

save(cf_cftr_familial,file = "<path to data>")


#### Identify distinct enrolids that do not have CF at some point ####
cftr_enrolids <- enrolids %>%
  anti_join(cf_enrolids_any)

cftr_familial <- cftr_familial %>%
  left_join(mutate(cftr_enrolids,carrier=1L),by="enrolid") %>%
  left_join(mutate(cf_enrolids_any,cf=1L),by="enrolid") %>%
  left_join(mutate(cftr_enrolids,
                   efamid=as.integer(enrolid %/% 100)) %>%
              distinct(efamid) %>%
              mutate(carrier_in_fam=1L),by="efamid") %>%
  left_join(mutate(cf_enrolids_any,
                   efamid=as.integer(enrolid %/% 100)) %>%
              distinct(efamid) %>%
              mutate(cf_in_fam=1L),by="efamid") %>%
  mutate_at(vars(carrier:cf_in_fam),funs(ifelse(is.na(.),0L,.)))

save(cftr_familial,file = "<path to data>")

carrier_enrolids <- cftr_enrolids

save(cftr_enrolids,cf_enrolids_any,cftr_familial,file = "<path to data>")

# build cases by joining with enrollment info
cases <- cftr_familial %>% filter(carrier==1)

#### compute distinct strata and number of observations per stata ####


## Matching on first age
strata <- cases %>% 
  count(first_age,sex,enrmon) %>% 
  rename(strata_size=n) %>% 
  mutate(strata=row_number())

#### MAKE PORENTIAL CONTROL COHORT ####

## Matching on first age
controls <- inner_join(ccae_enroll,
                       strata,
                       by=c("first_age","sex","enrmon")) %>% 
  anti_join(select(cftr_familial,enrolid),by="enrolid") %>%
  arrange(strata)

# seed to draw random controls
set.seed(1234)
controls <- controls %>% 
  mutate(rand=runif(nrow(controls))) %>% 
  group_by(strata) %>% 
  arrange(strata,rand) %>% 
  mutate(strata_num=row_number()) %>% 
  ungroup()

### draw a 1:5 match
temp <- controls %>% 
  filter(strata_num<=strata_size*5)

# final controls
controls <- temp %>% 
  mutate(match_id=((strata_num-1) %/% 5)+1) %>% 
  select(-rand,-strata_num,-strata_size)

## matching on first age
cases <- cases %>%
  inner_join(select(strata,-strata_size),by=c("first_age","sex","enrmon")) %>%
  arrange(strata) %>%
  group_by(strata) %>%
  mutate(match_id=row_number()) %>%
  ungroup()

study_cohort <- bind_rows(mutate(cases,case=1L),
                          mutate(controls,case=0L)) %>% 
  mutate_at(vars(carrier:cf_in_fam),funs(ifelse(is.na(.),0,.))) 

## update match_id so it is distinct 
new_match_id <- cases %>% 
  select(strata,match_id) %>% 
  arrange(strata,match_id) %>% 
  mutate(match_id_new=row_number())

study_cohort <- study_cohort %>% 
  inner_join(new_match_id,by=c("strata","match_id")) %>% 
  mutate(match_id=match_id_new) %>% 
  select(-match_id_new)

carrier_study_cohort <- study_cohort

save(carrier_study_cohort,file="<path>")

