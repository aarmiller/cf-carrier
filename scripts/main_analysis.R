

rm(list=ls()) 
library(bit64)
library(tidyverse)
library(icd)
library(survival)
library(icdplus)
source("functions/small_db_functions.R")
source("functions/analysis_functions.R")

#########################################
#### Connect to the cohort databases ####
#########################################

# Note: After running the build_mathced_cohorts.R database we extracted all enrollee ids
#       and added these enrollee's data to a SQLite database.

# carrier database
db_con_carrier <- src_sqlite("<base_path>/carrier.db")

# cf database
db_con_cf <- src_sqlite("<base_path>/carrier.db")

# collect carrier study cohort information
carrier_study_cohort <- tbl(db_con_carrier,"carrier_study_cohort") %>% collect(n=Inf)

# collect cf study cohort information
cf_study_cohort <- tbl(db_con_cf,"cf_study_cohort") %>% collect(n=Inf)


#####################################
#### load conditions of interest ####
#####################################

source("functions/condition_codes.R")
condition_codes <- get_cond_codes()

cond_names <- unlist(map(condition_codes,~.$include)) %>% 
  .[.] %>% 
  names()

cond_table <- tibble(name=cond_names) %>% 
  mutate(label=map_chr(name,~condition_codes[[.]]$label),
         group=map_chr(name,~condition_codes[[.]]$group)) %>% 
  mutate(icd9_codes=map(name,~condition_codes[[.]]$icd9_codes)) %>% 
  mutate(icd10_codes=map(name,~condition_codes[[.]]$icd10_codes)) 

#########################################################################
#### Create indicators for conditions of interest in carrier cohorts ####
#########################################################################

#### For carriers ####

# Get diagnosis keys for carriers
carrier_dx_keys <- get_dx_links(icd_9_codes = cond_table$icd9_codes %>% unlist() %>% unique(),
                                icd_10_codes = cond_table$icd10_codes %>% unlist() %>% unique(),
                                collect_tab = collect_table(dbs = "ccae"),
                                db.con = db_con_carrier)

# Label diagnosis keys for carriers
carrier_dx_keys <- bind_rows(cond_table %>% 
                               select(name,dx=icd9_codes) %>% 
                               unnest(),
                             cond_table %>% 
                               select(name,dx=icd10_codes) %>% 
                               unnest()) %>% 
  distinct() %>% 
  inner_join(carrier_dx_keys)

# Create indicators for carriers
carrier_dx_inds <- carrier_dx_keys %>% 
  select(-dx) %>% 
  distinct() %>% 
  mutate(ind=1L) %>% 
  spread(key=name,value=ind,fill = 0L)

#### For CF ####

# Get diagnosis keys for carriers
cf_dx_keys <- get_dx_links(icd_9_codes = cond_table$icd9_codes %>% unlist() %>% unique(),
                           icd_10_codes = cond_table$icd10_codes %>% unlist() %>% unique(),
                           collect_tab = collect_table(dbs = "ccae"),
                           db.con = db_con_cf)

# Label diagnosis keys for carriers
cf_dx_keys <- bind_rows(cond_table %>% 
                               select(name,dx=icd9_codes) %>% 
                               unnest(),
                             cond_table %>% 
                               select(name,dx=icd10_codes) %>% 
                               unnest()) %>% 
  distinct() %>% 
  inner_join(cf_dx_keys)

# Create indicators for carriers
cf_dx_inds <- cf_dx_keys %>% 
  select(-dx) %>% 
  distinct() %>% 
  mutate(ind=1L) %>% 
  spread(key=name,value=ind,fill = 0L)


########################################################################
#### Create longitudinal dataset of visits and condition indicators ####
########################################################################

#### For Carriers ####

# create longitudinal time maps
carrier_time_map <- make_time_map(db_con_carrier) %>% 
  left_join(distinct(carrier_study_cohort,enrolid,case),by="enrolid")

enrolid_cond_counts <- carrier_time_map %>% 
  inner_join(carrier_dx_inds) %>% 
  group_by(enrolid) %>% 
  summarise_at(vars(abdominal_pain:venous_thromboembolism),funs(sum)) %>% 
  mutate_at(vars(abdominal_pain:venous_thromboembolism),funs(ifelse(.>0,1,0)))

enrolid_cond_inds <- carrier_study_cohort %>% 
  left_join(enrolid_cond_counts) %>% 
  mutate_at(vars(abdominal_pain:venous_thromboembolism),funs(ifelse(is.na(.),0,.)))


#### For CF ####

# create longitudinal time maps
cf_time_map_ <- make_time_map(db_con_cf) %>% 
  left_join(distinct(cf_study_cohort,enrolid,case),by="enrolid")

enrolid_cond_counts_cf <- cf_time_map %>% 
  inner_join(cf_dx_inds) %>% 
  group_by(enrolid) %>% 
  summarise_at(vars(abdominal_pain:venous_thromboembolism),funs(sum)) %>% 
  mutate_at(vars(abdominal_pain:venous_thromboembolism),funs(ifelse(.>0,1,0)))

enrolid_cond_inds_cf <- cf_study_cohort %>% 
  left_join(enrolid_cond_counts_cf) %>% 
  mutate_at(vars(abdominal_pain:venous_thromboembolism),funs(ifelse(is.na(.),0,.)))


##############################################
#### Estimate relative risk / odds ratios ####
##############################################


res_carrier <- tibble(cond=enrolid_cond_inds %>% 
         select(abdominal_pain:venous_thromboembolism) %>% 
         names(.)) %>% 
  #slice(1:2) %>%  # uncomment to test
  mutate(or_est=map(cond,~get_paired_or(enrolid_cond_inds,.))) %>% 
  unnest()


res_cf <- tibble(cond=enrolid_cond_inds_cf %>% 
                        select(abdominal_pain:venous_thromboembolism) %>% 
                        names(.)) %>% 
  #slice(1:2) %>%  # uncomment to test
  mutate(or_est=map(cond,~get_paired_or(enrolid_cond_inds,.))) %>% 
  unnest()

# save results for plotting

save(res_carrier,res_cf,file = "<path to file>")



##################################################################
#### Compute counts of number of conditions and organ systems ####
##################################################################

# function to count the number of conditions per enrollee
get_cond_count_or <- function(ind_data,low_list,high_list){
  temp <- ind_data %>% 
    select_(.dots=c("enrolid","case","strata",paste0(cond_names))) %>% 
    gather(key=name,value = value,-enrolid,-case,-strata) %>% 
    group_by(enrolid,case,strata) %>% 
    summarise(n=sum(value==TRUE)) 
  
  temp_func <- function(low,high){
    temp_fit <- temp %>% 
      mutate(cond_ind=ifelse(between(n,low,high),1,0)) %>% 
      survival::clogit(cond_ind~case + strata(strata), data=.)
    
    return(tibble(or=exp(coef(temp_fit)[[1]]),
                  or_lo=exp(confint(temp_fit))[1],
                  or_hi=exp(confint(temp_fit))[2]))
  }
  
  out_dat <- tibble(low=low_list,
                    high=high_list) %>% 
    mutate(data=map2(low,high,~temp_func(.x,.y))) %>% 
    unnest()
  
  return(out_dat)
}

# count of conditions for carriers
cftr_cond_count_dat <- get_cond_count_or(enrolid_cond_inds,
                                         c(0,1,3,5,7,9,11),
                                         c(0,2,4,6,8,10,100))

# count of conditions for cf
cf_cond_count_dat <- get_cond_count_or(enrolid_cond_inds_cf,
                                       c(0,1,3,5,7,9,11),
                                       c(0,2,4,6,8,10,100))


# function to count the number of organ systems per enrollee
get_org_count_or<- function(ind_data,low_list,high_list){
  temp <- ind_data %>% 
    select_(.dots=c("enrolid","case","strata",paste0(cond_names))) %>% 
    gather(key=name,value = value,-enrolid,-case,-strata) %>% 
    mutate(name=str_replace(name,"_ind","")) %>% 
    inner_join(select(include_table,name,group),by="name") %>% 
    mutate(temp=ifelse(value==TRUE,group,NA)) %>% 
    group_by(enrolid,case,strata) %>% 
    summarise(n=n_distinct(temp,na.rm = TRUE)) 
  
  temp_func <- function(low,high){
    temp_fit <- temp %>% 
      mutate(cond_ind=ifelse(between(n,low,high),1,0)) %>% 
      survival::clogit(cond_ind~case + strata(strata), data=.)
    
    return(tibble(or=exp(coef(temp_fit)[[1]]),
                  or_lo=exp(confint(temp_fit))[1],
                  or_hi=exp(confint(temp_fit))[2]))
  }
  
  out_dat <- tibble(low=low_list,
                    high=high_list) %>% 
    mutate(data=map2(low,high,~temp_func(.x,.y))) %>% 
    unnest()
  
  return(out_dat)
}


cftr_org_count_dat <- get_org_count_or(enrolid_cond_inds,
                                       c(0:4,5),
                                       c(0:4,100))

cf_org_count_dat <- get_org_count_or(enrolid_cond_inds_cf,
                                     c(0:4,5),
                                     c(0:4,100))

