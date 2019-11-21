
#### Analysis Functions #####


sim_get_paired_or <- function(in_data,var_name){
  var_name <- paste0(var_name)
  
  in_data <- in_data %>% 
    rename(cond_ind := !!var_name)
  
  
  
  temp_mod <- tryCatch(survival::clogit(cond_ind~case + strata(strata), data=in_data),
                       warning = function(w) {NA},
                       error = function(e) {NA})
  
  
  mod_or <- tryCatch(exp(coef(temp_mod)[[1]]),
                     warning = function(w) {NA},
                     error = function(e) {NA})
  
  mod_or_p <- tryCatch(coef(summary(temp_mod))[,5],
                       warning = function(w) {NA},
                       error = function(e) {NA})

  
  return(tibble(or=mod_or,
                or_p_val=mod_or_p))
  
}


get_paired_or_est <- function(in_data,var_name){
  var_name <- paste0(var_name)
  
  in_data <- in_data %>% 
    rename(cond_ind := !!var_name)
  
  
  
  temp_mod <- tryCatch(survival::clogit(cond_ind~case + strata(strata), data=in_data),
                       warning = function(w) {NA},
                       error = function(e) {NA})
  
  temp_mod2 <- tryCatch(glm(cond_ind~case, data=in_data, family = poisson),
                        warning = function(w) {NA},
                        error = function(e) {NA})
  
  temp_conf <- tryCatch(exp(confint(temp_mod)),
                        warning = function(w) {c(NA,NA)},
                        error = function(e) {c(NA,NA)})
  
  temp_conf2 <- tryCatch(exp(confint(temp_mod2)),
                         warning = function(w) {matrix(NA,2,2)},
                         error = function(e) {matrix(NA,2,2)})
  
  
  temp_counts <- in_data %>%
    group_by(case) %>% 
    summarise(n=n(),
              exposure=sum(cond_ind),
              odds=exposure/(n-exposure),
              risk=exposure/n)
  
  mod_rr <- tryCatch(exp(coef(temp_mod2)[[2]]),
                     warning = function(w) {NA},
                     error = function(e) {NA})
  
  mod_or <- tryCatch(exp(coef(temp_mod)[[1]]),
                     warning = function(w) {NA},
                     error = function(e) {NA})
  
  mod_or_p <- tryCatch(coef(summary(temp_mod))[,5],
                       warning = function(w) {NA},
                       error = function(e) {NA})
  
  mod_rr_p <- tryCatch(coef(summary(temp_mod2))[2,4],
                       warning = function(w) {NA},
                       error = function(e) {NA})
  
  return(tibble(n_case=filter(temp_counts,case==1)$exposure,
                n_control=filter(temp_counts,case==0)$exposure,
                unpaired_or=filter(temp_counts,case==1)$odds/filter(temp_counts,case==0)$odds,
                unpaired_rr=filter(temp_counts,case==1)$risk/filter(temp_counts,case==0)$risk,
                or=mod_or,
                or_low=temp_conf[[1]],
                or_high=temp_conf[[2]],
                or_p_val=mod_or_p,
                rr=mod_rr,
                rr_low=temp_conf2[[2,1]],
                rr_high=temp_conf2[[2,2]],
                rr_p_val=mod_rr_p))
  
}

get_paired_or <- function(in_data,var_name){
  var_name <- paste0(var_name)
  
  in_data <- in_data %>% 
    rename(cond_ind := !!var_name)
  
  
  
  temp_mod <- tryCatch(survival::clogit(cond_ind~case + strata(strata), data=in_data),
                       warning = function(w) {NA},
                       error = function(e) {NA})
  
  temp_mod2 <- tryCatch(glm(cond_ind~case, data=in_data, family = poisson),
                        warning = function(w) {NA},
                        error = function(e) {NA})
  
  temp_conf <- tryCatch(exp(confint(temp_mod)),
                        warning = function(w) {c(NA,NA)},
                        error = function(e) {c(NA,NA)})
  
  temp_conf2 <- tryCatch(exp(confint(temp_mod2)),
                         warning = function(w) {matrix(NA,2,2)},
                         error = function(e) {matrix(NA,2,2)})
  
  
  temp_counts <- in_data %>%
    group_by(case) %>% 
    summarise(n=n(),
              exposure=sum(cond_ind),
              odds=exposure/(n-exposure),
              risk=exposure/n)
  
  mod_rr <- tryCatch(exp(coef(temp_mod2)[[2]]),
                     warning = function(w) {NA},
                     error = function(e) {NA})
  
  mod_or <- tryCatch(exp(coef(temp_mod)[[1]]),
                     warning = function(w) {NA},
                     error = function(e) {NA})
  
  mod_or_p <- tryCatch(coef(summary(temp_mod))[,5],
                       warning = function(w) {NA},
                       error = function(e) {NA})
  
  mod_rr_p <- tryCatch(coef(summary(temp_mod2))[2,4],
                       warning = function(w) {NA},
                       error = function(e) {NA})
  
  return(tibble(n_case=filter(temp_counts,case==1)$exposure,
                n_control=filter(temp_counts,case==0)$exposure,
                unpaired_or=filter(temp_counts,case==1)$odds/filter(temp_counts,case==0)$odds,
                unpaired_rr=filter(temp_counts,case==1)$risk/filter(temp_counts,case==0)$risk,
                or=mod_or,
                or_low=temp_conf[[1]],
                or_high=temp_conf[[2]],
                or_p_val=mod_or_p,
                rr=mod_rr,
                rr_low=temp_conf2[[2,1]],
                rr_high=temp_conf2[[2,2]],
                rr_p_val=mod_rr_p))
  
}

get_paired_or_icd <- function(icd9_list,icd10_list,cohort_data=study_cohort,time_map_data=time_map,con){
  temp <- get_dx_links(icd_9_codes = icd9_list,
                       icd_10_codes = icd10_list,
                       collect_tab = collect_table(dbs = "ccae"),
                       db.con = con)
  
  left_join(cohort_data,
            inner_join(time_map_data,
                       select(temp,key),
                       by="key") %>% 
              distinct(enrolid) %>% 
              mutate(dx_ind=1L),
            by="enrolid") %>% 
    mutate(dx_ind=ifelse(is.na(dx_ind),0L,dx_ind)) %>% 
    get_paired_or("dx_ind")
}

make_sim_data <- function(){
  cf_cases <- select(cf_study_cohort,enrolid,first_age,sex,enrmon,case) %>% 
    filter(case==1)
  
  cftr_cases <- filter(carrier_study_cohort,case==1)
  
  
  cftr_strata <- carrier_study_cohort %>% 
    filter(case==1) %>% 
    count(first_age,sex,enrmon) %>% 
    rename(n_cftr=n) %>% 
    arrange(first_age,sex,enrmon) %>% 
    mutate(strata=row_number())
  
  cf_strata <- cf_cases %>% 
    count(first_age,sex,enrmon) %>% 
    rename(n_cf=n)
  
  strata_match <- cftr_strata %>% 
    left_join(cf_strata,by=c("first_age","sex","enrmon")) %>% 
    mutate(n_cf=ifelse(is.na(n_cf),0,n_cf))
  
  
  #### Perfect_Match ###
  
  
  #### Case 1: cf_num >= cftr_num ####
  temp_case1 <- filter(strata_match,n_cftr<=n_cf)
  
  temp_case_draw1 <- filter(cf_study_cohort,case==1) %>% 
    select(-strata) %>% 
    inner_join(temp_case1, by = c("first_age","sex","enrmon")) %>% 
    select(strata,enrolid,first_age,sex,enrmon)
  
  temp_control_draw1 <- filter(cf_study_cohort,case==0) %>%
    select(-strata) %>%
    inner_join(temp_case1, by = c("first_age","sex","enrmon")) %>%
    select(strata,enrolid,first_age,sex,enrmon)
  
  
  ##### Case 2: cftr_num > cf_num > 0 ####
  
  temp_case2 <- filter(strata_match,n_cftr>n_cf & n_cf>0)
  
  ### Note: this approach resamples CF cases....we may need to select control cases without CF
  temp_case_draw2 <- temp_case2 %>% 
    inner_join(filter(cf_study_cohort,case==1) %>% 
                 select(enrolid,first_age,sex,enrmon),
               by = c("first_age","sex","enrmon")) %>% 
    mutate(dif=((n_cftr-n_cf) %/% n_cf)+1) %>% 
    group_by(strata,dif,n_cftr) %>% 
    nest() %>% 
    mutate(expand=map2(dif,data,~.y[rep(1:nrow(.y),each=(.x+1)),])) %>% 
    select(strata,n_cftr,expand) %>% 
    unnest() %>% 
    select(strata,enrolid,first_age,sex,enrmon)
  
  temp_control_draw2 <- temp_case2 %>%
    inner_join(filter(cf_study_cohort,case==0) %>%
                 select(enrolid,first_age,sex,enrmon),
               by = c("first_age","sex","enrmon")) %>%
    group_by(strata,n_cftr,n_cf) %>%
    nest() %>%
    mutate(n_control=n_cf*5,
           draws=ifelse(n_control>=n_cftr,1,((n_cftr-n_control) %/% n_control) + 2)) %>%
    mutate(expand=map2(draws,data,~.y[rep(1:nrow(.y),each=(.x)),])) %>%
    select(strata,n_cftr,expand) %>%
    unnest() %>%
    select(strata,enrolid,first_age,sex,enrmon)
  
  
  #### Case 3: cf_num==0 ####
  #### These cases we draw from controls ####
  temp_case3 <- filter(strata_match,n_cf==0)
  
  
  match_diff <- temp_case3 %>% 
    left_join(select(cf_cases, enrolid, sex, enrmon, first_age_cf=first_age)) %>% 
    mutate(yr_dif=abs(first_age-first_age_cf)) %>% 
    group_by(strata) %>% 
    summarise(match_dist=min(yr_dif)) %>% 
    inner_join(temp_case3)
  
  # threshold to increase the number of matches
  temp_thresh <- 0
  
  temp_match <- match_diff %>% 
    left_join(select(cf_cases, enrolid, sex, enrmon, first_age_cf=first_age)) %>% 
    filter(first_age_cf>=first_age-(match_dist+temp_thresh) & first_age_cf<=first_age+(match_dist+temp_thresh)) %>% 
    group_by(strata) %>% 
    mutate(n_cf_new=n()) %>% 
    ungroup() 
  
  # strata where match enough cases
  temp_case_draw3_1 <- temp_match %>% 
    filter(n_cf_new>=n_cftr) %>% 
    select(strata,enrolid,first_age=first_age_cf,sex,enrmon)
  
  temp_control_draw3_1 <- temp_case_draw3_1 %>% 
    distinct(strata,first_age,sex,enrmon) %>% 
    inner_join(filter(cf_study_cohort,case==0) %>%
                 select(-strata)) %>%
    select(strata,enrolid,first_age,sex,enrmon) 
  
  # strata where not enough
  temp_case_draw3_2 <- temp_match %>% 
    filter(n_cf_new<n_cftr) %>% 
    mutate(dif=((n_cftr-n_cf_new) %/% n_cf_new)+1) %>% 
    group_by(strata,dif,n_cftr) %>% 
    nest() %>% 
    mutate(expand=map2(dif,data,~.y[rep(1:nrow(.y),each=(.x+1)),])) %>% 
    select(strata,n_cftr,expand) %>% 
    unnest() %>% 
    select(strata,enrolid,first_age=first_age_cf,sex,enrmon)
  
  
  temp_control_draw3_2 <- temp_case_draw3_2 %>%
    select(-enrolid) %>% 
    left_join(filter(cf_study_cohort,case==0) %>%
                select(-strata)) %>%
    select(strata,enrolid,first_age,sex,enrmon)
  
  # strata where no match
  temp_case_draw3_3 <- match_diff %>% 
    filter(is.na(match_dist)) %>% 
    select(strata,first_age,sex,enrmon) %>% 
    inner_join(select(carrier_study_cohort,-strata)) %>% 
    filter(case==1) %>% 
    select(strata,enrolid,first_age,sex,enrmon)
  
  temp_control_draw3_3 <- match_diff %>% 
    filter(is.na(match_dist)) %>% 
    select(strata,first_age,sex,enrmon) %>% 
    inner_join(select(carrier_study_cohort,-strata)) %>% 
    filter(case==0) %>% 
    select(strata,enrolid,first_age,sex,enrmon)
  
  
  ##### Combine groups ####
  control_draw <- bind_rows(temp_control_draw1,
                            temp_control_draw2,
                            temp_control_draw3_1,
                            temp_control_draw3_2,
                            temp_control_draw3_3)
  
  case_draw <- bind_rows(temp_case_draw1,
                         temp_case_draw2,
                         temp_case_draw3_1,
                         temp_case_draw3_2,
                         temp_case_draw3_3)
  
  rm(temp_case_draw1,
     temp_case_draw2,
     temp_case_draw3_1,
     temp_case_draw3_2,
     temp_case_draw3_3,
     temp_control_draw1,
     temp_control_draw2,
     temp_control_draw3_1,
     temp_control_draw3_2,
     temp_control_draw3_3)
  
  rm(temp_case1,temp_case2,temp_case3)
  
  ##### Pops #####
  
  control_pop <- filter(carrier_study_cohort,case==0) %>% 
    select(enrolid,first_age,sex,enrmon,case) %>% 
    inner_join(cftr_strata,by = c("first_age", "sex", "enrmon"))
  
  cftr_pop <- filter(carrier_study_cohort,case==1) %>% 
    select(enrolid,first_age,sex,enrmon,case) %>% 
    inner_join(cftr_strata,by = c("first_age", "sex", "enrmon"))
  
  cftr_pop <- cftr_pop %>% 
    group_by(strata) %>% 
    mutate(strata_num=row_number()) %>% 
    ungroup()
  
  return(list(cftr_pop=cftr_pop,
              control_pop=control_pop,
              case_draw=case_draw,
              control_draw=control_draw))
}

#### Plotting Functions ####

capitalize <- function(string){
  string_len <- nchar(string)
  out <- paste0(toupper(str_sub(string,1,1)),tolower(str_sub(string,2,string_len)))
  out[out=="NANA"] <- NA
  return(out)
}

rem_z <- function(words){
  temp_func <- function(word){
    if (is.na(word)) {
      NA
    } else if (str_sub(word,1,1)=="z") {
      str_sub(word,2,1000)
    } else {
      word
    }
  }
  map(words,temp_func) %>% 
    unlist()
}


gen_plot1_log <- function(data){
  
  main_vals <- data %>% 
    arrange(group,label) %>% 
    mutate(n_case=ifelse(name!="infert_male",1000*(n_case/10380),1000*(n_case/1733)),
           n_control=ifelse(name!="infert_male",1000*(n_control/51898),1000*(n_control/8665))) %>% 
    select(label:n_control,or:or_p_val) %>% 
    mutate(or_p_val=as.character(sprintf("%.3f", round(or_p_val,3))),
           or_p_val=ifelse(or_p_val=="0.000","<0.001",or_p_val)) %>% 
    mutate(order=2)
  
  group_headings <- data  %>% 
    distinct(group) %>% 
    mutate(label=NA,
           n_case=NA,
           n_control=NA,
           or=NA,
           or_low=NA,
           or_high=NA,
           or_p_val=NA,
           order=1) 
  
  plot_data <-  bind_rows(main_vals,group_headings) %>% 
    arrange(group,order,label) %>% 
    mutate(label=rem_z(label)) %>% 
    rename(condition=label,
           p_val=or_p_val) %>% 
    mutate(group=ifelse(!is.na(condition),NA,group)) %>% 
    select(group,everything(),-order) %>%  # fix other code and remove this
    mutate(condition=as.factor(condition)) %>% 
    mutate(condition=capitalize(as.character(condition))) %>% 
    mutate(condition=ifelse(condition=="NA",NA,condition)) %>% 
    mutate_at(vars(n_case:or_high),funs(round(.,2))) %>%
    mutate(log_or=log(or),
           log_or_low=log(or_low),
           log_or_high=log(or_high),
           or_high_new=ifelse(log_or_high>5,5,log_or_high),
           or_text=paste0(sprintf("%.2f", round(or,2)),'  (',sprintf("%.2f", round(or_low,2)),'-',sprintf("%.2f", round(or_high,2)),')')) %>% 
    mutate(condition=ifelse(condition=="Cancer of colon, stomach, gi organs",
                            "Cancer of colon, stomach, GI organs",
                            ifelse(condition=="Gerd",
                                   "GERD",
                                   ifelse(condition=="Diabetes (type 1 or secondary)",
                                          "Diabetes (Type 1 or secondary)",
                                          condition))))
  
  rowseq <- seq(nrow(plot_data),1)
  par(mai=c(1,0,0,0), family = "Times")
  plot(plot_data$log_or, rowseq, pch=15,
       xlim=c(-9,9), ylim=c(1,nrow(plot_data)+1.5),
       xlab='', ylab='', yaxt='n', xaxt='n',
       bty='n')
  axis(1, seq(-1,5,by=1), cex.axis=1)
  
  segments(0,-1,0,nrow(plot_data), lty=3)
  segments(plot_data$log_or_low, rowseq, plot_data$or_high_new, rowseq)
  
  #mtext('Off-Pump\nCABG Better',1, line=2.5, at=0, cex=.5, font=2)
  #mtext('On-Pump\nCABG Better',1.5, line=2.5, at=2, cex=.5, font=2)
  
  text(-9,nrow(plot_data)+2.5, "Organ System & Condition", cex=1.2, font=2, pos=4)
  
  t1h <- ifelse(!is.na(plot_data$group), plot_data$group, '')
  text(-9,rowseq, t1h, cex=1.2, pos=4, font=2)
  t1 <- ifelse(!is.na(plot_data$condition), plot_data$condition, '')
  text(-8.5,rowseq, t1, cex=1.1, pos=4)
  
  text(-4.3,nrow(plot_data)+3.25, "Prevalence (per 1,000 enrollees)", cex=1.2, font=2, pos=4) 
  text(-4.3,nrow(plot_data)+2, "CF Carriers", cex=1.2, font=2, pos=4)
  text(-4.15,nrow(plot_data)+.85, "N = 10,380", cex=1.1, font=1, pos=4)
  t2 <- ifelse(!is.na(plot_data$n_case), format(plot_data$n_case,big.mark=","), '')
  text(-3.9, rowseq, t2, cex=1.1, pos=4)
  
  text(-2.6,nrow(plot_data)+2, "Matched Cohort", cex=1.2, font=2, pos=4)
  text(-2.25,nrow(plot_data)+1, "n = 51,898", cex=1.1, font=1, pos=4)
  t2 <- ifelse(!is.na(plot_data$n_control), format(plot_data$n_control,big.mark=","), '')
  text(-2.1, rowseq, t2, cex=1.1, pos=4)
  
  text(1,nrow(plot_data)+2.5, "LN Odds Ratio", cex=1.2, font=2, pos=4)
  
  text(4.9,nrow(plot_data)+2.5, "Odds Ratio (95% CI)", cex=1.2, font=2, pos=4)
  t3 <- ifelse(!is.na(plot_data$or),plot_data$or_text,'')
  text(5.1,rowseq, t3, cex=1.1, pos=4)
  
  text(7.8,nrow(plot_data)+2.5, "P-Value", cex=1.2, font=2, pos=4)
  t2 <- ifelse(!is.na(plot_data$p_val), format(plot_data$p_val,big.mark=","), '')
  text(8.9, rowseq, t2, cex=1.1, pos=2)
}


gen_all_cond_plot_log <- function(data,study_cohort_data){
  
  tab_counts <- study_cohort_data %>% 
    summarise(tot_case=sum(case==1),
              tot_control=sum(case==0),
              tot_male_case=sum(case==1 & sex==1),
              tot_male_control=sum(case==0 & sex==1))
  
  main_vals <- data %>% 
    arrange(group,label) %>% 
    mutate(n_case=ifelse(name!="infert_male",
                         1000*(n_case/tab_counts$tot_case),
                         1000*(n_case/tab_counts$tot_male_case)),
           n_control=ifelse(name!="infert_male",
                            1000*(n_control/tab_counts$tot_control),
                            1000*(n_control/tab_counts$tot_male_control))) %>% 
    select(label:n_control,or:or_p_val) %>% 
    mutate(or_p_val=as.character(sprintf("%.3f", round(or_p_val,3))),
           or_p_val=ifelse(or_p_val=="0.000","<0.001",or_p_val)) %>% 
    mutate(order=2)
  
  group_headings <- data  %>% 
    distinct(group) %>% 
    mutate(label=NA,
           n_case=NA,
           n_control=NA,
           or=NA,
           or_low=NA,
           or_high=NA,
           or_p_val=NA,
           order=1) 
  
  plot_data <-  bind_rows(main_vals,group_headings) %>% 
    arrange(group,order,label) %>% 
    mutate(label=rem_z(label)) %>% 
    rename(condition=label,
           p_val=or_p_val) %>% 
    mutate(group=ifelse(!is.na(condition),NA,group)) %>% 
    select(group,everything(),-order) %>%  # fix other code and remove this
    mutate(condition=as.factor(condition)) %>% 
    mutate(condition=capitalize(as.character(condition))) %>% 
    mutate(condition=ifelse(condition=="NA",NA,condition)) %>% 
    mutate_at(vars(n_case:or_high),funs(round(.,2))) %>%
    mutate(log_or=ifelse(log(or)>-1,log(or),-1),
           log_or_low=ifelse(log(or_low)>-1,log(or_low),-1),
           log_or_high=log(or_high),
           or_high_new=ifelse(log_or_high>5,5,log_or_high),
           or_text=paste0(sprintf("%.2f", round(or,2)),'  (',sprintf("%.2f", round(or_low,2)),'-',sprintf("%.2f", round(or_high,2)),')')) %>% 
    mutate(condition=ifelse(condition=="Cancer of colon, stomach, gi organs",
                            "Cancer of colon, stomach, GI organs",
                            ifelse(condition=="Gerd",
                                   "GERD",
                                   ifelse(condition=="Diabetes (type 1 or secondary)",
                                          "Diabetes (Type 1 or secondary)",
                                          condition))))
  
  rowseq <- seq(nrow(plot_data),1)
  par(mai=c(1,0,0,0), family = "Times")
  plot(plot_data$log_or, rowseq, pch=15,
       xlim=c(-9,8), ylim=c(1,nrow(plot_data)+1.5),
       xlab='', ylab='', yaxt='n', xaxt='n',
       bty='n')
  axis(1, seq(-1,4,by=1), cex.axis=1)
  
  segments(0,-1,0,nrow(plot_data), lty=3)
  segments(plot_data$log_or_low, rowseq, plot_data$or_high_new, rowseq)
  
  #mtext('Off-Pump\nCABG Better',1, line=2.5, at=0, cex=.5, font=2)
  #mtext('On-Pump\nCABG Better',1.5, line=2.5, at=2, cex=.5, font=2)
  
  text(-9,nrow(plot_data)+2.5, "Organ System & Condition", cex=1.2, font=2, pos=4)
  
  t1h <- ifelse(!is.na(plot_data$group), plot_data$group, '')
  text(-9,rowseq, t1h, cex=1.2, pos=4, font=2)
  t1 <- ifelse(!is.na(plot_data$condition), plot_data$condition, '')
  text(-8.5,rowseq, t1, cex=1.1, pos=4)
  
  text(-4.3,nrow(plot_data)+3.25, "Prevalence (per 1,000 enrollees)", cex=1.2, font=2, pos=4) 
  text(-4.3,nrow(plot_data)+2, "CF Carriers", cex=1.2, font=2, pos=4)
  text(-4.15,nrow(plot_data)+1, paste0("N = ",formatC(tab_counts$tot_case,big.mark=",")), cex=1.1, font=1, pos=4)
  t2 <- ifelse(!is.na(plot_data$n_case), format(plot_data$n_case,big.mark=","), '')
  text(-3.9, rowseq, t2, cex=1.1, pos=4)
  
  text(-2.6,nrow(plot_data)+2, "Matched Cohort", cex=1.2, font=2, pos=4)
  text(-2.25,nrow(plot_data)+1, paste0("N = ",formatC(tab_counts$tot_control,big.mark=",")), cex=1.1, font=1, pos=4)
  t2 <- ifelse(!is.na(plot_data$n_control), format(plot_data$n_control,big.mark=","), '')
  text(-2.1, rowseq, t2, cex=1.1, pos=4)
  
  #text(1,nrow(plot_data)+2.5, expression("log"[e]*" Odds Ratio"), cex=1.2, font=2, pos=4)
  text(1,nrow(plot_data)+2.5, "$log_{e}$ Odds Ratio", cex=1.2, font=2, pos=4)
  
  text(3.8,nrow(plot_data)+2.5, "Odds Ratio (95% CI)", cex=1.2, font=2, pos=4)
  t3 <- ifelse(!is.na(plot_data$or),plot_data$or_text,'')
  text(4.1,rowseq, t3, cex=1.1, pos=4)
  
  text(6.6,nrow(plot_data)+2.5, "P-Value", cex=1.2, font=2, pos=4)
  text(6.3,nrow(plot_data)+1.5, "(One-Sided)", cex=1.2, font=2, pos=4)
  t2 <- ifelse(!is.na(plot_data$p_val), format(plot_data$p_val,big.mark=","), '')
  text(7.6, rowseq, t2, cex=1.1, pos=2)
} 


gen_all_cond_plot_log_cf <- function(data,study_cohort_data){
  
  tab_counts <- study_cohort_data %>% 
    summarise(tot_case=sum(case==1),
              tot_control=sum(case==0),
              tot_male_case=sum(case==1 & sex==1),
              tot_male_control=sum(case==0 & sex==1))
  
  main_vals <- data %>% 
    arrange(group,label) %>% 
    mutate(n_case=ifelse(name!="infert_male",
                         1000*(n_case/tab_counts$tot_case),
                         1000*(n_case/tab_counts$tot_male_case)),
           n_control=ifelse(name!="infert_male",
                            1000*(n_control/tab_counts$tot_control),
                            1000*(n_control/tab_counts$tot_male_control))) %>% 
    select(label:n_control,or:or_p_val) %>% 
    mutate(or_p_val=as.character(round(or_p_val,3)),
           or_p_val=ifelse(or_p_val=="0","<0.001",or_p_val)) %>% 
    mutate(order=2)
  
  group_headings <- data  %>% 
    filter(!(group %in% c("Control","NEW"))) %>% 
    distinct(group) %>% 
    mutate(label=NA,
           n_case=NA,
           n_control=NA,
           or=NA,
           or_low=NA,
           or_high=NA,
           or_p_val=NA,
           order=1) 
  
  plot_data <-  bind_rows(main_vals,group_headings) %>% 
    arrange(group,order,label) %>%
    mutate(label=rem_z(label)) %>%
    rename(condition=label,
           p_val=or_p_val) %>% 
    mutate(group=ifelse(!is.na(condition),NA,group)) %>% 
    select(group,everything(),-order) %>%  # fix other code and remove this
    mutate(condition=as.factor(condition)) %>% 
    mutate(condition=capitalize(as.character(condition))) %>% 
    mutate(condition=ifelse(condition=="NA",NA,condition)) %>% 
    mutate_at(vars(n_case:or_high),funs(round(.,2))) %>%
    mutate(log_or=log(or),
           log_or_low=log(or_low),
           log_or_high=log(or_high),
           or_high_new=ifelse(log_or_high>9,9,log_or_high),
           or_text=paste0(sprintf("%.2f", round(or,2)),'  (',sprintf("%.2f", round(or_low,2)),'-',sprintf("%.2f", round(or_high,2)),')')) %>% 
    mutate(condition=ifelse(condition=="Cancer of colon, stomach, gi organs",
                            "Cancer of colon, stomach, GI organs",
                            ifelse(condition=="Gerd",
                                   "GERD",
                                   ifelse(condition=="Diabetes (type 1 or secondary)",
                                          "Diabetes (Type 1 or secondary)",
                                          condition))))
  
  rowseq <- seq(nrow(plot_data),1)
  par(mai=c(1,0,0,0), family = "Times")
  plot(plot_data$log_or, rowseq, pch=15,
       xlim=c(-11,14.5), ylim=c(1,nrow(plot_data)+1.5),
       xlab='', ylab='', yaxt='n', xaxt='n',
       bty='n')
  axis(1, seq(-1,9,by=1), cex.axis=1)
  
  segments(0,-1,0,nrow(plot_data), lty=3)
  segments(plot_data$log_or_low, rowseq, plot_data$or_high_new, rowseq)
  
  #mtext('Off-Pump\nCABG Better',1, line=2.5, at=0, cex=.5, font=2)
  #mtext('On-Pump\nCABG Better',1.5, line=2.5, at=2, cex=.5, font=2)
  
  text(-11.2,nrow(plot_data)+2.5, "Organ System & Condition", cex=1.2, font=2, pos=4)
  
  t1h <- ifelse(!is.na(plot_data$group), plot_data$group, '')
  text(-11.2,rowseq, t1h, cex=1.2, pos=4, font=2)
  t1 <- ifelse(!is.na(plot_data$condition), plot_data$condition, '')
  text(-10.7,rowseq, t1, cex=1.1, pos=4)
  
  text(-5,nrow(plot_data)+3, "Prevalence (per 1,000 enrollees)", cex=1.2, font=2, pos=4) 
  text(-4.84,nrow(plot_data)+2, "CF Cohort", cex=1.2, font=2, pos=4)
  text(-4.74,nrow(plot_data)+1, paste0("N = ",formatC(tab_counts$tot_case,big.mark=",")), cex=1.1, font=1, pos=4)
  t2 <- ifelse(!is.na(plot_data$n_case), format(plot_data$n_case,big.mark=","), '')
  text(-4.5, rowseq, t2, cex=1.1, pos=4)
  
  text(-2.94,nrow(plot_data)+2, "Matched Cohort", cex=1.2, font=2, pos=4)
  text(-2.48,nrow(plot_data)+1, paste0("N = ",formatC(tab_counts$tot_control,big.mark=",")), cex=1.1, font=1, pos=4)
  t2 <- ifelse(!is.na(plot_data$n_control), format(plot_data$n_control,big.mark=","), '')
  text(-2.2, rowseq, t2, cex=1.1, pos=4)
  
  text(4,nrow(plot_data)+2.5, "LN Odds Ratio", cex=1.2, font=2, pos=4)
  
  text(9.3,nrow(plot_data)+2.5, "Odds Ratio (95% CI)", cex=1.2, font=2, pos=4)
  t3 <- ifelse(!is.na(plot_data$or),plot_data$or_text,'')
  text(9.3,rowseq, t3, cex=1.1, pos=4)
  
  text(13.3,nrow(plot_data)+2.5, "P-Value", cex=1.2, font=2, pos=4)
  text(13,nrow(plot_data)+1.5, "(One-Sided)", cex=1.2, font=2, pos=4)
  t2 <- ifelse(!is.na(plot_data$p_val), format(plot_data$p_val,big.mark=","), '')
  text(15, rowseq, t2, cex=1.1, pos=2)
}

plot_main_cor <- function(conditions){
  temp_sort <- cftr_any_or %>% 
    filter(name %in% conditions) %>% 
    select(name) %>% 
    inner_join(cf_any_or) %>% 
    filter(group!="Control") %>% 
    select(name,label,group,or,or_low,or_high) %>% 
    arrange(desc(or)) %>% 
    mutate(row_n=row_number()) %>% 
    select(label,row_n)
  
  bind_rows(cf_any_or %>% 
              select(name,label,group,or,or_low,or_high) %>% 
              mutate(cohort="CF"),
            cftr_any_or %>% 
              select(name,label,group,or,or_low,or_high) %>% 
              mutate(cohort="CF Carrier")) %>% 
    mutate(cohort=fct_relevel(cohort,"CF Carrier")) %>% 
    mutate_at(vars(or:or_high),funs(log)) %>%
    inner_join(temp_sort) %>% 
    mutate(label=factor(label,levels=temp_sort$label)) %>% 
    ggplot(aes(x=label,y=or, color=cohort)) +
    geom_point(size=3) +
    geom_pointrange(aes(ymin=or_low,ymax=or_high),size=1) +
    coord_flip() +
    theme_minimal() +
    ylab("LN Odds Ratio") +
    xlab("") +
    scale_color_manual(values=c("royalblue3","red2")) +
    theme(text=element_text(family = "Times")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "grey")) +
    theme(legend.position="bottom") +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14)) +
    scale_color_manual(values=c("red","blue")) +
    theme(legend.title=element_blank())
}

plot_scatter <- function(conditions){
  bind_rows(cf_any_or %>% 
              select(name,label,group,or,or_low,or_high) %>% 
              mutate(cohort="cf"),
            cftr_any_or %>% 
              select(name,label,group,or,or_low,or_high) %>% 
              mutate(cohort="cftr")) %>% 
    filter(name %in% conditions) %>% 
    mutate_at(vars(or:or_high),funs(log)) %>%
    select(label,or,cohort) %>% 
    spread(key=cohort,value=or) %>% 
    ggplot(aes(cf,cftr)) +
    geom_point() +
    geom_smooth(method = "lm",se = FALSE,color="green4") +
    theme_minimal() +
    coord_fixed() +
    xlim(0,8.5) +
    ylim(0,8.5) +
    xlab("LN Odds Ratio for CF") +
    ylab("LN Odds Ratio for CF Carriers") +
    ggthemes::theme_few() +
    geom_abline(slope = 1, intercept = 0,linetype=2) +
    theme(text=element_text(family = "Times")) 
}

count_conditions <- function(condition_indicators,time_map,count_by_days=21){
  tmp_cond_names <- names(select(condition_indicators,-key))
  
  tmp_time_map <- time_map %>% 
    left_join(condition_indicators,by="key") %>% 
    mutate_at(vars(tmp_cond_names),funs(replace_na(.,0)))
  
  out_counts <- tmp_time_map %>% 
    mutate(admdate=admdate %/% count_by_days) %>% 
    select(enrolid,admdate,!!tmp_cond_names) %>% 
    group_by(enrolid,admdate) %>% 
    summarise_all(funs(sum)) %>% 
    ungroup() %>% 
    mutate_at(vars(-enrolid,-admdate),funs(ifelse(.>0,1,0))) %>% 
    select(-admdate) %>% 
    group_by(enrolid) %>% 
    summarise_all(funs(sum))
}


collapse_enrollment <- function(enrolid_list,vars,db.con,collect_tab=collect_table()){
  temp <- collect_tab %>%
    select(-setting) %>%
    mutate(data=map2(db,year,~tbl(db.con,paste0("enrollment_detail_",.x,"_",.y)) %>%
                       filter(enrolid %in% enrolid_list) %>%
                       select_(.dots = c("enrolid","dtstart","dtend",vars)) %>%
                       collect() %>%
                       mutate(enrolid=as.integer64(enrolid))))
  
  temp <- temp %>%
    select(data) %>%
    mutate(data=map(data,~mutate_at(.,vars(vars),funs(as.integer)))) %>%
    unnest()
  
  temp_strata <- temp %>%
    distinct_(.dots=c("enrolid",vars)) %>%
    mutate(strata=row_number())
  
  temp <- temp %>%
    inner_join(temp_strata, by = c("enrolid",vars))
  
  out <- temp %>%
    arrange(enrolid,dtstart) %>%
    group_by(enrolid) %>%
    mutate(gap=((dtstart-lag(dtend))>1) | strata!=lag(strata),
           gap=ifelse(is.na(gap),FALSE,gap)) %>%
    mutate(period=cumsum(gap)) %>%
    group_by_(.dots=c("enrolid","period",vars)) %>%
    summarise(dtstart=min(dtstart),
              dtend=max(dtend)) %>%
    ungroup()
  
  
  return(out)
}

