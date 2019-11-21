


###### New Functions to use with small databases #####

collect_table <- function(setting="inpatient",years=NA,dbs=c("ccae","mdcr"),HCUP=FALSE) {
  if (is.na(years)){
    years <- str_pad(c(1:17),2,pad="0")
  }
  if (is.na(setting)) {
    tibble(db=rep(dbs,each=length(years)),
           year=rep(years,times=length(dbs)))
  } else {
    tibble(setting=setting,
           db=rep(dbs,each=length(years)),
           year=rep(years,times=length(dbs)))
  }
}

make_time_map <- function(db=db_con){
  dat <- rbind(db %>% 
                 tbl("outpatient_keys") %>%
                 collect(n=Inf) %>%
                 mutate(disdate=svcdate,
                        inpatient=0L) %>% 
                 select(key,year,ccae,enrolid,admdate=svcdate,disdate,inpatient,stdplac),
               db %>% 
                 tbl("inpatient_keys") %>% 
                 select(-caseid) %>% 
                 collect(n=Inf) %>% 
                 mutate(inpatient=1L,stdplac=-1L) %>% 
                 select(key,everything())) %>% 
    arrange(enrolid,admdate,inpatient)
  
  return(dat)
}

get_dx_links <- function(icd_9_codes,icd_10_codes,db.con=db_con,collect_tab=collect_table()){
  
  ##### get inpatient diagnoses ####
  in_temp1 <- collect_tab %>% 
    filter(as.integer(year)<15) %>% 
    mutate(data=map2(db,year,
                     ~tbl(db.con,paste0("inpatient_dx_",.x,"_",.y)) %>% 
                       filter(dx %in% icd_9_codes) %>% 
                       distinct(caseid,dx) %>% 
                       inner_join(tbl(db.con,paste0("inpatient_core_",.x,"_",.y)) %>% 
                                    select(caseid,enrolid),
                                  by="caseid") %>% 
                       distinct(enrolid,caseid,dx) %>% 
                       collect(n=Inf) %>% 
                       mutate(enrolid=as.integer64(enrolid))))
  
  in_temp2 <- collect_tab %>% 
    filter(as.integer(year)>14) %>% 
    mutate(data=map2(db,year,
                     ~tbl(db.con,paste0("inpatient_dx9_",.x,"_",.y)) %>% 
                       filter(dx %in% icd_9_codes) %>% 
                       distinct(caseid,dx) %>%
                       inner_join(tbl(db.con,paste0("inpatient_core_",.x,"_",.y)) %>% 
                                    select(caseid,enrolid),
                                  by="caseid") %>% 
                       distinct(enrolid,caseid,dx) %>% 
                       collect(n=Inf) %>% 
                       mutate(enrolid=as.integer64(enrolid))))
  
  in_temp3 <- collect_tab %>% 
    filter(as.integer(year)>14) %>% 
    mutate(data=map2(db,year,
                     ~tbl(db.con,paste0("inpatient_dx10_",.x,"_",.y)) %>% 
                       filter(dx %in% icd_10_codes) %>% 
                       distinct(caseid,dx) %>%
                       inner_join(tbl(db.con,paste0("inpatient_core_",.x,"_",.y)) %>% 
                                    select(caseid,enrolid),
                                  by="caseid") %>% 
                       distinct(enrolid,caseid,dx) %>% 
                       collect(n=Inf) %>% 
                       mutate(enrolid=as.integer64(enrolid))))
  
  in_temp <- bind_rows(in_temp1,in_temp2,in_temp3) %>% 
    select(-setting) %>% 
    unnest() %>% 
    group_by(db,year) %>% 
    nest()
  
  rm(in_temp1,in_temp2,in_temp3)
  
  #### get outpatient diagnoses ####
  out_temp1 <- collect_tab %>% 
    filter(as.integer(year)<15) %>% 
    mutate(data=map2(db,year,
                     ~tbl(db.con,paste0("outpatient_dx_",.x,"_",.y)) %>% 
                       filter(dx %in% icd_9_codes) %>% 
                       distinct(seqnum_o,enrolid,svcdate,dx) %>% 
                       inner_join(tbl(db.con,paste0("outpatient_core_",.x,"_",.y)) %>% 
                                    select(seqnum_o,stdplac),
                                  by="seqnum_o") %>% 
                       distinct(enrolid,svcdate,stdplac,dx) %>% 
                       collect(n=Inf) %>% 
                       mutate(enrolid=as.integer64(enrolid))))
  
  out_temp2 <- collect_tab %>% 
    filter(as.integer(year)>14) %>% 
    mutate(data=map2(db,year,
                     ~tbl(db.con,paste0("outpatient_dx9_",.x,"_",.y)) %>% 
                       filter(dx %in% icd_9_codes) %>% 
                       distinct(seqnum_o,enrolid,svcdate,dx) %>% 
                       inner_join(tbl(db.con,paste0("outpatient_core_",.x,"_",.y)) %>% 
                                    select(seqnum_o,stdplac),
                                  by="seqnum_o") %>% 
                       distinct(enrolid,svcdate,stdplac,dx) %>% 
                       collect(n=Inf) %>% 
                       mutate(enrolid=as.integer64(enrolid))))
  
  out_temp3 <- collect_tab %>% 
    filter(as.integer(year)>14) %>% 
    mutate(data=map2(db,year,
                     ~tbl(db.con,paste0("outpatient_dx10_",.x,"_",.y)) %>% 
                       filter(dx %in% icd_10_codes) %>% 
                       distinct(seqnum_o,enrolid,svcdate,dx) %>% 
                       inner_join(tbl(db.con,paste0("outpatient_core_",.x,"_",.y)) %>% 
                                    select(seqnum_o,stdplac),
                                  by="seqnum_o") %>% 
                       distinct(enrolid,svcdate,stdplac,dx) %>% 
                       collect(n=Inf) %>% 
                       mutate(enrolid=as.integer64(enrolid))))
  
  
  out_temp <- bind_rows(out_temp1,out_temp2,out_temp3) %>% 
    select(-setting) %>% 
    unnest() %>% 
    group_by(db,year) %>% 
    nest()
  rm(out_temp1,out_temp2,out_temp3)
  #### Merge in the keys ####
  ## Merge inpatient keys
  in_dx_keys <- db.con %>% 
    tbl("inpatient_keys") %>% 
    collect(n=Inf) %>% 
    select(ccae,year,caseid,key) %>% 
    inner_join(in_temp %>%
                 mutate(ccae=ifelse(db=="ccae",1L,0L)) %>% 
                 unnest(),
               by = c("ccae", "year", "caseid")) %>% 
    select(dx,key) 
  
  ## Merge outpatient keys
  out_dx_keys <- db.con %>% 
    tbl("outpatient_keys") %>% 
    collect(n=Inf) %>% 
    select(enrolid,stdplac,svcdate,key) %>% 
    inner_join(out_temp %>% 
                 select(data) %>% 
                 unnest(),
               by = c("enrolid", "stdplac", "svcdate")) %>% 
    select(key,dx) 
  
  dx_keys <- bind_rows(in_dx_keys,out_dx_keys) %>% 
    distinct()
  
  #### Return ####
  return(dx_keys)
}

get_dx_links_primary_hosp <- function(icd_9_codes,icd_10_codes,db.con=db_con,collect_tab=collect_table()){
  
  ##### get inpatient diagnoses ####
  in_temp1 <- collect_tab %>% 
    filter(as.integer(year)<15) %>% 
    mutate(data=map2(db,year,
                     ~tbl(db.con,paste0("inpatient_dx_",.x,"_",.y)) %>% 
                       filter(dx %in% icd_9_codes) %>% 
                       filter(dx_num==1) %>% 
                       distinct(caseid,dx) %>% 
                       inner_join(tbl(db.con,paste0("inpatient_core_",.x,"_",.y)) %>% 
                                    select(caseid,enrolid),
                                  by="caseid") %>% 
                       distinct(enrolid,caseid,dx) %>% 
                       collect(n=Inf) %>% 
                       mutate(enrolid=as.integer64(enrolid))))
  
  in_temp2 <- collect_tab %>% 
    filter(as.integer(year)>14) %>% 
    mutate(data=map2(db,year,
                     ~tbl(db.con,paste0("inpatient_dx9_",.x,"_",.y)) %>% 
                       filter(dx %in% icd_9_codes) %>%
                       filter(dx_num==1) %>% 
                       distinct(caseid,dx) %>%
                       inner_join(tbl(db.con,paste0("inpatient_core_",.x,"_",.y)) %>% 
                                    select(caseid,enrolid),
                                  by="caseid") %>% 
                       distinct(enrolid,caseid,dx) %>% 
                       collect(n=Inf) %>% 
                       mutate(enrolid=as.integer64(enrolid))))
  
  in_temp3 <- collect_tab %>% 
    filter(as.integer(year)>14) %>% 
    mutate(data=map2(db,year,
                     ~tbl(db.con,paste0("inpatient_dx10_",.x,"_",.y)) %>% 
                       filter(dx %in% icd_10_codes) %>% 
                       filter(dx_num==1) %>% 
                       distinct(caseid,dx) %>%
                       inner_join(tbl(db.con,paste0("inpatient_core_",.x,"_",.y)) %>% 
                                    select(caseid,enrolid),
                                  by="caseid") %>% 
                       distinct(enrolid,caseid,dx) %>% 
                       collect(n=Inf) %>% 
                       mutate(enrolid=as.integer64(enrolid))))
  
  in_temp <- bind_rows(in_temp1,in_temp2,in_temp3) %>% 
    select(-setting) %>% 
    unnest() %>% 
    group_by(db,year) %>% 
    nest()
  
  rm(in_temp1,in_temp2,in_temp3)
  
  #### Merge in the keys ####
  ## Merge inpatient keys
  in_dx_keys <- db.con %>% 
    tbl("inpatient_keys") %>% 
    collect(n=Inf) %>% 
    select(ccae,year,caseid,key) %>% 
    inner_join(in_temp %>%
                 mutate(ccae=ifelse(db=="ccae",1L,0L)) %>% 
                 unnest(),
               by = c("ccae", "year", "caseid")) %>% 
    select(dx,key) 
  
  ## Merge outpatient keys
  
  dx_keys <- in_dx_keys %>% 
    distinct()
  
  #### Return ####
  return(dx_keys)
}


get_proc_links <- function(icd_9_codes,icd_10_codes,collect_tab=collect_table()) {
  
}

get_rx <- function(ndc_codes,db.con=db_con,collect_tab=collect_table()){
  collect_tab %>% 
    select(-setting) %>% 
    mutate(data=map2(db,year,~tbl(db.con,paste0("rx_core_",.x,"_",.y)) %>% 
                       filter(ndcnum %in% ndc_codes) %>% 
                       select(enrolid,ndcnum,svcdate) %>% 
                       collect() %>% 
                       mutate(enrolid=as.integer64(enrolid))))
}
