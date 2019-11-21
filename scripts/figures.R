rm(list=ls()) 
library(bit64)
library(tidyverse)
library(icd)
library(survival)
library(icdplus)
source("functions/small_db_functions.R")
source("functions/analysis_functions.R")

# load in results data from main_analysis.R
load("<path to file>")

carrier_or <- res_carrier
cf_or <- res_cf

db_con_carrier <- src_sqlite("<base_path>/carrier.db")
db_con_cf <- src_sqlite("<base_path>/carrier.db")

##### Build Timemap ####
carrier_study_cohort <- tbl(db_con_carrier,"carrier_study_cohort") %>% collect(n=Inf)
cf_study_cohort <- tbl(db_con_cf,"cf_study_cohort") %>% collect(n=Inf)


gen_all_cond_plot_log_new <- function(data,study_cohort_data){
  
  tab_counts <- study_cohort_data %>% 
    summarise(tot_case=sum(case==1),
              tot_control=sum(case==0),
              tot_male_case=sum(case==1 & sex==1),
              tot_male_control=sum(case==0 & sex==1))
  
  main_vals <- data %>% 
    arrange(group,label) %>% 
    mutate(n_case=ifelse(name!="infert_male",
                         100*(n_case/tab_counts$tot_case),
                         100*(n_case/tab_counts$tot_male_case)),
           n_control=ifelse(name!="infert_male",
                            100*(n_control/tab_counts$tot_control),
                            100*(n_control/tab_counts$tot_male_control))) %>% 
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
    mutate_at(vars(n_case,n_control),funs(round(.,3))) %>%
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
  
  text(-3.7,nrow(plot_data)+3.25, "Prevalence Rate (%)", cex=1.2, font=2, pos=4) 
  text(-4.3,nrow(plot_data)+2, "CF Carriers", cex=1.2, font=2, pos=4)
  text(-4.15,nrow(plot_data)+1, paste0("N = ",formatC(tab_counts$tot_case,big.mark=",")), cex=1.1, font=1, pos=4)
  t2 <- ifelse(!is.na(plot_data$n_case), format(plot_data$n_case,big.mark=","), '')
  text(-3.9, rowseq, t2, cex=1.1, pos=4)
  
  text(-2.6,nrow(plot_data)+2, "Matched Cohort", cex=1.2, font=2, pos=4)
  text(-2.25,nrow(plot_data)+1, paste0("N = ",formatC(tab_counts$tot_control,big.mark=",")), cex=1.1, font=1, pos=4)
  t2 <- ifelse(!is.na(plot_data$n_control), format(plot_data$n_control,big.mark=","), '')
  text(-2.1, rowseq, t2, cex=1.1, pos=4)
  
  #text(1,nrow(plot_data)+2.5, expression("log"[e]*" Odds Ratio"), cex=1.2, font=2, pos=4)
  text(.1,nrow(plot_data)+2.5, "Natural Log Odds Ratio", cex=1.2, font=2, pos=4)
  
  text(3.8,nrow(plot_data)+2.5, "Odds Ratio (95% CI)", cex=1.2, font=2, pos=4)
  t3 <- ifelse(!is.na(plot_data$or),plot_data$or_text,'')
  text(4.1,rowseq, t3, cex=1.1, pos=4)
  
  text(6.6,nrow(plot_data)+2.5, "P-Value", cex=1.2, font=2, pos=4)
  #text(6.3,nrow(plot_data)+1.5, "(One-Sided)", cex=1.2, font=2, pos=4)
  t2 <- ifelse(!is.na(plot_data$p_val), format(plot_data$p_val,big.mark=","), '')
  text(7.6, rowseq, t2, cex=1.1, pos=2)
} 


pdf("~/Desktop/temp_figs/new_carrier_table.pdf",height = 16,width = 14,family = "Times")
gen_all_cond_plot_log_new(data=carrier_or %>% 
                        mutate(or_p_val=or_p_val/2),
                      study_cohort_data = carrier_study_cohort)
text(-9,73+5.5, "Figure 2", cex=1.4, font=2, pos=4)
dev.off()


#### Correlation plot #####
load("../data/temp_results/cf_or.RData")
cftr_any_or <- carrier_or
cftr_any_or <- cftr_any_or %>% 
  mutate(label=rem_z(label)) %>%
  mutate(label=capitalize(as.character(label))) %>%
  mutate(label=ifelse(label=="Cancer of colon, stomach, gi organs",
                      "Cancer of colon, stomach, GI organs",
                      ifelse(label=="Gerd",
                             "GERD",
                             ifelse(label=="Diabetes (Type 1 or Secondary)",
                                    "Diabetes (Type 1 or secondary)",
                                    label))))
cf_any_or <- cf_or
cf_any_or <- cf_any_or %>% 
  mutate(label=rem_z(label)) %>%
  mutate(label=capitalize(as.character(label))) %>%
  mutate(label=ifelse(label=="Cancer of colon, stomach, gi organs",
                      "Cancer of colon, stomach, GI organs",
                      ifelse(label=="Gerd",
                             "GERD",
                             ifelse(label=="Diabetes (Type 1 or Secondary)",
                                    "Diabetes (Type 1 or secondary)",
                                    label))))


corr_conds <- cftr_any_or %>%
  .$name


p1 <- plot_main_cor(corr_conds) +
  geom_hline(yintercept = 0,linetype=2) +
  theme(axis.text.y = element_text(size = 15,color="black"),
        axis.text.x = element_text(size = 16,color="black"),
        axis.title = element_text(size=16,color="black")) +
  theme(axis.line = element_line(color="black"))


plot_scatter_new <- function(conditions){
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
    xlab("CF Cases") +
    ylab("CF\nCarriers") +
    ggthemes::theme_few() +
    geom_abline(slope = 1, intercept = 0,linetype=2) +
    theme(text=element_text(family = "Times")) +
    scale_x_continuous(breaks = c(log(1),log(7),log(50),log(400),log(3000)),
                       labels=c(1,7,50,400,3000)) +
    scale_y_continuous(breaks = c(log(1),log(4),log(12)),
                       labels=c(1,4,12), limits = c(0,2.8)) +
    theme(axis.title.y = element_text(angle = 0, vjust = 1.2,
                                      hjust = 1))
}

p2 <- plot_scatter_new(corr_conds) +
  ggtitle("Odds Ratios (ln scale axis)") +
  theme(plot.title = element_text(size = 18,hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 16,color="black"),
        axis.text.y = element_text(size = 16,color="black"),
        axis.title = element_text(size=16,color="black"))  +
  theme(axis.line = element_line(color="black"))

p2

#### Plot combined correlation plot
g <- ggplotGrob(p2)
pdf("~/Desktop/temp_figs/new_corr_fig2_1.pdf",height=16,width = 14)
p1 +
  theme(legend.text=element_text(size=16),
        legend.position = c(.85,.25),
        legend.background = element_rect(linetype = "solid")) +
  annotation_custom(
    grob = g,
    xmin = 33,
    xmax = 60,
    ymin = 2.1,
    ymax = 9.3
  ) +
  theme(plot.title = element_text(face="bold",size=18,hjust = -.38,vjust = 1))
dev.off()


g <- ggplotGrob(p2)
pdf("~/Desktop/temp_figs/new_corr_fig2_2.pdf",height=16,width = 14)
p1 +
  theme(legend.text=element_text(size=16),
        legend.background = element_rect(linetype = "solid")) +
  annotation_custom(
    grob = g,
    xmin = 33,
    xmax = 60,
    ymin = 2.1,
    ymax = 9.3
  ) +
  theme(plot.title = element_text(face="bold",size=18,hjust = -.38,vjust = 1))
dev.off()
  
p1 + theme(legend.position = c(.8,.3))
