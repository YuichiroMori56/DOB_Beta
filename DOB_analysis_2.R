library(tidyverse)
library(data.table)
library(bit64)
library(lubridate)
library(tableone)
library(survival)
library(modelr)
library(survminer)
library(missRanger)
library(doParallel)
library(MatchIt)
library(cobalt)
library(stddiff)

rm(list = ls())

# Read dataset for analysis -----------------------------------------------
dataset.for.analysis = fread('**********/dataset.for.analysis.csv') %>% dplyr::select(-V1) %>% 
  mutate(across(all_of(c("sex", "hpscalecode", "smokingexp", "hfbloodpressure", "Diure.firstday", "Vasodil.firstday", "Aline.firstday", "outcome", "destination")),  as.factor),
         across(all_of(c("ID", "AST", "CK", "LD", "WBC", "ALT", "ggt")), as.double),
         sex = factor(as.character(sex), levels = c("male", "female"))) 

# Define predictor variables ----------------------------------------------
RFimpVariables = c("age","sex","bmi","Hb", "hpscalecode",
                   "eGFR","AST","NTproBNP","Na", "K","Alb", "smokingexp",
                   "hfbloodpressure", "CRP",
                   "CK", "LD", "BUN", "Plt", "WBC", "ALT", "ggt", "bil",
                   "Diure.firstday", "Vasodil.firstday", "Aline.firstday")

PredictorVariables = c("age","sex","bmi", "Hb","Na", "K", "Alb", "hfbloodpressure",
                       "eGFR","NTproBNP_log", "BUN_log")

# RFimp --------------------------------------------------------------
cores <- detectCores(logical = FALSE)
registerDoParallel(cores = cores)
RFimp <- function(data){
  imp.df <- missRanger(data, . -ID - patientid - outcome - destination - firstday_beta - beta - death - fu_end - beta.initiated.with.dob
                       ~ . -ID - patientid - outcome - destination - firstday_beta - beta - death - fu_end - beta.initiated.with.dob,
                       maxiter = 10,
                       num.trees = 100, verbose = 0)
  return(imp.df)
}

# Clone and censor----------------------------------------------

#strategy = 0: Before day n, do not initiate beta with dob
#strategy = 1: Initiate beta with dob before day n, 

clone.censor = function(data, graceperiod, obs.time, discharged.alive=FALSE){
  nday = graceperiod
  t = data.frame(strategy = rep(0:1, each = nrow(data)),  ID = rep(data$ID, 2))
  dataset_cloned = data %>% 
    left_join(t, by = "ID") %>%
    #delta.w: administrative censoring (0: no, 1: yes)
    #time.w: the date of delta.w determined
    #delta.o: outcome indicator (0: alive, 1: death)
    #time.o: followup time in the emulated trial
    mutate(delta.w = case_when(strategy == 1 & beta == 0 ~ ifelse(fu_end <= nday, 0, 1),
                               strategy == 1 & beta == 1 ~ ifelse(beta.initiated.with.dob == 1 & firstday_beta <= nday, 0,
                                                                  ifelse(fu_end <= nday, 0, 1)),
                               strategy == 0 & beta == 0 ~ 0,
                               strategy == 0 & beta == 1 ~ ifelse(beta.initiated.with.dob == 1 & firstday_beta <= nday, 1, 0)),
           time.w = case_when(strategy == 1 & beta == 0 ~ ifelse(fu_end <= nday, fu_end, nday),
                              strategy == 1 & beta == 1 ~ ifelse(beta.initiated.with.dob == 1 & firstday_beta <= nday, 
                                                                 firstday_beta, 
                                                                 ifelse(fu_end <= nday, fu_end, nday)),
                              strategy == 0 & beta == 0 ~ ifelse(fu_end <= nday, fu_end, nday),
                              strategy == 0 & beta == 1 ~ ifelse(beta.initiated.with.dob == 1 & firstday_beta <= nday, 
                                                                 firstday_beta,
                                                                 ifelse(fu_end <= nday, fu_end, nday))),
           delta.o = ifelse(delta.w == 0, death, 0),
           time.o = ifelse(delta.w == 0, fu_end, time.w),
           #when admission duration > obs.time
           delta.w = ifelse(time.o > obs.time, 0, delta.w),
           time.w = ifelse(time.o > obs.time, obs.time,time.w),
           delta.o = ifelse(time.o > obs.time, 0, delta.o),
           time.o = ifelse(time.o > obs.time, obs.time,time.o)
    ) 
  if(discharged.alive==TRUE){
    #for sensitivity analysis#
    #treat discharged patients alive until day 60 if there is a record of "cured" at discharge
    dataset_cloned = dataset_cloned %>%
      mutate(time.o = ifelse(time.o < 60 & outcome %in% c("cure", "recovery") &
                          destination == "discharge" & delta.w == 0, 60,time.o))
  }
  return(dataset_cloned)
}

# IPCW --------------------------------------------------------
weight.d = function(data){
  d.weighted = data.table()
  for (i in 0:1) { #conduct IPCW in each strategy
    d2 = data %>%
      dplyr::filter(strategy == i) %>%
      mutate(Tstart = 0)
    t_events = sort(unique(d2$time.o))
    times = data.frame("tevent"=t_events, "ID_t"=seq(1:length(t_events)))
    d.long.cens = survSplit(d2, cut = t_events, end = "time.o", start = "Tstart", event = "delta.w", id = "id") %>%
      arrange(id, time.o)
    d.long.event = survSplit(d2, cut = t_events, end = "time.o", start = "Tstart", event = "delta.o", id = "id") %>%
      arrange(id, time.o) %>%
      dplyr::select(-delta.w) %>%
      left_join(d.long.cens %>% dplyr::select(delta.w, id, Tstart), by = c("id","Tstart")) %>%
      mutate(Tstop = time.o) %>%
      left_join(times, by = c("Tstart" = "tevent")) %>%
      arrange(id, time.o) %>%
      replace_na(replace = list(ID_t = 0)) 
    ms_cens = coxph(formula(paste("Surv(Tstart, Tstop, delta.w) ~",
                                  paste(PredictorVariables, collapse = "+"))), data = d.long.event)
    cens  = coxph(Surv(Tstart, Tstop, delta.w) ~ 1, data = d.long.event)
    d.long.event2 = d.long.event %>% 
      add_predictions(ms_cens, var = "lin_pred") %>%
      left_join(basehaz(ms_cens,centered=F) %>%
                  data.frame() %>% 
                  left_join(times, by = c("time" = "tevent")), by = "ID_t") %>%
      left_join(basehaz(cens) %>%
                  data.frame() %>%
                  rename(hazard.unadj = hazard) %>%  
                  left_join(times, by = c("time" = "tevent")), by = c("ID_t", "time")) %>%
      arrange(id, time.o) %>%
      replace_na(replace = list(hazard = 0, hazard.unadj = 0)) %>%
      mutate(P_uncens = exp(-(hazard)*exp(lin_pred)),
             P_uncens.unadj = exp(-(hazard.unadj)),
             s.weight_Cox =  P_uncens.unadj/P_uncens,
             s.weight_trunc = pmin(pmax(s.weight_Cox, quantile(s.weight_Cox, .025)), 
                                        quantile(s.weight_Cox, .975)))
    d.weighted = rbind(d.weighted, d.long.event2)
  }
  return(d.weighted)
}


#Point estimate for each outcome indicator
analysis = function(data, cutoffday = 30){
  d.w = data
  if (max(d.w$s.weight_trunc) >= Inf){
    res = c(NA,NA)
    names(res)<-c("Diff_surv.early.minus.cons","Diff_RMST.early.minus.cons")
  }
  else{
    emul_Cox <- survfit(Surv(Tstart, Tstop, delta.o) ~ strategy,
                        data= d.w, weights = iptw.stab_trunc * s.weight_trunc)
    Diff_surv.early.minus.cons = min(emul_Cox[strata = "strategy=1"]$surv) - min(emul_Cox[strata = "strategy=0"]$surv)
    fit.table.dob <- summary(emul_Cox[strata = "strategy=0"], rmean=cutoffday)$table
    fit.table.beta <- summary(emul_Cox[strata = "strategy=1"], rmean=cutoffday)$table
    Diff_RMST.early.minus.cons = fit.table.beta["rmean"] - fit.table.dob["rmean"]
    # res = c(Diff_surv.early.minus.cons,Diff_RMST.early.minus.cons)
    # names(res)<-c("Diff_surv.early.minus.cons","Diff_RMST.early.minus.cons") 
    res = c(1-min(emul_Cox[strata = "strategy=1"]$surv), 1-min(emul_Cox[strata = "strategy=0"]$surv), fit.table.beta["rmean"],fit.table.dob["rmean"])
    names(res)<-c("Mort-early","Mort-cons","RMST-early","RMST-cons")
  }
  return(res)
}


# Analysis ----------------------------------------------------------------
set.seed(1234)
Impdf = RFimp(dataset.for.analysis) %>% mutate(ID = row_number()) %>% 
  mutate(beta.initiated.with.dob = ifelse(beta.initiated.with.dob == 1 & firstday_beta <= 7, 1,0)) 

PSM_formula = formula(paste("beta.initiated.with.dob ~",
                            paste(PredictorVariables, collapse = "+")))
m.out = matchit(PSM_formula, data = Impdf, method = "nearest", caliper = 0.10)
summary(m.out, standardize = T)
balance =  summary(m.out, standardize = T)
cbind(balance$sum.all[,1:3], balance$sum.matched[,1:3])
balance_df = cbind(balance$sum.all[,1:3], balance$sum.matched[,1:3]) %>% 
  as.data.frame() %>% 
  round(3)
write.csv(balance_df, file = "*********/balance.csv", row.names = T)
love.plot(m.out, 
          threshold = 0.1, 
          abs = T, 
          drop.distance = T, 
          shapes = c(1, 16), 
          sample.names = c("Unmatched","Matched"),
          var.order = "unadjusted")
plot(m.out, type = "jitter", interactive = F)
plot(m.out, type = "hist")      

newdf = match.data(m.out) %>% 
  mutate(m = mean(Impdf$beta.initiated.with.dob)) %>% 
  group_by(beta.initiated.with.dob) %>% 
  mutate(iptw = ifelse(beta.initiated.with.dob == 1, 1/distance, 1/(1-distance)),
         iptw.stab = ifelse(beta.initiated.with.dob == 1, m/distance, (1-m)/(1-distance)),
         iptw.stab_trunc = pmin(pmax(iptw.stab, quantile(iptw.stab, .025)), 
                                quantile(iptw.stab, .975)),
         subclass = as.integer(subclass)) %>% 
  ungroup(beta.initiated.with.dob)
newdf.clone.censor = clone.censor(newdf, graceperiod = 7, obs.time = 30) 
newdf.clone.censor.w = weight.d(newdf.clone.censor)

analysis(newdf.clone.censor.w)

RF.emul_Cox <- survfit(Surv(Tstart, Tstop, delta.o) ~ strategy,
                       data= newdf.clone.censor.w, weights = iptw.stab_trunc * s.weight_trunc)

ggsurvplot(RF.emul_Cox,
           xlab = "Days from admission",
           size = 1.2, 
           linetype =  c("twodash", "solid"),
           risk.table = TRUE,
           legend.title = "Strategy",
           legend.labs = c("Conservative","Early initiation")) %>% 
  ggpar(font.x = c(25,"bold"),
        font.y = c(25,"bold"),
        font.tickslab = c(25,"bold"))


boot = function(d)
{
  return(d[sample(1:dim(d)[1], replace=TRUE), ])
}
set.seed(1234)
system.time(
  RD <- replicate(1000, boot(newdf %>% distinct(subclass)) %>%
                    left_join(newdf, by ="subclass") %>% 
                    clone.censor(graceperiod = 7, obs.time = 30) %>% 
                    weight.d() %>% 
                    analysis())
)
result.main = apply(RD, 1, quantile, c(0.025, 0.5,  0.975), na.rm = TRUE)
result.main 

Cox = coxph(Surv(Tstart, Tstop, delta.o) ~ strategy,
            data= newdf.clone.censor.w, weights = iptw.stab_trunc * s.weight_trunc)
summary(Cox)
