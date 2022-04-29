#!/usr/bin/Rscript
# this program calculates outcomes from interpolated datasets
# and saves them as an .rds file

setwd("/gpfs/fs0/scratch/g/gtomlins/cyarnell/ObservedCriteria")
library(dplyr)
library(tidyr)
library(ggplot2)
source("/gpfs/fs0/scratch/g/gtomlins/cyarnell/ObservedCriteria/models/check_thresholds_obsonly.R")

# colours
c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

c_blue = "cornflower blue"

logit <- function(x){log(x/(1-x))}
inv_logit <- function(x){exp(x)/(1+exp(x))}

# take in the patient id here

# load the interpolation

dat <- readRDS("./data/interpolation_fixedGP_dat.rds")

pdx_b <- readRDS("data/oc_pdx_b.rds")

sample = 100

# completed patients
completed <- list.files("outcomes2/")
completed <- data.frame(filename = completed) %>% 
  mutate(stay_id = as.numeric(gsub(".rds","", gsub("outcomes_","",filename))))

todo <- dat$patients

outcomes_only <- function(pt_stay_id){

interp <- readRDS(paste0("interpolations2/interpolations_",pt_stay_id,".rds"))

df <- interp$df %>%
    # no rules use sbp
    filter(variable != "sbp") %>%
    select(stay_id, iter, time, variable, value)
pt_bl <- pdx_b %>% filter(stay_id == pt_stay_id)

# boil it down to the outcomes

# number of thresholds
Nthresh <- 17

outcomes <- data.frame(iter = rep(1:sample, each = Nthresh),
                       threshold = NA,
                       met = NA,
                       followed = NA,
                       fol8 = NA,
                       fol24 = NA,
                       fol48 = NA,                       
                       fol72 = NA,
                       fol120 = NA,
                       death3 = NA,
                       niv3 = NA,
                       hfnc3 = NA,
                       niv = NA,
                       hfnc = NA,
                       wob = NA,
                       hour = NA)

imv <- pt_bl$imv_ever
death <- pt_bl$death_obs

pos = 0
for (i in 1:sample){
  df_temp <- filter(df, iter == i)  
  outcomes[(pos+1):(pos+Nthresh),-1] <- check_thresholds(df_temp, imv, death)
  pos = pos + Nthresh
}

output_list <- list(stay_id = pt_stay_id,
                    pt_bl = pt_bl,
                    outcomes = outcomes)

saveRDS(output_list,
        file = paste0("outcomes2/outcomes_",
                      pt_stay_id,
                      ".rds"))

print(paste0("finished patient ", pt_stay_id))
}

library(doParallel)

print("start parallelization")

ncores = Sys.getenv("SLURM_CPUS_PER_TASK") 
registerDoParallel(cores=80)

output <- foreach(i = 1:nrow(todo)) %dopar% outcomes_only(todo$stay_id[i])