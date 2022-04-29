#!/usr/bin/Rscript
# this program calculates outcomes from interpolated datasets
# and saves them as an .rds file
# for AMDS cohort

setwd("/gpfs/fs0/scratch/g/gtomlins/cyarnell/ObservedCriteria")
library(dplyr)
library(tidyr)
library(ggplot2)
source("/gpfs/fs0/scratch/g/gtomlins/cyarnell/ObservedCriteria/models/amd_check_thresholds.R")

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

dat <- readRDS("./data/amd_interpolation_dat.rds")

amd_b <- readRDS("./data/oc_amd_b.rds")

sample = 100


todo <- dat$patients

outcomes_only <- function(pt_stay_id){

interp <- readRDS(paste0("amd_interpolations/amd_interpolations_",pt_stay_id,".rds"))

df <- interp$df %>%
    select(stay_id, iter, time, var.y, value)
pt_bl <- amd_b %>% filter(admissionid == pt_stay_id)

# boil it down to the outcomes

# number of thresholds
Nthresh <- 16

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
                       niv = NA,
                       hour = NA)

imv <- !is.na(pt_bl$imv_time)
    # some date of death are before eligibletime
    # because date of death is a day-resolution number
    # so if date of death is positive 
    # and dateofdeath < 7200 min after enrolment 
    # we take the last observation to be the moment of death
death <- ((pt_bl$dateofdeath < 7200) %in% TRUE) & (is.na(pt_bl$imv_time))

pos = 0
for (i in 1:sample){
  df_temp <- filter(df, iter == i)  
  outcomes[(pos+1):(pos+Nthresh),-1] <- check_thresholds(df_temp, imv, pt_bl$imv_time/60, death)
  pos = pos + Nthresh
}

output_list <- list(admissionid = pt_stay_id,
                    pt_bl = pt_bl,
                    outcomes = outcomes)

saveRDS(output_list,
        file = paste0("amd_outcomes/amd_outcomes_",
                      pt_stay_id,
                      ".rds"))

print(paste0("finished patient ", pt_stay_id))
}

library(doParallel)

print("start parallelization")

ncores = Sys.getenv("SLURM_CPUS_PER_TASK") 
registerDoParallel(cores=80)

output <- foreach(i = 1:nrow(todo)) %dopar% outcomes_only(todo$admissionid[i])

files <- list.files("amd_outcomes/")
stem = "amd_outcomes/"
filenames <- sapply(stem, paste0, files)

df_all <- lapply(filenames, readRDS)

df_id <- list()
df_outcomes <- list()
df_bl <- readRDS("data/oc_amd_b.rds")

for (n in 1:length(df_all)){
    df_id[[n]] <- df_all[[n]]$admissionid
    df_outcomes[[n]] <- df_all[[n]]$outcomes
}

df_id <- unlist(df_id)
df_outcomes <- bind_rows(df_outcomes)

rows_per_patient <- nrow(df_outcomes)/length(df_id)

df_outcomes$id <- rep(df_id, each = rows_per_patient)
df_outcomes <- df_outcomes[,c(14,1:13)]

saveRDS(df_outcomes, "amd_obscrit_outcome_df_all.rds")