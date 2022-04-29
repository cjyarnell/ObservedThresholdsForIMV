# check each interpolation against the proposed thresholds

check_thresholds <- function(df, imv, death){

# pivot makes it easier
  
df <- df %>%
  spread(key = "variable",
         value = "value")

if(!("ph" %in% names(df))){df$ph = NA}
if(!("wob" %in% names(df))){df$wob = NA}
if(!("po2" %in% names(df))){df$po2 = NA}
if(!("pco2" %in% names(df))){df$pco2 = NA}
if(!("GCS" %in% names(df))){df$GCS = NA}
if(!("niv" %in% names(df))){df$niv = NA}
if(!("stdo2" %in% names(df))){df$stdo2 = NA}
if(!("hfnc" %in% names(df))){df$hfnc = NA}
    
endtime <- max(df$time)
  
# thresholds

thresholds <- c("PF150",
                "PF100",
                "PF80",
                "SF180",
                     "SF120",
                     "SF90",
                     "ROX5",
                     "ROX4",
                     "ROX3",
                     "HACOR5",
                     "HACOR10",
                     "FLORALI",
                     "DeMontBauer",
                     "GCS12",
                     "GCS9",
                     "Pressor",
                     "Darreau")

Nthresh = length(thresholds)

threshold_names <- thresholds

outcomes <- data.frame(threshold = thresholds,
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

# functions for thresholds

pf <- function(df){100*df$po2/df$fio2}

sf <- function(df){df$spo2*100/df$fio2}

rox <- function(df){df$spo2*100/df$fio2/df$resp_rate}

  hacor <- function(df){
    hr_points <- ifelse(df$heart_rate > 120,1,0)
    ph_points <- ifelse(is.na(df$ph), 0, ifelse((df$ph >= 7.35) %in% TRUE, 0,
                        ifelse((df$ph >= 7.30) %in% TRUE, 2,
                               ifelse((df$ph >= 7.25) %in% TRUE, 3, 4))))
    gcs_points <- ifelse(is.na(df$GCS), 0, ifelse((df$GCS == 15) %in% TRUE, 0,
                         ifelse((df$GCS >= 13)  %in% TRUE, 2,
                                ifelse((df$GCS >= 11) %in% TRUE, 5, 10))))
    pf_points <- ifelse(is.na(pf(df)), 0, ifelse((pf(df) > 200) %in% TRUE, 0,
                        ifelse((pf(df) > 175) %in% TRUE, 2,
                               ifelse((pf(df) > 150) %in% TRUE, 3,
                                      ifelse((pf(df) > 125) %in% TRUE, 4,
                                             ifelse((pf(df) > 100) %in% TRUE, 5, 6))))))
    rr_points <- ifelse(is.na(df$resp_rate), 0, ifelse((df$resp_rate <= 30) %in% TRUE, 0,
                        ifelse((df$resp_rate <= 35) %in% TRUE, 1,
                               ifelse((df$resp_rate <= 40) %in% TRUE, 2, 
                                      ifelse((df$resp_rate <= 45) %in% TRUE, 3, 4)))))
    
    hacor_points <- hr_points + ph_points + gcs_points + pf_points + rr_points
    hacor_points
  }
    
florali <- function(df){
  f_sat <- (df$spo2 < 90) & (df$fio2 >= 90)
  f_rr  <- df$resp_rate > 40
  f_wob <- df$wob == 1
  f_ph  <- (df$ph < 7.35)  %in% TRUE
  f_gcs <- df$GCS < 12
  f_prs <- (df$pressor == 1) 
  
  f_resp <- (f_sat + f_rr + f_wob + f_ph)>=2
  f_all <- f_resp | f_gcs | f_prs
  f_all
}

demontb <- function(df){
  dmb_sf <- sf(df) < 120
  dmb_rr <- df$resp_rate > 30
  dmb_wob<- df$wob == 1
  dmb_all <- dmb_sf & dmb_rr & dmb_wob
  dmb_all
}

darreau <- function(df){
  f_sat <- (df$spo2 < 90) & (df$fio2 >= 90)
  f_rr  <- df$resp_rate > 35
  f_wob <- df$wob == 1
  f_ph  <- ((df$ph < 7.35) & (df$pco2 >45)) %in% TRUE
  f_gcs <- df$GCS < 10
  f_prs <- (df$pressor == 1)
  
  f_resp <- (f_sat + f_rr + f_wob + f_ph)>=2
  f_all <- (f_resp | f_gcs) & f_prs
  f_all
}
    

thresh_followed <- function(temp, endtime, imv){
  c(ifelse(temp[1] == 0, NA,
         ((endtime - temp[2]) < 3) & (imv == 1)),
    ifelse(temp[1] == 0, NA,
         ((endtime - temp[2]) < 8) & (imv == 1)),
    ifelse(temp[1] == 0, NA,
         ((endtime - temp[2]) < 24) & (imv == 1)),
    ifelse(temp[1] == 0, NA,
         ((endtime - temp[2]) < 48) & (imv == 1)),
    ifelse(temp[1] == 0, NA,
         ((endtime - temp[2]) < 72) & (imv == 1)),
    ifelse(temp[1] == 0, NA,
         ((endtime - temp[2]) < 120) & (imv == 1))
   )
}

check_thresh <- function(binvec, df){
    if(max(binvec, na.rm = T) != 1){
        return(c(met = 0, time = NA, row = NA))
    } else {
        return(c(
        met = max(binvec, na.rm = T),
        time = df$time[which.max(binvec)],
        row = which.max(binvec)))
    }
}
    
calc_outcomes <- function(temp, endtime, imv, death, df){
    met = temp[1]
    if (met == 0){return(c(0,rep(NA, 13)))} else {
    followed <- thresh_followed(temp, endtime, imv)
    death3 <- ifelse((endtime-temp[2]) < 3 & death == 1,1,0)
    
    if(endtime == temp[2]){
        niv3 = 0
        hfnc3 = 0
    } else {
        tmpdf <- filter(df, time >= temp[2], time < (temp[2]+3))
        niv3 <-  ifelse(df$niv[temp[3]] != 1 & max(tmpdf$niv == 1), 1, 0)
        hfnc3 <- ifelse(df$hfnc[temp[3]] != 1 & max(tmpdf$hfnc == 1), 1, 0)
    }
    tmpniv <- df$niv[temp[3]]
    tmpstdo2 <- df$hfnc[temp[3]]
    tmpwob <- df$wob[temp[3]]
    c(met, followed, death3, niv3, hfnc3, tmpniv, tmpstdo2, tmpwob, temp[2])
    }
}
      
  # pf 
  temp <- check_thresh(((pf(df) < 150) %in% TRUE 
                        & (df$niv == 1 | df$hfnc == 1 | df$stdo2 == 1)),
                      df)
  outcomes[1,2:15] <- calc_outcomes(temp, endtime, imv, death, df)
  
  temp <- check_thresh(((pf(df) < 100) %in% TRUE 
                        & (df$niv == 1 | df$hfnc == 1 | df$stdo2 == 1)),
                      df)
  outcomes[2,2:15] <- calc_outcomes(temp, endtime, imv, death, df)

    
  temp <- check_thresh(((pf(df) < 80) %in% TRUE 
                        & (df$niv == 1 | df$hfnc == 1 | df$stdo2 == 1)),
                      df)
  outcomes[3,2:15] <- calc_outcomes(temp, endtime, imv, death, df)

  # sf
  temp <- check_thresh(((sf(df) < 180) %in% TRUE 
                        & (df$niv == 1 | df$hfnc == 1 | df$stdo2 == 1)),
                      df)
  outcomes[4,2:15] <- calc_outcomes(temp, endtime, imv, death, df)

  temp <- check_thresh(((sf(df) < 120) %in% TRUE 
                        & (df$niv == 1 | df$hfnc == 1 | df$stdo2 == 1)),
                      df)
  outcomes[5,2:15] <- calc_outcomes(temp, endtime, imv, death, df)

  temp <- check_thresh(((sf(df) < 90) %in% TRUE 
                        & (df$niv == 1 | df$hfnc == 1 | df$stdo2 == 1)),
                      df)
  outcomes[6,2:15] <- calc_outcomes(temp, endtime, imv, death, df)

  # rox
  temp <- check_thresh(((rox(df) < 5) %in% TRUE 
                        & (df$niv == 1 | df$hfnc == 1 | df$stdo2 == 1)),
                      df)
  outcomes[7,2:15] <- calc_outcomes(temp, endtime, imv, death, df)

  temp <- check_thresh(((rox(df) < 4) %in% TRUE 
                        & (df$niv == 1 | df$hfnc == 1 | df$stdo2 == 1)),
                      df)
  outcomes[8,2:15] <- calc_outcomes(temp, endtime, imv, death, df)

  temp <- check_thresh(((rox(df) < 3) %in% TRUE 
                        & (df$niv == 1 | df$hfnc == 1 | df$stdo2 == 1)),
                      df)
  outcomes[9,2:15] <- calc_outcomes(temp, endtime, imv, death, df)

  # hacor
  temp <- check_thresh(((hacor(df) > 5) %in% TRUE 
                        & (df$niv == 1 | df$hfnc == 1 | df$stdo2 == 1)),
                      df)
  outcomes[10,2:15] <- calc_outcomes(temp, endtime, imv, death, df)

  temp <- check_thresh(((hacor(df) > 10) %in% TRUE 
                        & (df$niv == 1 | df$hfnc == 1 | df$stdo2 == 1)),
                      df)
  outcomes[11,2:15] <- calc_outcomes(temp, endtime, imv, death, df)
  
  # florali
  temp <- check_thresh((florali(df) %in% TRUE 
                        & (df$niv == 1 | df$hfnc == 1 | df$stdo2 == 1)),
                      df)
  outcomes[12,2:15] <- calc_outcomes(temp, endtime, imv, death, df)
  
  # deMontBauer
  temp <- check_thresh((demontb(df) %in% TRUE 
                        & (df$niv == 1 | df$hfnc == 1 | df$stdo2 == 1)),
                      df)
  outcomes[13,2:15] <- calc_outcomes(temp, endtime, imv, death, df)
  
  # GCS
  temp <- check_thresh(((df$GCS < 12) %in% TRUE 
                        & (df$niv == 1 | df$hfnc == 1 | df$stdo2 == 1)),
                      df)
  outcomes[14,2:15] <- calc_outcomes(temp, endtime, imv, death, df)

  temp <- check_thresh(((df$GCS < 9) %in% TRUE 
                        & (df$niv == 1 | df$hfnc == 1 | df$stdo2 == 1)),
                      df)
  outcomes[15,2:15] <- calc_outcomes(temp, endtime, imv, death, df)

  # pressor
  temp <- check_thresh(((df$pressor == 1) %in% TRUE 
                        & (df$niv == 1 | df$hfnc == 1 | df$stdo2 == 1)),
                      df)
  outcomes[16,2:15] <- calc_outcomes(temp, endtime, imv, death, df)
 
 # darreau
  temp <- check_thresh((darreau(df) %in% TRUE 
                        & (df$niv == 1 | df$hfnc == 1 | df$stdo2 == 1)),
                      df)
  outcomes[17,2:15] <- calc_outcomes(temp, endtime, imv, death, df)
    
outcomes

}
