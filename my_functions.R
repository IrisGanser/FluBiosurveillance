epi_start <- function(df, cutoff){
  for(j in 2:13){ # first few rows without looking back for previous outbreak
    if(df$bcp.postprob[j] >= cutoff & df$bcp.postprob[j-1] < cutoff &
       df$count_smooth[j] > df$count_smooth[j-1]){
      bcp_start[j] <- df$date[j]
    } else {
      bcp_start[j] <- NA
    }
  }
  for(j in 10:nrow(df)){
    if(df$bcp.postprob[j] >= cutoff & df$bcp.postprob[(j-1)] < cutoff & # transition from non-epidemic to epidemic
       df$count_smooth[j] > df$count_smooth[j-1] & # running mean to ensure that curve is rising (beware of spikes!)
       sum(!is.na(bcp_start[(j-10):(j-1)])) == 0){ # No outbreak flagged during the previous 10 weeks
      bcp_start[j] <- df$date[j]
    } else {
      bcp_start[j] <- NA
    }
  }
  return(bcp_start)
}

epi_end <- function(df, cutoff){ 
  for(j in (nrow(df)-2):2){# run in reverse because otherwise, third criterion cannot be recognized
    if(df$bcp.postprob[j-1] >= cutoff & df$bcp.postprob[j] < cutoff & # transition from non-epidemic to epidemic
       df$count_smooth[j] > df$count_smooth[j+1] & # running mean to ensure that curve is rising (beware of spikes!)
       sum(!is.na(bcp_end[(j+1):(j+15)])) == 0){ # No outbreak flagged during the previous 15 weeks
      bcp_end[j] <- df$date[j]
    } else {
      bcp_end[j] <- NA
    }
  }
  return(bcp_end)
}

a.far <- function(df, truth) {
  state = truth$epidemic
  alarm = df$epidemic
  
  # Summing alarm vector when state is FALSE gives false alerts
  fa = sum(alarm[state==FALSE & !is.na(alarm)])
  # Length of alarm vector when state is FALSE gives number of days with an alarm value
  far = fa / length(alarm[state == FALSE & !is.na(alarm)])
  
  return(far)
  
}

o.detected <- function(df, truth) {
  return(sum(df$epidemic[truth$epidemic == TRUE]) > 0)
}


o.prevented <- function(df, truth) {
  state = truth$epidemic
  alarm = df$epidemic
  length = sum(truth$epidemic)
  prevented = 0
  
  detect = ifelse(sum(alarm[state == TRUE], na.rm = TRUE) > 0, TRUE, FALSE)
  if (detect) {
    first.alarm = min(which(alarm[state==TRUE]==TRUE))
    prevented = (length - first.alarm) / length
  } # if
  
  return(prevented)
}




