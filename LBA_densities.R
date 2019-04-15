lba_pdf = function(t, b, A, v, s){
  #PDF of the LBA model
  
  b_A_tv_ts = (b - A - t*v)/(t*s)
  b_tv_ts = (b - t*v)/(t*s)
  term_1 = v*pnorm(b_A_tv_ts)
  term_2 = s*dnorm(b_A_tv_ts,0,1)
  term_3 = v*pnorm(b_tv_ts)
  term_4 = s*dnorm(b_tv_ts,0,1)
  pdf = (1/A)*(-term_1 + term_2 + term_3 - term_4)
  
  return(pdf)
}

lba_cdf = function(t, b, A, v, s){
  #CDF of the LBA model
  
  b_A_tv = b - A - t*v;
  b_tv = b - t*v;
  ts = t*s;
  term_1 = b_A_tv/A * pnorm(b_A_tv/ts)	
  term_2 = b_tv/A   * pnorm(b_tv/ts)
  term_3 = ts/A     * dnorm(b_A_tv/ts,0,1)
  term_4 = ts/A     * dnorm(b_tv/ts,0,1)
  cdf = 1 + term_1 - term_2 + term_3 - term_4
  
  return(cdf)
  
}

lba_log_likelihood = function(parms,data){
  print(parms)
  #unpack parameters
  mean_drift_rate = parms[1:2]
  sd_drift_rate = c(1,1) #fixed to 1 for scaling purposes
  threshold = parms[3:4]
  max_start_point = parms[5]
  non_decision_time = parms[6]
  
  raw_threshold = max_start_point + threshold
  
  likelihood = rep(NA,nrow(data))
  
  for (i in 1:nrow(data)){
    decision_time = data[i,2] - non_decision_time
    #only consider responses where predicted decision time is greater than 0
    if(decision_time > 0){
        cdf = 1;
        for(j in 1:length(mean_drift_rate)){
        if(data[i,1] == j){
          pdf = lba_pdf(decision_time, raw_threshold[j] , max_start_point, mean_drift_rate[j],  sd_drift_rate[j])
        }else{
          cdf = (1-lba_cdf(decision_time, raw_threshold[j], max_start_point, mean_drift_rate[j],  sd_drift_rate[j])) * cdf;
        }
      }
      prob_neg = 1;
      for(j in 1:length(mean_drift_rate)){
        prob_neg = pnorm(-mean_drift_rate[j]/ sd_drift_rate[j] ) * prob_neg;
      }
      likelihood[i] = pdf*cdf;
      likelihood[i] = likelihood[i]/(1-prob_neg);
      if(likelihood[i] < 1e-10){
        likelihood[i] = 1e-10;
      }
      
    }else{
      likelihood[i] = 1e-10;
    }
  }
  out = -sum(log(likelihood));
  return(out);
}
