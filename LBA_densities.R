
#PDF of the LBA model
lba_pdf = function(t, b, A, v, s){
  
  b_A_tv_ts = (b - A - t*v)/(t*s)
  b_tv_ts = (b - t*v)/(t*s)
  term_1 = v*pnorm(b_A_tv_ts)
  term_2 = s*dnorm(b_A_tv_ts,0,1)
  term_3 = v*pnorm(b_tv_ts)
  term_4 = s*dnorm(b_tv_ts,0,1)
  pdf = (1/A)*(-term_1 + term_2 + term_3 - term_4)
  
  return(pdf)
}

#CDF of the LBA model
lba_cdf = function(t, b, A, v, s){
  
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

#Function to calculate the negative summed log-likelihood for LBA
lba_log_likelihood = function(parms,data){
  
  #print parameters being considered
  print(parms)
  
  #unpack parameters
  mean_drift_rate = parms[1:2]
  sd_drift_rate = 1  #fixed to 1 for scaling purposes
  threshold = parms[3]
  max_start_point = parms[4]
  non_decision_time = parms[5]
  
  #calculate raw threshold
  raw_threshold = max_start_point + threshold
  
  #initialise variable to store likelihoods
  likelihood = rep(NA,nrow(data))
  
  #loop through data, calculating likelihood for each observation at a time
  for (i in 1:nrow(data)){
    
    #calculate decision time by subtracting non_decision_time parameter from observed RT
    decision_time = data[i,2] - non_decision_time
    
    #only consider responses where predicted decision time is greater than 0. If less than 0,
    #we assign a very low likelihood (see below)
    if(decision_time > 0){
      
      #evaluate the pdf of the observed response, and the product of 1-cdf for all other responses. 
      cdf = 1;
      for(j in 1:length(mean_drift_rate)){
        if(data[i,1] == j){
          pdf = lba_pdf(decision_time, raw_threshold , max_start_point, mean_drift_rate[j],  sd_drift_rate)
        }else{
          cdf = (1-lba_cdf(decision_time, raw_threshold , max_start_point, mean_drift_rate[j],  sd_drift_rate)) * cdf;
        }
      }
      #calculate the probability of all responses having a negative drift rate (in this case the RT is undefined)
      prob_neg = 1;
      for(j in 1:length(mean_drift_rate)){
        prob_neg = pnorm(-mean_drift_rate[j]/ sd_drift_rate ) * prob_neg;
      }
      
      #calculate the likelihood of the observed response and response time by multiplying the pdf of the observed
      #response with the product of 1-cdf of all the other responses.
      likelihood[i] = pdf*cdf;
      
      #normalise the likelihood by the probability that at least one accumulator has a positive drift rate
      likelihood[i] = likelihood[i]/(1-prob_neg);
      
      #set lower bound on the likelihood to avoid underflow issues
      likelihood[i] = max(likelihood[i],1e-10)
      
    }else{
      #if decision time is less than 0, likelihood is set to very low value.
      likelihood[i] = 1e-10;
    }
  }
  
  #calculate the negative summed log likelihood by summing the log of each likelihood and multiplying the result by -1
  out = -sum(log(likelihood));
  
  return(out);
}

#The function below is identical to the one above but contains code to plot the predictions of the model under
#the parameters being considered

lba_log_likelihood_w_plot = function(parms,data){
  print(parms)
  #unpack parameters
  mean_drift_rate = parms[1:2]
  sd_drift_rate = 1  #fixed to 1 for scaling purposes
  threshold = parms[3]
  max_start_point = parms[4]
  non_decision_time = parms[5]
  
  raw_threshold = max_start_point + threshold
  
  likelihood = rep(NA,nrow(data))
  
  for (i in 1:nrow(data)){
    decision_time = data[i,2] - non_decision_time
    #only consider responses where predicted decision time is greater than 0
    if(decision_time > 0){
      cdf = 1;
      for(j in 1:length(mean_drift_rate)){
        if(data[i,1] == j){
          pdf = lba_pdf(decision_time, raw_threshold , max_start_point, mean_drift_rate[j],  sd_drift_rate)
        }else{
          cdf = (1-lba_cdf(decision_time, raw_threshold , max_start_point, mean_drift_rate[j],  sd_drift_rate)) * cdf;
        }
      }
      prob_neg = 1;
      for(j in 1:length(mean_drift_rate)){
        prob_neg = pnorm(-mean_drift_rate[j]/ sd_drift_rate ) * prob_neg;
      }
      likelihood[i] = pdf*cdf;
      likelihood[i] = likelihood[i]/(1-prob_neg);
      likelihood[i] = max(likelihood[i],1e-10)
      
    }else{
      likelihood[i] = 1e-10;
    }
  }
  out = -sum(log(likelihood));
  
  #quick way to randomly generate observations
  vs = mean_drift_rate
  s = sd_drift_rate
  A = max_start_point
  n = nrow(data)
  b = threshold + max_start_point
  t0 = non_decision_time
  st0 = 0
  
  n.with.extras=ceiling(n*(1+3*prod(pnorm(-vs))))
  drifts=matrix(rnorm(mean=vs,sd=s,n=n.with.extras*length(vs)),ncol=length(vs),byrow=TRUE)
  
  repeat {
    drifts=rbind(drifts,matrix(rnorm(mean=vs,sd=s,n=n.with.extras*length(vs)),ncol=length(vs),byrow=TRUE))
    tmp=apply(drifts,1,function(x) any(x>0))
    drifts=drifts[tmp,]
    if (nrow(drifts)>=n) break
  }
  
  drifts=drifts[1:n,]
  drifts[drifts<0]=0
  starts=matrix(runif(min=0,max=A,n=n*length(vs)),ncol=length(vs),byrow=TRUE)
  ttf=t((b-t(starts)))/drifts
  RT=apply(ttf,1,min)+t0+runif(min=-st0/2,max=+st0/2,n=n)
  resp=apply(ttf,1,which.min)
  predicted = data.frame(resp=resp,RT=RT)
  
  data$source = "observed"
  predicted$source = "predicted"
  
  plot = bind_rows(data,predicted) %>%
    filter(RT<5) %>%
    ggplot() +
    geom_density(aes(x=RT,group=source,colour=source),alpha=0.1) +
    facet_grid(.~resp)
  
  print(plot)
  return(out);
}
