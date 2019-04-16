# Simulate choice RT data using drift diffusion model (quick and dirty method)
# written K. Garner, Feb 2019, free to share and use
# based on matlab code by David Sewell 
rm(list = ls())
# code structure:
# 1. install packages (if necessary), load packages
# 2. define functions that simulate the model
# 3. set intital parameters and run
# 4. plot output RTs as histograms
#--------------------------------------------------------------------------

# install packages (if necessary), load packages
# install required packages
# install.packages("tidyverse")
library(tidyverse)
#--------------------------------------------------------------------------

# DEFINE FUNCTIONS
#--------------------------------------------------------------------------
generate_rts <- function(max_start_point, mean_drift_rates, sd_drift_rates, threshold, non_decision_time, nTrials){
  # generates choices and response times for N trials
  
  
  get.trials <- function(max_start_point, mean_drift_rates, sd_drift_rates, threshold, non_decision_time){
    #this function simulates the evidence accumulation process according to the LBA.
    
    #first, get number of response alternatives, which can two or more.
    n_responses = length(mean_drift_rates)
    
    #sample starting point from uniform distribution bounded between 0 and the max starting point
    start_point = runif(n=1,min=0,max=max_start_point)
    
    #sample drift rates for each accumulator from normal distribution.
    #need to ensure at least one drift rate is positive, because otherwise the response time is infinite
    #so we repeat the sampling process until at least one drift rate is positive.
    repeat{
      #sample drift rates
      drift_rates = rnorm(n=n_responses,mean=mean_drift_rates,sd=sd_drift_rates)
      #check to see if at least one is positive
      if(any(drift_rates > 0) ){
        break
      }
    }
    
    #calculate raw threshold. Threshold is typically expressed as the raw threshold minus
    #the maximum starting point. So we need to calculate the actual threshold (i.e., the distance
    #between 0 evidence and the threshold) by summing the threshold and max starting point. 
    raw_threshold = threshold + max_start_point 
    
    #calculate evidence required to reach threshold. This is just the difference between the raw_threshold
    #and the randomly sampled starting point
    evidence_required = raw_threshold - start_point
    
    #calculate time required for each acculator to reach threshold by dividing the evidence required 
    #(i.e. the distance) by the rate of evidence accumulation
    time_required = evidence_required / drift_rates
    
    #identify the minimum time time required. This value becomes the decision time, which represents the
    #time taken for the first accumulator to hit the threshold. Here we disregard values that are negative,
    #because these are only generated when the drift rate is negative. In this situation, the evidence 
    #never breaches the threshold
    decision_time = min(time_required[time_required > 0])
    
    #compute the response time by summing the decision time and the non decision time.
    RT = decision_time + non_decision_time
    
    #identify the response that breaches threshold
    resp = which(time_required == decision_time)
    
    #assign the response and RT variables to a data frame
    out = data.frame(resp = resp, RT = RT) 
      
    return(out)
  }
  
  # initialise a vector to collect responses
  max_start_point = rep(max_start_point, times = nTrials)
  tmp = lapply(max_start_point, get.trials, mean_drift_rates = mean_drift_rates, 
               sd_drift_rates = sd_drift_rates, threshold = threshold, non_decision_time = non_decision_time) # apply get.trials function nTrials times
  data = do.call(rbind, tmp) # make into a nice data frame
  return(data)
}


# GENERATE DATA FROM THE LBA
#--------------------------------------------------------------------------
# first, set the values for the core parameters for the LBA
mean_drift_rates        = c(0.7,1.3)  # drift rate for each accumlator (here we assume there are two responses)
sd_drift_rates          = 1           # standard deviation of drift rates
threshold               = 1.2         # response threshold
max_start_point         = 0.6         # maximum starting point
non_decision_time       = 0.2         # non-decision time

# set the number of trials
nTrials = 10000
#
# Start the clock! (if you want to time, uncomment line below and beneath 'Stop the clock')
# ptm <- proc.time()
predicted.data = generate_rts(mean_drift_rates = mean_drift_rates, 
                              sd_drift_rates = sd_drift_rates, 
                              threshold = threshold,
                              max_start_point = max_start_point,
                              non_decision_time = non_decision_time, 
                              nTrials = nTrials)

# Stop the clock
# proc.time() - ptm
# _____________________________________________________________________________________________________
# load the previously generated data
load("observed_data_lba.RData")

# PLOT OUTPUT DATAs AS DENSITIES
plot.observed_vs_predicted_RTs <- function(observed, predicted){
  
  observed$source = "observed"
  predicted$source = "predicted"
  
  bind_rows(observed,predicted) %>%
    filter(RT<5) %>%
    ggplot() +
    geom_density(aes(x=RT,group=source,colour=source),alpha=0.1) +
    facet_grid(.~resp)
  
}

plot.observed_vs_predicted_RTs(observed.data, predicted.data)

# IMPLEMENT MAXIMUM LIKELIHOOD PARAMETER ESTIMATION
#-------------------------------------------------------------------
#source functions needed to calculate the likelihood
source("LBA_densities.R")

#set initial parameter values
initial_parms = c(mean_drift_rate = c(1,1),
          threshold = 1,
          max_start_point = 0.1,
          non_decision_time = 0.3)

#optimise to find parameter values that are most likely given that data
optim(par = initial_parms,
      fn = lba_log_likelihood,
      data = observed.data,
      lower = c(0,0, 0, 0.001, 0.05),
      upper = c(5,5, 5, 5, 1),
      control = list(trace=6)) 

