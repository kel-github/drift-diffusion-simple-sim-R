# Simulate choice RT data using drift diffusion model (quick and dirty method)
# written K. Garner, Feb 2019, free to share and use
# based on matlab code by David Sewell 
rm(list = ls())
# code structure:
# 1. install packages (if necessary), load packages
# 2. define diffusion model functions - one for single accumulator model, one for two independent accumulators
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
generate_rts <- function(max_start_point, mean_drift_rates, sd_drift_rates, thresholds, non_decision_time, nTrials){
  # use Euler's method to generate sample paths of the diffusion process
  # each trial of the simulation is constructed as a loop that will generate simulated 
  # sample paths (i.e. accumulations of evidence). This loop is then applied nTrials times
  
  
  get.trials <- function(max_start_point, mean_drift_rates, sd_drift_rates, thresholds, non_decision_time){
    # the following loop controls the evidence accumulation process. At each time_step,
    # the current evidence total (T_evidene) is compared againt the value of the two absorbing 
    # boundaries (boundary and 0). Sample steps occur until the evidence total crosses
    # one of the boundaries. Each sample of info changes the evidence total by
    # an amount that is fixed across time steps (drift_rate * time_steps),
    # plus a scaled random draw from a normal distribution with mean 0 and sd noise_sdev 
    # (a value of 0.1 or 1 are often used). Once 
    # evidence total crosses a boundary, the time steps and the boundary crossed are 
    # recorded, the non-decision time is then added to the time steps to make the predicted RT
    # initialise evidence
    
    #get number of responses
    n_responses = length(mean_drift_rates)
    
    #sample starting point
    start_point = runif(n=1,min=0,max=max_start_point)
    
    #sample drift rates
    #Need to ensure at least one drift rate is positive, because otherwise the response time is infinite
    repeat{
      drift_rates = rnorm(n=n_responses,mean=mean_drift_rates,sd=sd_drift_rates)
      
      if(any(drift_rates > 0) ){
        break
      }
      
    }
    
    #calculate evidence required to reach threshold
    evidence_required = thresholds - start_point
    
    #calculate time required for each acculator to reach threshold
    time_required = evidence_required / drift_rates
    
    decision_time = min(time_required[time_required > 0])
    
    RT = decision_time + non_decision_time
    
    resp = which(time_required == decision_time)
    
    out = data.frame(resp = resp, RT = RT) # assign the variables to a list
      
    return(out)
  }
  
  # initialise a vector to collect responses
  max_start_point = rep(max_start_point, times = nTrials)
  tmp = lapply(max_start_point, get.trials, mean_drift_rates = mean_drift_rates, 
               sd_drift_rates = sd_drift_rates, thresholds = threshold, non_decision_time = non_decision_time) # apply get.trials function nTrials times
  data = do.call(rbind, tmp) # make into a nice data frame
  return(data)
}

# GENERATE DATA FROM SINGLE DRIFT MODEL
#--------------------------------------------------------------------------
# first, set the values for the core parameters in the diffusion model (with single drift implementation):
mean_drift_rates        = c(0.7,1.3)  # drift rate
sd_drift_rates          = c(1,1)  # standard deviation of drift rate
threshold               = c(1.3,2)
max_start_point         = 0.5
non_decision_time       = 0.2 # non-decision time

# set the number of trials
nTrials = 10000
#
# Start the clock! (if you want to time, uncomment line below and beneath 'Stop the clock')
# ptm <- proc.time()
predicted.data = generate_rts(mean_drift_rates = mean_drift_rates, 
                              sd_drift_rates = sd_drift_rates, 
                              thresholds = thresholds,
                              max_start_point = max_start_point,
                              non_decision_time = non_decision_time, 
                              nTrials = nTrials)
# Stop the clock
# proc.time() - ptm
# _____________________________________________________________________________________________________
# load the previously made distribution
load("observed_data_lba.RData")
# PLOT OUTPUT DATAs AS HISTOGRAMS
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
#source functions needed to calculate chi-squared
source("LBA_densities.R")

initial_parms = c(mean_drift_rate = c(1,1),
          threshold = c(1,1),
          max_start_point = 0.5,
          non_decision_time = 0.3)

lba_log_likelihood(initial_parms,predicted.data)

optim(par = initial_parms,
      fn = lba_log_likelihood,
      data = observed.data,
      lower = c(-10,-10, 0, 0, 0.001, 0.05),
      upper = c(10,10, 10, 10, 10, 1),
      control = list(trace=6)) 
      

parms=c(10.00,            -10.00   ,           0.00  ,           10.00      ,        0.00  ,            0.05)
