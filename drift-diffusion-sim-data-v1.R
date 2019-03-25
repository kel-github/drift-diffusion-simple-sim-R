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
# install.packages("ggplot2")
library(ggplot2)
#--------------------------------------------------------------------------

# DEFINE FUNCTIONS
#--------------------------------------------------------------------------
generate_rts_w_single_drift <- function(drift_rate, boundary, non_decision_time, start_point, noise_sdev, time_step, nTrials, t0){
  # use Euler's method to generate sample paths of the diffusion process
  # each trial of the simulation is constructed as a loop that will generate simulated 
  # sample paths (i.e. accumulations of evidence). This loop is then applied nTrials times


  get.trials <- function(start_point, t0, boundary, non_decision_time, drift_rate, time_step, noise_sdev){
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
    T_evidence = start_point
    t = t0
    repeat{
      T_evidence = T_evidence + drift_rate*time_step + rnorm(1, mean=0, sd = noise_sdev)*sqrt(time_step)
      t = t + time_step
      if (T_evidence > boundary | T_evidence < 0){
        break
      }
    }
    # a threshold has been reached (hence termination of while loop, so assign response and get RT)
    if ( T_evidence > boundary )
      resp = 1
    else if ( T_evidence < 0 ){
      resp = 0
    }
    RT = t + non_decision_time
    
    out = data.frame(resp = resp, RT = RT) # assign the variables to a list
    out # output the list
  }
  
  # initialise a vector to collect responses
  start_point = rep(start_point, times = nTrials)
  tmp = lapply(start_point, get.trials, t0 = t0, boundary = boundary, 
               non_decision_time = non_decision_time, drift_rate = drift_rate, 
               time_step = time_step, noise_sdev = noise_sdev) # apply get.trials function nTrials times
  data = do.call(rbind, tmp) # make into a nice data frame
  data
}
  
# GENERATE DATA FROM SINGLE DRIFT MODEL
#--------------------------------------------------------------------------
# first, set the values for the core parameters in the diffusion model (with single drift implementation):
drift_rate        = 0.4 # drift rate
boundary          = .04 # boundary separation
non_decision_time = 0.3 # non-decision time
start_point       = .02 #

# next, set the standard deviation of within trial noise, traditionally set to 0.1 
# (Ratcliff, 1978), but others set noise_sdev = 1. Either is fine, as this parameter “scales” the 
# other parameters in the model (i.e., diffusion model parameters are only defined on a 
# ratio scale, and so you can get the same predictions from different combinations of 
# parameter values---if one parameter were to be mulitplied by a constant, the same 
# predictions would obtain if all other parameters were multiplied (or divided) by the same value).
noise_sdev = 0.1 # diffusion coefficient
time_step  = .001 # Time steps in seconds

# set the number of trials
nTrials = 10000
t0 = 0 # time zero
#
# Start the clock! (if you want to time, uncomment line below and beneath 'Stop the clock')
# ptm <- proc.time()
predicted.data = generate_rts_w_single_drift(drift_rate = drift_rate, boundary = boundary,
                                             non_decision_time = non_decision_time, 
                                             start_point = start_point, noise_sdev = noise_sdev, 
                                             time_step = time_step, nTrials = nTrials, t0=t0)
# Stop the clock
# proc.time() - ptm
# _____________________________________________________________________________________________________
# load the previously made distribution
load("observed_data.Rda")
# PLOT OUTPUT DATAs AS HISTOGRAMS
plot.observed_vs_predicted_RTs <- function(observed, predicted){
  
par(mfrow=c(2,1), mar=c(3, 3, 1, 1))
 
  get.plot <- function(observed, predicted, response){
    tmp = with(predicted, density(RT[resp==response]))
    with(observed, hist(RT[resp==response], prob=TRUE, main = paste("resp = ", response, sep=""),
                        xlab="RT", xlim=c(0,1), ylim=c(0,max(tmp$y)), col=scales::alpha('skyblue',.5)))
    with(predicted, lines(density(RT[resp==response]), col="#F24D29", lwd=2))
  }
  
  get.plot(observed = observed, predicted = predicted, response = 0)
  get.plot(observed = observed, predicted = predicted, response = 1)
}

plot.observed_vs_predicted_RTs(observed.data, predicted.data)
