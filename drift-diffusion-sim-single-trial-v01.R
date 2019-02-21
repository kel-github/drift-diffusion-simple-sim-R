# Simulate single trial from single-accumulator drift diffusion model
# written T. Ballard, Feb 2018, free to share and use
# based on R code by Kelly Garner and matlab code by David Sewell 

### clear workspace ###

rm(list = ls())

### load packages ###

library(tidyverse) #for first time installation: "install.packages('tidyverse')"

### set parameter values for simulation ###

drift_rate_mean = .2      #mean change in evidence per second (positive values favour upper bound, negative values favour lower bound)
drift_rate_sd = 1         #sd of change in evidence per second (aka diffusion coefficient)
boundary_separation = 2   #distance between response thresholds, and value of upper threshold (lower threshold has value of 0)
start_point = 1           #level of evidence at start of trial (when half of boundary separation, no bias is present)
ter = 0.2                 #non-decision time (component of response time used for encoding stimulus and generating response)
time_step = 0.001         #length of time between simulation updates (in seconds)

### preallocate objects to store evidence and time information ###

seconds_stored = 100 #length of time (in seconds) for which info is stored. Set this to be much longer than you expect decisions to take

evidence = rep(NA,seconds_stored/time_step) #object used to store level of evidence
time = rep(NA,seconds_stored/time_step)     #object used to store time variable

### set initial levels of variable ###

evidence[1] = start_point   #initial level of evidence
time[1] = ter/2             #time at which decision process starts (assumes encoding process takes 1/2 of non-decision time)

### run simulation ###

ctr=1 #counter used for indexing

#repeat sampling process until one boundary is breached
repeat{
  
  #increment counter by 1
  ctr=ctr+1
  
  #increment evidence by mean rate and random sample of noise. Update is scaled by time step. 
  evidence[ctr] = evidence[ctr-1] + drift_rate_mean*time_step + drift_rate_sd*rnorm(1)*sqrt(time_step)
  
  #update time variable
  time[ctr] = time[ctr-1] + time_step
  
  #check to see if evidence has breached either boundary
  if (evidence[ctr] > boundary_separation | evidence[ctr] < 0){
    break
  }
}

### plot trial ###

tibble(time,evidence) %>%
  filter(!is.na(time)) %>%
  ggplot() +
  geom_line(aes(x=time,y=evidence)) +
  geom_hline(yintercept=0) +
  geom_hline(yintercept=boundary_separation)



