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
generate_rts_w_single_drift <- function(v, a, Ter, z, s, dt, nTrials, t0){
  # use Euler's method to generate sample paths of the diffusion process
  # each trial of the simulation is constructed as a loop that will generate simulated 
  # sample paths. This loop is then applied nTrials times


  get.trials <- function(z, t0, a, Ter, v, dt){
    # the following loop controls the evidence accumulation process. At each time step,
    # the current evidence total (x) is compared againt the value of the two absorbing 
    # boundaries (a and 0). Sample steps occur until the evidence total crosses
    # one of the boundaries. Each sample of info changes the evidence total by
    # an amount that is fixed across time steps (v * dt, mu drift rate * time step),
    # plus a scaled random draw from a normal distribution with mean 0 and sd s. Once 
    # evidence total crosses a boundary, the time steps and the boundary crossed are 
    # recorded
    # initialise evidence
    x = z
    t = t0
    repeat{
      x = x + v*dt + s*rnorm(1)*sqrt(dt)
      t = t + dt
      if (x > a | x < 0){
        break
      }
    }
    # a threshold has been reached (hence termination of while loop, so assign response and get RT)
    if ( x > a )
      resp = 1
    else if ( x < 0 ){
      resp = 0
    }
    RT = t + Ter
    
    out = data.frame(resp = resp, RT = RT) # assign the variables to a list
    out # output the list
  }
  
  # initialise a vector to collect responses
  z = rep(z, times = nTrials)
  tmp = lapply(z, get.trials, t0 = t0, a = a, Ter = Ter) # apply get.trials function nTrials times
  data = do.call(rbind, tmp) # make into a nice data frame
  data
}
  

generate_rts_w_independent_drifts <- function(v_i, v_j, a, Ter, z, s, dt, nTrials, t0){
  # use Euler's method to generate sample paths of the diffusion process
  # each trial of the simulation is constructed as a loop that will generate simulated 
  # sample paths, one for each independent accumulator (i & j). 
  # This loop is then applied nTrials times
  
  
  get.trials <- function(z, t0, a, Ter, v_i, v_j, dt){
    # the following loop controls the evidence accumulation process. At each time step,
    # the current evidence total (x) is compared againt the value of the absorbing 
    # boundaries for each accumulator (i: a_i, j: a_j). Sample steps occur until the evidence total crosses
    # one of the boundaries. Each sample of info changes the evidence total by
    # an amount that is fixed across time steps (v_y * dt, mu drift rate * time step),
    # plus a scaled random draw from a normal distribution with mean 0 and sd s. Once 
    # evidence total crosses a boundary, the time steps and the boundary crossed are 
    # recorded
    # initialise evidence
    x_i = z
    x_j = z
    t = t0
    repeat{
      x_i = x_i + v_i*dt + s*rnorm(1)*sqrt(dt)
      x_j = x_j + v_j*dt + s*rnorm(1)*sqrt(dt)
      t = t + dt
      if (x_i > a | x_j > a){
        break
      }
    }
    # a threshold has been reached (hence termination of while loop, so assign response and get RT)
    if ( x_i > a )
      resp = 1
    else {
      resp = 0
    }
    RT = t + Ter
    
    out = data.frame(resp = resp, RT = RT) # assign the variables to a list
    out # output the list
  }
  
  # initialise a vector to collect responses
  z = rep(z, times = nTrials)
  tmp = lapply(z, get.trials, t0 = t0, a = a, Ter = Ter, v_i = v_i, v_j = v_j, dt = dt) # apply get.trials function nTrials times
  data = do.call(rbind, tmp) # make into a nice data frame
  data
}


# GENERATE DATA FROM SINGLE DRIFT MODEL
#--------------------------------------------------------------------------
# first, set the values for the core parameters in the diffusion model (with single drift implementation):
v = 0.2 # drift rate
a = .04 # boundary separation
Ter = 0.2 # non-decision time
z = a/2 # for an unbiased start point

# next, set the standard deviation of within trial noise, traditionally set to 0.1 
# (Ratcliff, 1978), but others set s = 1. Either is fine, as this parameter “scales” the 
# other parameters in the model (i.e., diffusion model parameters are only defined on a 
# ratio scale, and so you can get the same predictions from different combinations of 
# parameter values---if one parameter were to be mulitplied by a constant, the same 
# predictions would obtain if all other parameters were multiplied (or divided) by the same value).
s = 0.1 # diffusion coefficient
dt = .001 # Time steps in seconds

# set the number of trials
nTrials = 10000
t0 = 0 # time zero
#
# now run the functions to get simulated data
sing.accum.dat = generate_rts_w_single_drift(v = v, a = a, Ter = Ter, z = z, s = s, dt = dt, nTrials = nTrials, t0=t0)

# _____________________________________________________________________________________________________

# GENERATE DATA FROM INDEPENDENT DRIFTS MODEL
#--------------------------------------------------------------------------
# first, set the values for the core parameters in the diffusion model (with independent drift implementation):
v_i = 0.19 # first of two independent drift rates (drift i)
v_j = 0.15 # same for drift j
a = .04 # boundary separation 
Ter = 0.2 # non-decision time
z = a/2 # for an unbiased start point

ind.accum.dat = generate_rts_w_independent_drifts(v_i = v_i, v_j = v_j, a = a, Ter = Ter, z = z, s = s, dt = dt, nTrials = nTrials, t0 = t0)
# (v_i, v_j, a, Ter, z, s, dt, nTrials, t0)
#------------------------------------------------------------------------------

# PLOT OUTPUT DATAs AS HISTOGRAMS
plot.RTs <- function(data){
  
  plot <- ggplot(data, aes(x=RT, color=as.factor(resp))) +
    geom_histogram(binwidth=.02) +
    facet_wrap(~as.factor(resp))
  plot
}

plot.RTs(sing.accum.dat)
plot.RTs(ind.accum.dat)
