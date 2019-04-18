# Simulate choice RT data using linear ballistic accumulator (LBA) model
#
# Written by T. Ballard, Apr 2019, free to share and use.
#------------------------------------------------------------------------------

rm(list = ls())

library(tidyverse)


# DEFINE FUNCTIONS
#------------------------------------------------------------------------------
do.trial <- function(max_start_point, mean_drift_rates, sd_drift_rates,
                     threshold, non_decision_time) {
  #' Simulates the evidence accumulation process according to the LBA model
  #'
  #' @param max_start_point Each accumulator starts at a uniformly random point
  #'   between zero and this value
  #' @param mean_drift_rates A vector of drift rates, one for each accumulator.
  #'   The length of this parameter determines the number of accumulators.
  #' @param sd_drift_rates The standard deviation of the drift rates.
  #' @param threshold The decision threshold for each accumulator.
  #' @param non_decision_time Constant time added for non-decision processes.

  # First, get number of response alternatives, which can be two or more.
  n_responses <- length(mean_drift_rates)

  # Sample starting point from uniform distribution bounded between 0 and the
  # max starting point.
  start_point <- runif(n = 1, min = 0, max = max_start_point)

  # Sample drift rates for each accumulator from normal distribution. Need to
  # ensure at least one drift rate is positive (otherwise the response time is
  # infinite) so repeat the sampling process until at least one drift rate is
  # positive.
  repeat {
    drift_rates <- rnorm(n = n_responses, mean = mean_drift_rates,
                         sd = sd_drift_rates)
    if (any(drift_rates > 0)) {
      break
    }
  }

  # Calculate raw threshold. Threshold is typically expressed as the raw
  # threshold minus the maximum starting point. So we need to calculate the
  # actual threshold (i.e., the distance between 0 evidence and the threshold)
  # by summing the threshold and max starting point.
  raw_threshold <- threshold + max_start_point

  # Calculate evidence required to reach threshold. This is just the difference
  # between the raw_threshold and the randomly sampled starting point.
  evidence_required <- raw_threshold - start_point

  # Calculate time required for each accumulator to reach threshold by dividing
  # the evidence required (i.e. the distance) by the rate of evidence
  # accumulation.
  time_required <- evidence_required / drift_rates

  # Identify the minimum time required. This value becomes the decision time,
  # which represents the time taken for the first accumulator to hit the
  # threshold. Here we disregard values that are negative, because these are
  # only generated when the drift rate is negative. In this situation, the
  # evidence never breaches the threshold.
  decision_time <- min(time_required[time_required > 0])

  # Compute the response time by summing the decision time and the non decision
  # time.
  RT <- decision_time + non_decision_time

  # Identify the response that breaches threshold.
  resp <- which(time_required == decision_time)

  # Assign the response and RT variables to a data frame.
  out <- data.frame(resp = resp, RT = RT)

  return(out)
}


collect.trials <- function(n_trials, max_start_point, mean_drift_rates,
                           sd_drift_rates, threshold, non_decision_time) {
  #' Generates choices and response times for N trials
  #'
  #' @param n_trials The number of trials to run.

  data <- replicate(n_trials,
    do.trial(max_start_point = max_start_point,
             mean_drift_rates = mean_drift_rates,
             sd_drift_rates = sd_drift_rates,
             threshold = threshold,
             non_decision_time = non_decision_time),
    simplify=TRUE)
  return(as.data.frame(t(data)))
}


plot.observed_vs_predicted_RTs <- function(observed, predicted) {
  #' Plot output data as densities

  observed$source = "observed"
  predicted$source = "predicted"

  bind_rows(observed, predicted) %>%
    filter(RT < 5) %>%
    ggplot() +
    geom_density(aes(x = RT, group = source, colour = source), alpha=0.1) +
    facet_grid(.~resp)

}


# GENERATE DATA FROM THE LBA
#------------------------------------------------------------------------------
mean_drift_rates        = c(0.7, 1.3)  # drift rate for two accumulators
sd_drift_rates          = 1            # standard deviation of drift rates
threshold               = 1.2          # response threshold
max_start_point         = 0.6
non_decision_time       = 0.2

n_trials = 10000
time.trials = FALSE

# Start the clock!
if (time.trials) {
  ptm <- proc.time()
}
predicted.data <- collect.trials(
  n_trials,
  max_start_point = max_start_point,
  mean_drift_rates = mean_drift_rates,
  sd_drift_rates = sd_drift_rates,
  threshold = threshold,
  non_decision_time = non_decision_time)

# Stop the clock
if (time.trials) {
  proc.time() - ptm
}


# LOAD AND PLOT THE PREVIOUSLY GENERATED DATA
#------------------------------------------------------------------------------
load("observed_data_lba.RData")

plot.observed_vs_predicted_RTs(observed.data, predicted.data)


# IMPLEMENT MAXIMUM LIKELIHOOD PARAMETER ESTIMATION
#------------------------------------------------------------------------------
# Source functions needed to calculate the likelihood.
source("LBA_densities.R")

# Set initial parameter values.
initial_parms = c(mean_drift_rate = c(1, 1),
                  threshold = 1,
                  max_start_point = 0.1,
                  non_decision_time = 0.3)

# Optimise to find parameter values that are most likely given that data.
optim(par = initial_parms,
      fn = lba_log_likelihood,
      data = observed.data,
      lower = c(0,0, 0, 0.001, 0.05),
      upper = c(5,5, 5, 5, 1),
      control = list(trace=6))

