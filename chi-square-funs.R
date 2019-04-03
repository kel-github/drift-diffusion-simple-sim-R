data2RTQ = function(data) {
  #Convert empirical data from a single condition to quantiles for diffusion
  #model analysis.
  
  #Assumes data are organized as an N x 2 matrix, where N indexes the number
  #of trials in this condition. Column 1 of the data denotes choice (usually
  #correct or error). Column 2 of the data denotes RT (currently in seconds).
  
  #We can use the "choiceRT" output from simDiff as our data.
  
  Qs = c(.1, .3, .5, .7, .9) #RT quantiles of interest for diffusion modeling
  
  #STEP 1: Setting up the data
  #The first thing we need to do with the data is to separate correct and
  #error responses, as the RT distributions for these will be different.
  
  #Loop that goes through the data and separates correct and error responses.
  dataC = numeric(0)
  dataE = numeric(0)
  
  for (i in 1:dim(data)[1]) {
    #for Trial #1 - #N
    if (data[i, 1] == 1) {
      #Correct response
      dataC = c(dataC,data[i, 2]) #only take the second column entry, since we no longer need accuracy information
    } else {
      #Error response
      dataE = c(dataE,data[i, 2])
    }
  }
  
  #STEP 2: Identifying empirical RT quantiles for each type of response
  #We now need to figure out where the [0.1, 0.3, 0.5, 0.7, and 0.9] RT
  #quantiles are along the time axis for both Correct and Error responses.
  #These quantiles will form the edges of RT "bins" in the data, which are
  #used for parameter estimation using the chi-square method.
  
  #Matlab has an in-built function for extracting quantiles, which I use
  #here, but you can also do this manually by sorting RTs from fastest to
  #slowest (for, say, corrects), and then figuring out where the 10th, 30th,
  #50th, 70th, and 90th percentiles of the data are along the time axis.
  
  binedge_c = quantile(dataC, Qs)
  binedge_e = quantile(dataE, Qs)
  
  #STEP 3: Compute the proportion of observations within each RT bin
  #Estimating parameters using the chi-square method means that we need a way
  #of aggregating the data so that we can compare predictions against
  #observations.
  
  #Due to how quantiles are defined, we know the proportion of observations
  #that are within each bin. For example, the proportion of observations
  #that are between the 0.1 and 0.3 quantiles is equal to 0.2 (i.e., 0.3 -
  #0.1 = 0.2). It follows that the data for any experiment, when binned
  #according to the RT quantiles we are interested in, will be:
  #[0.1, 0.2, 0.2, 0.2, 0.2, 0.1]
  #However, we need to scale these values by the observed correct/error rates
  #so that the proportion of observations sum to 1.
  
  pC = length(dataC) / dim(data)[1] #proportion of corrects
  pE = length(dataE) / dim(data)[1] #proportion of errors
  
  bindata_c = pC * c(.1, .2, .2, .2, .2, .1)
  bindata_e = pE * c(.1, .2, .2, .2, .2, .1)
  
  
  #What we now need to do is compare diffusion model predictions to our
  #empirical data. To do so, we will need to compare the predicted proportion
  #of responses that fall within the empirically defined RT quantiles (i.e.,
  #the bin edges defined by "binedge_c" and "binedge_e"). If the model
  #perfectly matches the data, the expected proportions will equal the
  #observed proportions (i.e., the data vectors defined by "bindata_c" and
  #"bindata_e").
  
  return(list(binedge_c, binedge_e, bindata_c, bindata_e))
}




preds2RTQ = function(preds, binedge_c, binedge_e) {
  #Convert model predictions for a single condition into RT Quantiles
  
  #Assumes preds are organized as an N x 2 matrix, where N indexes the number
  #of simulated trials in this condition. Conventions are otherwise the same
  #as for "data2RTQ"
  
  #We can use the "choiceRT" output from simDiff as our model predictions.
  
  Qc = binedge_c
  #Empirical RT quantiles for correct responses.
  Qe = binedge_e
  #Empirical RT quantiles for error responses.
  
  #We now need to convert our model predictions (here, based on simulation),
  #to the same RT-binned format as the data. The procedure is similar to how
  #we dealt with data, but with a few changes.
  
  #Loop that goes through the predictions and separates correct and error
  #trials.
  
  predC = numeric(0)
  predE = numeric(0)
  
  for (i in 1:dim(preds)[1]) {
    #for Trial #1 - #N
    if (preds[i, 1] == 1) {
      #Correct response
      predC = c(predC,preds[i, 2]) #only take the second column entry, since we no longer need accuracy information
    } else {
      #Error response
      predE = c(predE,preds[i, 2])
    }
  }
  
  #Convert frequencies to proportions (for later).
  predPC = length(predC) / dim(preds)[1] #Predicted P(correct)
  predPE = length(predE) / dim(preds)[1] #Predicted P(error)
  
  #We now identify the proportion of *predicted* trials that fall within the
  #*empirically-defined* bin edges for both correct and error trials.
  
  #First, we need to sort the predicted RTs to make life easier.
  predC = sort(predC)
  predE = sort(predE)
  
  #Now we see how many simulated trials fell between each of the bin edges.
  #This is a matter of comparing the numbers of trials that fall between
  #successive bin edges.
  
  bindata_pc = c(
    sum(predC <= Qc[1]),
    sum(predC <= Qc[2]) - sum(predC <= Qc[1]),
    sum(predC <= Qc[3]) - sum(predC <= Qc[2]),
    sum(predC <= Qc[4]) - sum(predC <= Qc[3]),
    sum(predC <= Qc[5]) - sum(predC <= Qc[4]),
    sum(predC > Qc[5]) / dim(preds)[1]
  )
  bindata_pe = c(
    sum(predE <= Qe[1]),
    sum(predE <= Qe[2]) - sum(predE <= Qe[1]),
    sum(predE <= Qe[3]) - sum(predE <= Qe[2]),
    sum(predE <= Qe[4]) - sum(predE <= Qe[3]),
    sum(predE <= Qe[5]) - sum(predE <= Qe[4]),
    sum(predE > Qe[5]) / dim(preds)[1]
  )
  
  return(list(bindata_pc, bindata_pe))
  
}



chisqFit = function(data, preds, N) {
  #Computes a chi-square statistic comparing diffusion model predictions
  #against data. N denotes the number of trials per condition.
  
  #The data are assumed to be an N x 12 matrix, where N denotes the number of
  #experimental conditions (we're just doing 1 condition for our demo), and
  #the RTQ-binned proportions of observations for correct responses followed
  #by error responses.
  
  #Model predictions ("preds") are assumed to be formatted in the same way.
  #The simulated output from the Diffusion Model loop can be converted to
  #quantiles using the "data2RTQ" function (just as we did for the data).
  
  #I have included a double-summation to allow for the possibility of
  #multiple conditions.
  
  #Set minimum value of predicted counts to a small constant
  #to avoid dividing by 0
  preds = pmax(preds,1e-10)
  
  #calculate chi-square
  x2 = N * sum( sum( ( (data - preds ) ^ 2) / preds ) )
  return(x2)
}
