lba_pdf = function(t, b, A, v, s){
  #PDF of the LBA model
  
  b_A_tv_ts = (b - A - t*v)/(t*s)
  b_tv_ts = (b - t*v)/(t*s)
  term_1 = v*pnorm(b_A_tv_ts)
  term_2 = s*exp(log(dnorm(b_A_tv_ts,0,1)))
  term_3 = v*pnorm(b_tv_ts)
  term_4 = s*exp(log(dnorm(b_tv_ts,0,1)))
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
  term_3 = ts/A     * exp(log(dnorm(b_A_tv/ts,0,1)))
  term_4 = ts/A     * exp(log(dnorm(b_tv/ts,0,1))) 
  cdf = 1 + term_1 - term_2 + term_3 - term_4
  
  return(cdf)
  
}