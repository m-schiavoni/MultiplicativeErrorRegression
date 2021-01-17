# function to simulate multiplicative errors with mean=1 and specified CV
library(truncnorm)
library(triangle)
sim_errors = function(err, cv, N, n) {
  
  if (err == 'lognormal') {
    location = log(1 / sqrt(cv^2 + 1))
    shape = sqrt(log(1 + cv^2))
    rand_error = matrix(rlnorm(N*n, meanlog=location, sdlog=shape), nrow=n, ncol=N)
    
  } else if (err == 'gamma') {
    if (cv > 0.8) stop('Max reasonable CV for gamma distribution with mean=1 is 0.8')
    rate = 1/cv^2
    rand_error = matrix(rgamma(N*n, shape=rate, rate=rate), nrow=n, ncol=N)
    
  } else if (err == 'weibull') {
    if (cv > 0.8) stop('Max reasonable CV for weibull distribution with mean=1 is 0.8')
    minfun = function(w) {
      mu_calc = w[2]*gamma(1 + 1/w[1])
      var_calc = w[2]^2 * (gamma(1+2/w[1]) - (gamma(1+1/w[1]))^2)
      return((mu_calc-1)^2 + ((var_calc-cv^2)/cv^2)^2)
    }
    weib = optim(c(2.7, 1.1), minfun)
    rand_error = matrix(rweibull(N*n, shape=weib$par[1], scale=weib$par[2]), nrow=n, ncol=N)
    
  } else if (err == 'truncnorm') {
    if (cv > 0.7) stop('Max reasonable CV for truncated normal distribution with mean=1 is 0.7')
    minfun = function(tn) {
      mu_calc = etruncnorm(a=0, mean=tn[1], sd=tn[2])
      var_calc = vtruncnorm(a=0, mean=tn[1], sd=tn[2])
      return((mu_calc-1)^2 + ((var_calc-cv^2)/cv^2)^2)
    }
    trn = optim(c(1, cv), minfun)
    rand_error = matrix(rtruncnorm(N*n, a=0, mean=trn$par[1], sd=trn$par[2]), nrow=n, ncol=N)
    
  } else if (err == 'triangle') {
    if (cv > 0.7) stop('Max reasonable CV for triangular distribution with mean=1 is 0.7')
    minfun = function(tr) {
      mu_calc = (tr[1]+tr[2]+tr[3])/3
      var_calc = (tr[1]^2 + tr[2]^2 + tr[3]^2 - tr[1]*tr[2] - tr[1]*tr[3] - tr[2]*tr[3]) / 18
      return((mu_calc-1)^2 + ((var_calc-cv^2)/cv^2)^2)
    }
    ineq_con = function(tr) return(c(tr[2]-tr[3], tr[3]-tr[1]))
    triang = suppressMessages(auglag(x0=c(0.2, 2.2, 0.5), fn=minfun, lower=c(0,0,0), hin=ineq_con, localsolver='LBFGS'))
    rand_error = matrix(rtriangle(N*n, a=triang$par[1], b=triang$par[2], c=triang$par[3]), nrow=n, ncol=N)
    
  } else {
    stop('ERROR: invalid error model.')
  }
  
  return(rand_error)
}