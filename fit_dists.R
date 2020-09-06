# This function fits distributions to a *positive-valued* univariate sample.
# It returns a named boolean vector indicating which (if any) of the distributions 
# passed the Cramer-von Mises test at the specified significance level.
# 
fit_dists = function(x, alpha=0.1) {
  require(fitdistrplus)
  require(extraDistr)
  require(triangle)
  
  mymin = min(x);  mymax=max(x)
  if (mymin <= 0) stop('Sorry, this function only accepts positive values.')
  
  fit_ln = fitdist(x, "lnorm", method='mle')
  fit_g  = fitdist(x, "gamma", method='mle')
  fit_w  = fitdist(x, "weibull", method='mle')
  fit_tn = fitdist(x, "tnorm", method='mle', fix.arg=list(a=0, b=Inf),
                       start=list(mean=mean(x), sd=sd(x)), lower=c(0, 0))
  fit_tr = suppressWarnings(fitdist(x, "triangle", method='mle', 
                   start=list(a=0.8*mymin, b=1.1*mymax, c=0.8*mymin + 0.4*(1.1*mymax - 0.8*mymin)), 
                   lower=c(0, 1.01*mymax, mymin), upper=c(0.99*mymin, Inf, mymax)))
  fit_u = suppressWarnings(fitdist(x, "unif", method='mle', 
                                   start=list(min=0.95*mymin, max=1.05*mymax), 
                                   lower=c(0.9*mymin, mymax), upper=c(mymin, 1.1*mymax)))
  
  gofs = gofstat(list(fit_ln, fit_g, fit_w, fit_tn, fit_tr, fit_u), 
                 fitnames=c('Lognormal', 'Gamma', 'Weibull', 'TruncNorm', 
                            'Triangle', 'Uniform'))
  return(gofs$cvm < alpha)
}
