# This function plots density curves for multiplicative error distributions with
# mean=1 and specified CV for demonstration/visualization purposes.
#
plot_dists = function(cv=0.6) {
  require(truncnorm)
  require(triangle)
  require(nloptr)
  
  xvec = seq(0.001, 3.501, length.out=350)
  
  # Lognormal
  location = log(1 / sqrt(cv^2 + 1))
  shape = sqrt(log(1 + cv^2))
  ylnorm = dlnorm(xvec, location, shape)
  
  # Gamma
  rate = 1/cv^2
  ygamma = dgamma(xvec, shape=rate, rate=rate)
  
  # Weibull
  minfun = function(par) {
    mu_calc = par[2]*gamma(1 + 1/par[1])
    cv_calc = par[2] * sqrt(gamma(1+2/par[1]) - (gamma(1+1/par[1]))^2)
    return((mu_calc-1)^2 + ((cv_calc-cv)/cv)^2)
  }
  weib = optim(c(2.7,1.1), minfun)
  yweibull = dweibull(xvec, shape=weib$par[1], scale=weib$par[2])
  
  # Truncated Normal
  minfun = function(par) {
    mu_calc = etruncnorm(a=0, mean=par[1], sd=par[2])
    sigma_calc = sqrt(vtruncnorm(a=0, mean=par[1], sd=par[2]))
    return((mu_calc-1)^2 + ((sigma_calc-cv)/cv)^2)
  }
  tnorm = optim(c(1,cv), minfun)
  ytnorm = dtruncnorm(xvec, a=0, mean=tnorm$par[1], sd=tnorm$par[2])
  
  # Triangular
  minfun = function(par) {
    mu_calc = (par[1]+par[2]+par[3])/3
    sigma_calc = sqrt( (par[1]^2 + par[2]^2 + par[3]^2 - par[1]*par[2] - par[1]*par[3] - par[2]*par[3])/18 )
    return((mu_calc-1)^2 + ((sigma_calc-cv)/cv)^2)
  }
  ineq_con = function(par) return(c(par[2]-par[3], par[3]-par[1]))
  triang = suppressMessages(auglag(x0=c(0.2, 2.2, 0.5), fn=minfun, lower=c(0,0,0), hin=ineq_con, localsolver='LBFGS'))
  ytriang = dtriangle(xvec, a=triang$par[1], b=triang$par[2], c=triang$par[3])
  
  # Bootstrap Colors: teal, purple, orange, blue, red, yellow, cyan, pink, indigo, green
  cols = c('#20c997','#6f42c1','#fd7e14','#007bff','#dc3545','#ffc107','#17a2b8','#e83e8c','#6610f2','#28a745')
  par(mar=c(2.25,4,2,1), oma=c(0,0,0,0), fig=c(0,1,0,1))
  plot(xvec, ylnorm, xlab='', ylab='', type='n', xlim=c(0,3), ylim=c(0, max(ylnorm)))
  rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col='grey90',border=NA); grid(col='white',lty=1); box()
  mtext('Density', side=2, line=2.5)
  mtext('Multiplicative Error Distributions with mean=1', side=3, line=0.75)
  lines(xvec, ylnorm, col=cols[1], lwd=2)
  lines(xvec, ygamma, col=cols[2], lwd=2, lty=2)
  lines(xvec, yweibull, col=cols[3], lwd=2)
  lines(xvec, ytnorm, col=cols[4], lwd=2, lty=3)
  lines(xvec, ytriang, col=cols[5], lwd=1)
  legend('topright', legend=c('Lognormal','Gamma','Weibull','TruncNorm','Triangle'), title=paste0("CV = ",100*cv,"%"),
         col=cols[1:5], lty=c(1,2,1,3,1), lwd=c(2,2,2,2,1), cex=0.9)
}