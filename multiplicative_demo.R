# This function compares different regression methods for multiplicative error models.
# INPUTS:
#   n = sample size for each iteration
#   err = error distribution ('lognormal', 'gamma', 'weibull', 'truncnorm', 'triangle')
#   cv = Coefficient of Variation of random error
#   N = number of iterations. N=1 yields a scatterplot; N>1 yields Monte Carlo plots
#   fixed = if TRUE, population parameters are same for each iteration
#           if FALSE, population parameters are selected randomly each iteration
#   seed = optional integer to set seed for random number generator
#
multiplicative_demo = function(n=30, err='lognormal', cv=0.6, N=1, fixed=FALSE, seed=NA) {
  require(minpack.lm)
  require(truncnorm)
  require(triangle)
  require(nloptr)
  require(stats4)
  
  # initialization ###################################
  if (N == 1) {
    scatter = TRUE
  } else {
    scatter = FALSE
    if (N < 200) {
      N = 200
      print('N increased to 200 as minimum value for Monte Carlo results.')
    }
  }
  if (is.na(seed)) seed = round(as.numeric(Sys.time()) - as.numeric(as.Date('2020-01-01'))*24*60*60)  # seconds since Jan 1, 2020
  set.seed(seed);  print(noquote(paste0('seed=', seed)))
  my_df = as.data.frame(matrix(nrow=N, ncol=3))
  colnames(my_df) = c('var_hat','b0','b1')
  my_list = list('n'=n, 'err'=err, 'cv'=cv, 'N'=N, 'fixed'=fixed, 'seed'=seed, 
                 'MPE'=my_df, 'LOLS'=my_df,'Gold'=my_df,  'PING'=my_df, 
                 'ZMPE'=my_df, 'MUPE'=my_df, 'MRLN'=my_df, 'b0'=0, 'b1'=0)
  
  # generate simulated data #####################################################
  xmin = 1000;  xmax = 8000;  X = matrix(runif(n*N, xmin, xmax), nrow=n, ncol=N)
  if (fixed) {  # same parameter values for each iteration
    b0 = 75;  my_list$b0 = b0
    b1 = 0.8;  my_list$b1 = b1
    b0_vec = rep(b0, N)
    b1_vec = rep(b1, N)
  } else {  # random parameter values for each iteration
    b0_vec = runif(N, 20, 200)
    b1_vec = runif(N, 0.6, 1.5)
  }
  B0 = matrix(rep(b0_vec, each=n), nrow=n, ncol=N)
  B1 = matrix(rep(b1_vec, each=n), nrow=n, ncol=N)
  Y = B0 * X^B1 * sim_errors(err, cv, N, n)
  p = 2  # number of parameters
  
  # run simulation ########################################
  timer = 0
  mult = 0
  for (i in 1:N) {
    if (N > 1) {
      if (i > mult*N/10) {
        print(paste0(timer,'%'), quote=FALSE)
        timer = timer + 10
        mult = mult + 1
      }
    }
    
    # subset data for this iteration
    x = X[,i];  y = Y[,i]
    
    # scatterplot ##########################################
    if (scatter) {
      # Bootstrap Colors: yellow, purple, teal, blue, orange, green, red, cyan, pink, indigo
      cols = c('#ffc107','#6f42c1','#20c997','#007bff','#fd7e14','#28a745','#dc3545','#17a2b8','#e83e8c','#6610f2')
      par(mar=c(2,2,2.5,1))
      plot(x,y,xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(xmin,xmax),ylim=c(min(y),max(y)*1.1),type='n')
      rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4], col='grey90', border=NA)
      grid(col='white',lty=1);  box();  points(x,y,pch=1)
      mtext('Driver',side=1,line=0.6);  mtext('Response',side=2,line=0.6)
      mtext('Multiplicative Error Method Comparison',side=3,line=0.75,cex=1.05)
      xvec = seq(xmin,xmax,length.out=500)
      legend('topleft', lwd=c(1,1,2,2,2,2,2,1), lty=c(2,3,1,3,1,4,2,1), cex=0.8,
             col=c('black','black',cols[1:5],'red'),
             legend=c('MPE','LOLS','Gold','PING','ZMPE','MUPE','MRLN','truth'), ncol=1)
    }
    
    ## log error ######################################
    logx = log(x);  logy = log(y)
    lm1 = lm(logy ~ logx)
    b0 = as.numeric(exp(lm1$coef[1]))
    b1 = as.numeric(lm1$coef[2]) 
    y_hat = b0*x^b1
    var_hat0 = calc_spe2(y, y_hat, n, p)
    my_list$LOLS[i,] = c(var_hat0, b0, b1)
    if (scatter) {yvec = b0*xvec^b1;  lines(xvec, yvec, lty=3)}
    
    ## Goldberger correction to LOLS ##########################
    meanlogx = mean(logx)
    y_hat = b0*x^b1 * exp((1 - 1/n - ((logx - meanlogx)^2) / sum((logx-meanlogx)^2))*var_hat0/2)
    var_hat = calc_spe2(y, y_hat, n, p)
    my_list$Gold[i,] = c(var_hat, NA, NA)
    if (scatter) {
      logxv = log(xvec)
      yvec = b0*xvec^b1 * exp((1 - 1/n - ((logxv - meanlogx)^2) / sum((logxv-meanlogx)^2))*var_hat0/2)
      lines(xvec, yvec, col=cols[1], lwd=2, lty=1)
    }
    
    ## PING correction to LOLS ##################################
    b0 = b0 * exp((1 - p/n)*var_hat0/2)
    y_hat = b0*x^b1
    var_hat = calc_spe2(y, y_hat, n, p)
    my_list$PING[i,] = c(var_hat, b0, b1)
    if (scatter) {yvec = b0*xvec^b1;  lines(xvec, yvec, col=cols[2], lwd=2, lty=3)}
    
    # set starting point for nonlinear optimization techniques
    start = c(b0=b0, b1=b1);  v_start = var_hat  # PING starting point
    # start = c(b0=10, b1=1);  v_start = 0.25  # dummy starting point
    
    
    ## MPE ######################################
    min_mpe = function(par) {  # this minimization function is only valid for y=b0*x^b1
      y_hat = par[1]*x^par[2]  # predicted values
      return((y_hat-y)/y_hat)  # nls.lm() minimizes the sum square of this vector
    }
    mpe = try(suppressWarnings(nls.lm(par=start, fn=min_mpe)), silent=TRUE)
    if (class(mpe) != 'try-error') {
      b0 = as.numeric(mpe$par[1]);  b1 = as.numeric(mpe$par[2])
      y_hat = b0*x^b1
      var_hat = calc_spe2(y, y_hat, n, p)
      my_list$MPE[i,] = c(var_hat, b0, b1)
    }
    if (scatter) {yvec = b0*xvec^b1;  lines(xvec, yvec, lty=2)}
    
    
    ## ZMPE ################################
    min_zmpe = function(par) {          # this minimization function is only valid for y=b0*x^b1
      y_hat = par[1] * x ^ par[2]       # predicted values
      return(sum(((y_hat-y)/y_hat)^2))  # sum of squared percent errors
    }
    # express equality constraint as pair of inequalities for COBYLA algorithm
    ineq_con = function(par) {     # this constraint is only valid for y=b0*x^b1
      y_hat = par[1] * x ^ par[2]  # predicted values
      con = sum((y_hat-y)/y_hat)   # sum of percent errors
      return(c(con, -con))
    }
    zmpe = try(suppressMessages(cobyla(x0=start, fn=min_zmpe, hin=ineq_con, control=list(maxeval=2000))))
    if (class(zmpe) != 'try-error') {
      b0 = as.numeric(zmpe$par[1]);  b1 = as.numeric(zmpe$par[2])
      y_hat = b0*x^b1
      var_hat = calc_spe2(y, y_hat, n, p, c=0)  # 1 constraint
      my_list$ZMPE[i,] = c(var_hat, b0, b1)
    }
    if (scatter) {yvec = b0*xvec^b1;  lines(xvec, yvec, col=cols[3], lwd=2, lty=1)}
    
    
    ## MUPE #######################################################
    pbeta = start;  wt = rep(1, n);  conv = 1.0;  outer = 0
    while (conv > 1e-5) {
      model = try(suppressWarnings(nlsLM(formula=y ~ b0 * x^b1, start=pbeta, weights=wt, 
                                         control=list(maxiter=10))), silent=TRUE)
      if (class(model) == 'try-error') break
      beta = model$m$getAllPars()         # solution of current iteration
      conv = max(abs((beta-pbeta)/beta))  # maximum fractional change in any one parameter
      wt = 1 / model$m$pred()^2           # reset weights
      pbeta = beta                        # reset prior beta
      outer = outer + 1;  if (outer > 200) break
    }
    if (class(model) != 'try-error') {
      b0 = beta[1];  b1 = beta[2]
      y_hat = b0*x^b1
      var_hat = calc_spe2(y, y_hat, n, p)
      my_list$MUPE[i,] = c(var_hat, b0, b1)
    }
    if (scatter) {yvec = b0*xvec^b1;  lines(xvec, yvec, col=cols[4], lwd=2, lty=4)}
    
    
    ## MRLN ######################################
    min_mrln = function(par) {  # this minimization function is only valid for y=b0*x^b1
      b0 = par[1];  b1 = par[2];  th = par[3]
      return(n*log(th)/2 + 1/(2*th)*sum((log(y) - log(b0) - b1*log(x) + th/2)^2))
    }
    mrln = try(suppressWarnings(optim(par=c(start, v_start), 
                                      fn=min_mrln, method='BFGS')), silent=TRUE)
    if (class(mrln) != 'try-error') {
      b0 = as.numeric(mrln$par[1]);  b1 = as.numeric(mrln$par[2])
      y_hat = b0*x^b1
      var_hat = calc_spe2(y, y_hat, n, p)
      my_list$MRLN[i,] = c(var_hat, b0, b1)
    }
    if (scatter) {yvec = b0*xvec^b1;  lines(xvec, yvec, col=cols[5], lwd=2, lty=2)}
  }
  
  
  if (scatter) {  # overlay true population curve
    yvec = b0_vec * xvec^b1_vec
    lines(xvec, yvec, col='white', lwd=2)
    lines(xvec, yvec, col='red')
    box()
  } else {  # plot Monte Carlo sim results
    print('100%', quote=FALSE)
    plot_mc(my_list, 'mean')
    # plot_mc(my_list, 'median')
  }
  return(my_list)
}


# function to calculate squared standard percent error
calc_spe2 = function(y, y_hat, n, p, c=0) {return(sum(((y-y_hat)/y_hat)^2)/(n-p-c))}


# function to simulate multiplicative errors with mean=1 and specified CV
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
    minfun = function(par) {
      mu_calc = (par[1]+par[2]+par[3])/3
      var_calc = (par[1]^2 + par[2]^2 + par[3]^2 - par[1]*par[2] - par[1]*par[3] - par[2]*par[3]) / 18
      return((mu_calc-1)^2 + ((var_calc-cv^2)/cv^2)^2)
    }
    ineq_con = function(par) return(c(par[2]-par[3], par[3]-par[1]))
    triang = suppressMessages(auglag(x0=c(0.2, 2.2, 0.5), fn=minfun, lower=c(0,0,0), hin=ineq_con, localsolver='LBFGS'))
    rand_error = matrix(rtriangle(N*n, a=triang$par[1], b=triang$par[2], c=triang$par[3]), nrow=n, ncol=N)
    
  } else {
    stop('ERROR: invalid error model.')
  }
  
  return(rand_error)
}


# function to plot Monte Carlo simulation results
plot_mc = function(list_in, func='mean') {
  if (list_in$N == 1) stop('Input cannot be a single iteration.')
  if (func == 'mean') {
    lab = 'Mean '
    central_tendency = function(x) mean(x, na.rm=TRUE)
  } else {
    lab = 'Median '
    central_tendency = function(x) median(x, na.rm=TRUE)
  }
  N = list_in$N
  twos = 2^(7:17)
  ivec = c(twos[twos<N], N)
  mylen = length(ivec)
  
  # Bootstrap Colors: yellow, purple, teal, blue, orange, green, red, cyan, pink, indigo
  cols = c('#ffc107','#6f42c1','#20c997','#007bff','#fd7e14','#28a745','#dc3545','#17a2b8','#e83e8c','#6610f2')
  
  if (list_in$fixed) {
    b0 = data.frame(PING=vector('numeric',mylen), ZMPE=vector('numeric',mylen), 
                    MUPE=vector('numeric',mylen), MRLN=vector('numeric',mylen), 
                    row.names=ivec)
    b1 = b0
    for (i in 1:mylen) {
      b0$PING[i] = central_tendency(list_in$PING$b0[1:ivec[i]])
      b0$ZMPE[i] = central_tendency(list_in$ZMPE$b0[1:ivec[i]])
      b0$MUPE[i] = central_tendency(list_in$MUPE$b0[1:ivec[i]])
      b0$MRLN[i] = central_tendency(list_in$MRLN$b0[1:ivec[i]])
      b1$PING[i] = central_tendency(list_in$PING$b1[1:ivec[i]])
      b1$ZMPE[i] = central_tendency(list_in$ZMPE$b1[1:ivec[i]])
      b1$MUPE[i] = central_tendency(list_in$MUPE$b1[1:ivec[i]])
      b1$MRLN[i] = central_tendency(list_in$MRLN$b1[1:ivec[i]])
    }
    
    par(mfrow=c(1,2), oma=c(0,0,3,0), mar=c(4,4,0,1))
    matplot(b0, type='n', xaxt='n', xlab='', ylab=paste(lab, 'of b0 Estimates'), cex.axis=0.9, cex.lab=0.9, 
            ylim=c(0, max(list_in$b0, max(b0, na.rm=TRUE))), las=1, new=TRUE)
    rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col='grey90',border=NA); grid(col='white',lty=1); box()
    abline(h=list_in$b0, col='red2')
    matlines(b0, xaxt='n', yaxt='n', col=cols[2:5], lwd=2, lty=c(3,1,4,2), new=TRUE)
    axis(side=1, at=1:mylen, labels=ivec, las=2, cex.axis=0.9)
    mtext(side=1, line=3, '# Monte Carlo Iterations', cex=0.9)
    
    matplot(b1, type='n', xaxt='n', xlab='', ylab=paste(lab, 'of b1 Estimates'), cex.axis=0.9, cex.lab=0.9, 
            ylim=c(min(0.95*list_in$b1, min(b1, na.rm=TRUE)), max(1.05*list_in$b1, max(b1, na.rm=TRUE))), las=1, new=TRUE)
    rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col='grey90',border=NA); grid(col='white',lty=1); box()
    abline(h=list_in$b1, col='red2')
    matlines(b1, xaxt='n', yaxt='n', col=cols[2:5], lwd=2, lty=c(3,1,4,2), new=TRUE)
    axis(side=1, at=1:mylen, labels=ivec, las=2, cex.axis=0.9)
    mtext(side=1, line=3, '# Monte Carlo Iterations', cex=0.9)
    
    par(fig=c(0,1,0,1), oma=c(0,0,1,0), mar=c(0,0,0,0), new=TRUE)
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    mtext(side=3, line=0, paste0('Monte Carlo Simulation Convergence Plots [n=', list_in$n, 
                                 ';  ', 'cv=', list_in$cv, ';  ', 'err=', list_in$err, ']'), cex=1.0)
    myleg = c('PING', 'ZMPE', 'MUPE', 'MRLN')
    legend('top', legend=myleg, col=cols[2:5], lwd=3, lty=c(3,1,4,2), 
           xpd=TRUE, horiz=TRUE, cex=0.9, seg.len=2.5, bty='n', text.width=1.2*strwidth(myleg))
  } else {
    var_hat = data.frame(Gold=vector('numeric',mylen), PING=vector('numeric',mylen), 
                         ZMPE=vector('numeric',mylen), MUPE=vector('numeric',mylen), 
                         MRLN=vector('numeric',mylen), row.names=ivec)
    for (i in 1:mylen) {
      var_hat$Gold[i] = central_tendency(list_in$Gold$var_hat[1:ivec[i]])
      var_hat$PING[i] = central_tendency(list_in$PING$var_hat[1:ivec[i]])
      var_hat$ZMPE[i] = central_tendency(list_in$ZMPE$var_hat[1:ivec[i]])
      var_hat$MUPE[i] = central_tendency(list_in$MUPE$var_hat[1:ivec[i]])
      var_hat$MRLN[i] = central_tendency(list_in$MRLN$var_hat[1:ivec[i]])
    }
    par(mfrow=c(1,2), oma=c(0,0,3,0), mar=c(4,4,0,1))
    var = list_in$cv^2
    matplot(var_hat, type='n', xaxt='n', xlab='', ylab=paste(lab, 'of Variance Estimates'), cex.axis=0.9, cex.lab=0.9, 
            ylim=c(min(0.9*var, min(var_hat, na.rm=TRUE)), max(1.1*var, max(var_hat, na.rm=TRUE))), las=1, new=TRUE)
    rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col='grey90',border=NA); grid(col='white',lty=1); box()
    abline(h=var, col='red2')
    matlines(var_hat, xaxt='n', yaxt='n', col=cols, lwd=2, lty=c(1,3,1,4,2), new=TRUE)
    axis(side=1, at=1:mylen, labels=ivec, las=2, cex.axis=0.9)
    mtext(side=1, line=3, '# Monte Carlo Iterations', cex=0.9)
    
    df = data.frame(Gold=list_in$Gold$var_hat, PING=list_in$PING$var_hat, 
                    ZMPE=list_in$ZMPE$var_hat, MUPE=list_in$MUPE$var_hat, MRLN=list_in$MRLN$var_hat)
    boxplot(df, col=cols, las=2, ylab='Variance Estimates', cex.axis=0.9, cex.lab=0.9)
    rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col='grey95',border=NA); box()
    boxplot(df, col=cols, las=2, xlab='', ylab='', xaxt='n', yaxt='n', add=TRUE)
    
    par(fig=c(0,1,0,1), oma=c(0,0,1,0), mar=c(0,0,0,0), new=TRUE)
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    mtext(side=3, line=0, paste0('Monte Carlo Simulation Variance Plots [n=', list_in$n, 
                                 ';  ', 'cv=', list_in$cv, ';  ', 'err=', list_in$err, ']'), cex=1.0)
    myleg = c('Gold', 'PING', 'ZMPE', 'MUPE', 'MRLN')
    legend('top', legend=myleg, col=cols, lwd=3, lty=c(1,3,1,4,2), 
           xpd=TRUE, horiz=TRUE, cex=0.9, seg.len=2.5, bty='n', text.width=1.2*strwidth(myleg))
  }
  par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(4.5,4.5,2.5,1))
}
