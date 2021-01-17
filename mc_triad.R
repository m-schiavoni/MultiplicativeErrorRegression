# This function compares different regression methods for multiplicative error models.
# INPUTS:
#   n = sample size for each iteration
#   err = error distribution ('lognormal', 'gamma', 'weibull', 'truncnorm', 'triangle')
#   cv = Coefficient of Variation of random error
#   N = number of iterations. N=1 yields a scatterplot; N>1 yields Monte Carlo plots
#   seed = optional integer to set seed for random number generator
#
source("./sim_errors.R")
source("./plot_triad.R")
library(minpack.lm)
library(nloptr)
library(stats4)
mc_triad = function(n=15, err='gamma', cv=0.6, N=1, seed=NA) {
  
  if (N == 1) {
    scatter = TRUE
  } else {
    scatter = FALSE
    if (N < 500) {
      N = 500
      print('N increased to 500 as minimum value for Monte Carlo results.')
    }
  }
  
  # generate simulated data ##############################################
  b_0 = 30000;  b_1 = 25;  b_2 = 1.2          # true population parameters
  start_coeffs = c(b0=10000, b1=10, b2=1)     # hard-coded starting guess for triad model
  p = length(start_coeffs)                    # number of parameters
  
  if (is.na(seed)) seed = round(as.numeric(Sys.time()) - as.numeric(as.Date('2020-01-01'))*24*60*60)  # seconds since Jan 1, 2020
  set.seed(seed);  print(noquote(paste0('seed=', seed)))
  X = matrix(rweibull(n*N, shape=1.4, scale=1200), nrow=n, ncol=N)
  b0_vec = rep(b_0, N);  B0 = matrix(rep(b0_vec, each=n), nrow=n, ncol=N)
  b1_vec = rep(b_1, N);  B1 = matrix(rep(b1_vec, each=n), nrow=n, ncol=N)
  b2_vec = rep(b_2, N);  B2 = matrix(rep(b2_vec, each=n), nrow=n, ncol=N)
  Y = (B0 + B1*X^B2) * sim_errors(err, cv, N, n)
  
  # initialize outputs ##############################################
  my_df = as.data.frame(matrix(nrow=N, ncol=p))
  colnames(my_df) = c('b0','b1','b2','b3','b4','b5')[1:p]
  my_list = list('LogErr'=my_df, 'PING'=my_df, 'MUPE'=my_df, 'ZMPE'=my_df, 'GRMLN'=my_df, 
                 'n'=n, 'err'=err, 'cv'=cv, 'N'=N, 'seed'=seed, 'b0'=b_0, 'b1'=b_1, 'b2'=b_2)
  
  # define model form #####################################################
  mupe_form = as.formula("y ~ b0 + b1*x^b2")
  y_hat_func = function(par, x) return(par[1] + par[2]*x^par[3])
  
  # functions to be minimized ############################################
  min_logErr = function(par, x){
    y_hat = y_hat_func(par, x)
    return(log(y) - log(y_hat))  # nls.lm minimizes the sum square of this vector
  }
  min_grmln = function(par, x, p) {
    y_hat = y_hat_func(par[1:p], x)
    th = par[p+1]
    return(n*log(th)/2 + 1/(2*th)*sum((log(y) - log(y_hat) + th/2)^2))
  }
  min_zmpe = function(par, x) {
    y_hat = y_hat_func(par, x)
    return(sum(((y_hat-y)/y_hat)^2))  # sum of squared percent errors
  }
  ineq_con = function(par, x) {
    y_hat = y_hat_func(par, x)
    con = sum((y_hat-y)/y_hat)  # sum of percent errors
    return(c(con, -con))  # express equality constraint as pair of inequalities for COBYLA
  }
  
  # function to calculate squared standard percent error
  calc_spe2 = function(y, y_hat, n, p, c=0) return(sum(((y-y_hat)/y_hat)^2)/(n-p-c))
  
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
      # Bootstrap Colors: yellow, purple, orange, blue, teal, green, red, cyan, pink, indigo
      cols = c('#ffc107','#6f42c1','#fd7e14','#007bff','#20c997','#28a745','#dc3545','#17a2b8','#e83e8c','#6610f2')
      xmin = min(X);  xmax=max(X)
      par(mar=c(2,2,2.5,1))
      plot(x,y,xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(xmin,xmax),ylim=c(min(y),max(y)*1.1),type='n')
      rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4], col='grey90', border=NA)
      grid(col='white',lty=1);  box();  points(x,y,pch=1)
      mtext('Driver',side=1,line=0.6);  mtext('Response',side=2,line=0.6)
      mtext('Multiplicative Error Method Comparison',side=3,line=0.75,cex=1.05)
      xvec = seq(xmin,xmax,length.out=500)
      legend('topleft', lwd=c(2,2,2,2,2,1), lty=c(1,3,1,4,2), cex=0.8,
             col=c(cols[1:5],'red'),
             legend=c('LogErr','PING','GRMLN','MUPE','ZMPE','truth'), ncol=1)
    }
    
    ## log error ######################################
    logErr = try(suppressWarnings(nls.lm(par=start_coeffs, lower=rep(0, p), fn=min_logErr, x=x, 
                                         control=list(maxiter=200))), silent = TRUE)
    if (class(logErr) != 'try-error') {
      if (logErr$info != 0 & logErr$info != 5 & logErr$info != 9) {
        y_hat = y_hat_func(logErr$par, x)
        var_hat = calc_spe2(y, y_hat, n, p)
        my_list$LogErr[i,] = logErr$par
        
        ## PING correction to LogErr ##################################
        ping = exp((1 - p/n)*var_hat/2)
        ping_sol = c(ping*logErr$par[1], ping*logErr$par[2], logErr$par[3])
        my_list$PING[i,] = ping_sol
        if (scatter) {
          yvec = y_hat_func(logErr$par, xvec)
          lines(xvec, yvec, col=cols[1], lwd=2, lty=1)
          yvec = y_hat_func(ping_sol, xvec)
          lines(xvec, yvec, col=cols[2], lwd=2, lty=3)
        }
      }
    }
    
    ## GRMLN ######################################
    grmln = try(suppressWarnings(optim(par=c(start_coeffs, 0.3), fn=min_grmln, x=x, p=p, 
                                       method='BFGS', control=list(maxit=200))), silent=TRUE)
    if (class(grmln) != 'try-error') {
      if (grmln$convergence == 0) {
        my_list$GRMLN[i,] = grmln$par[1:p]
        if (scatter) {
          yvec = y_hat_func(grmln$par[1:p], xvec)
          lines(xvec, yvec, col=cols[3], lwd=2, lty=1)
        }
      }
    }
    
    ## MUPE #######################################################
    pbeta = start_coeffs;  wt = rep(1, n);  conv = 1.0;  outer = 0
    while (conv > 1e-5) {
      mupe = try(suppressWarnings(nlsLM(formula=mupe_form, start=pbeta, weights=wt,
                                        control=list(maxiter=10))), silent=TRUE)
      if (class(mupe) == 'try-error') break
      beta = mupe$m$getAllPars()          # solution of current iteration
      conv = max(abs((beta-pbeta)/beta))  # maximum fractional change in any one parameter
      wt = 1 / mupe$m$pred()^2            # reset weights
      pbeta = beta                        # reset prior beta
      outer = outer + 1;  if (outer == 200) break
    }
    if (class(mupe) != 'try-error' & outer < 200) {
      my_list$MUPE[i,] = beta
      if (scatter) {
        yvec = y_hat_func(beta, xvec)
        lines(xvec, yvec, col=cols[4], lwd=2, lty=4)
      }
    }
    
    ## ZMPE #####################################################
    zmpe = try(suppressMessages(cobyla(x0=start_coeffs, fn=min_zmpe, hin=ineq_con, x=x, 
                                       control=list(maxeval=5000))))
    if (class(zmpe) != 'try-error') {
      if (zmpe$convergence > 0) {
        my_list$ZMPE[i,] = zmpe$par
        if (scatter) {
          yvec = y_hat_func(zmpe$par, xvec)
          lines(xvec, yvec, col=cols[5], lwd=2, lty=2)
        }
      }
    }
  }
  
  if (scatter) {  # overlay true population curve
    yvec = y_hat_func(c(b_0, b_1, b_2), xvec)
    lines(xvec, yvec, col='white', lwd=2)
    lines(xvec, yvec, col='red')
    box()
  } else {  # plot Monte Carlo sim results
    print('100%', quote=FALSE)
    plot_triad(my_list)
  }
  return(my_list)
}
