# This function compares different regression methods for multiplicative error models.
# INPUTS:
#   n = sample size for each iteration
#   err = error distribution ('lognormal', 'gamma', 'weibull', 'truncnorm', 'triangle')
#   cv = Coefficient of Variation of random error
#   N = number of iterations. N=1 yields a scatterplot; N>1 yields Monte Carlo plots
#   seed = optional integer to set seed for random number generator
#
source("./sim_errors.R")
source("./plot_multivariate.R")
library(minpack.lm)
library(nloptr)
library(stats4)
library(MASS)  # using mvrnorm()
mc_multivariate = function(n=30, err='gamma', cv=0.6, N=1, seed=NA) {
  
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
  b_0 = 35;  b_1 = 0.75;  b_2 = 0.85;  b_3 = 1.3    # true population parameters
  p = 4                                             # number of parameters
  rho_target = 0.25    # target correlation between x1 and x2 (will vary somewhat from sample to sample)
  
  if (is.na(seed)) seed = round(as.numeric(Sys.time()) - as.numeric(as.Date('2020-01-01'))*24*60*60)  # seconds since Jan 1, 2020
  set.seed(seed);  print(noquote(paste0('seed=', seed)))
  
  Sigma = matrix(rho_target, nrow=2, ncol=2) + diag(2)*(1-rho_target)  # covariance matrix
  x3_mat = matrix(rbinom(n=N*n, size=1, prob=0.5), nrow=n, ncol=N)  # binomial dummy variable
  err_mat = sim_errors(err, cv, N, n)
  
  # initialize outputs ##############################################
  my_df = as.data.frame(matrix(nrow=N, ncol=p))
  colnames(my_df) = c('b0','b1','b2','b3','b4','b5')[1:p]
  my_list = list('LogErr'=my_df, 'PING'=my_df, 'MUPE'=my_df, 'ZMPE'=my_df, 'GRMLN'=my_df, 
                 'n'=n, 'err'=err, 'cv'=cv, 'N'=N, 'seed'=seed, 'b0'=b_0, 'b1'=b_1, 'b2'=b_2, 'b3'=b_3,
                 'sample_corr'=vector('numeric',n))
  
  # define model form #####################################################
  mupe_form = as.formula("y ~ b0 * x1^b1 * x2^b2 * b3^x3")
  y_hat_func = function(par, X) return(par[1] * X[,1]^par[2] * X[,2]^par[3] * par[4]^X[,3])
  
  # functions to be minimized ############################################
  min_grmln = function(par, X, p) {
    y_hat = y_hat_func(par[1:p], X)
    th = par[p+1]
    return(n*log(th)/2 + 1/(2*th)*sum((log(y) - log(y_hat) + th/2)^2))
  }
  min_zmpe = function(par, X) {
    y_hat = y_hat_func(par, X)
    return(sum(((y_hat-y)/y_hat)^2))  # sum of squared percent errors
  }
  ineq_con = function(par, X) {
    y_hat = y_hat_func(par, X)
    con = sum((y_hat-y)/y_hat)  # sum of percent errors
    return(c(con, -con))  # express equality constraint as pair of inequalities for COBYLA
  }
  
  # function to calculate squared standard percent error
  calc_spe2 = function(y, y_hat, n, p, c=0) return(sum(((y-y_hat)/y_hat)^2)/(n-p-c))
  
  # Actuals vs. Fitted plot
  act_vs_fit = function(y, y_hat, leg) {
    lims = c(0, max(y, y_hat))
    plot(y~y_hat, xlab='', ylab='', xlim=lims, ylim=lims)
    mtext('Fitted Values', side=1, line=2.5, cex=0.9)
    mtext('Actual Values', side=2, line=2.5, cex=0.9)
    legend('topleft', leg)
  }
  
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
    
    # simulate correlated x1, x2
    # do NOT relocate this code outside the loop; doing so will increase variation of correlation around the target
    bivar_norm = mvrnorm(n=n, mu=c(0,0), Sigma=Sigma, empirical=TRUE)  # bivariate normal random variables
    probs = pnorm(bivar_norm)                                          # distribution function probabilities
    myweib = qweibull(probs, shape=1.4, scale=1200)                    # bivariate Weibull random variables
    x1 = myweib[,1];  x2 = myweib[,2];  x3 = x3_mat[,i]
    X_mat = matrix(c(x1,x2,x3), nrow=n, ncol=3)
    my_list$sample_corr[i] = cor(x1,x2)
    
    # simulate response values
    y = (b_0 * x1^b_1 * x2^b_2 * b_3^x3) * err_mat[,i]
    
    # scatterplot ##########################################
    if (scatter) {
      # Bootstrap Colors: yellow, purple, orange, blue, teal, green, red, cyan, pink, indigo
      cols = c('#ffc107','#6f42c1','#fd7e14','#007bff','#20c997','#28a745','#dc3545','#17a2b8','#e83e8c','#6610f2')
      pairs(~y+x1+x2+x3, lower.panel=NULL)
      par(mfrow=c(2,2), oma=c(0,0,2.5,0), mar=c(4.5,4,0.5,1))
    }
    
    ## log error ######################################
    logy = log(y);  logx1 = log(x1);  logx2 = log(x2)
    lols = lm(logy ~ logx1 + logx2 + x3)
    b0 = as.numeric(exp(lols$coef[1]))
    b1 = as.numeric(lols$coef[2])
    b2 = as.numeric(lols$coef[3])
    b3 = as.numeric(exp(lols$coef[4]))
    my_list$LogErr[i,] = c(b0, b1, b2, b3)
    
    ## PING correction to LogErr ##################################
    y_hat = y_hat_func(c(b0,b1,b2,b3), X_mat)
    var_hat = calc_spe2(y, y_hat, n, p)
    b0 = b0 * exp((1 - p/n)*var_hat/2)
    my_list$PING[i,] = c(b0, b1, b2, b3)
    if (scatter) {
      y_hat = y_hat_func(c(b0,b1,b2,b3), X_mat)
      act_vs_fit(y, y_hat, 'PING')
      abline(a=0, b=1, col=cols[2], lwd=2, lty=3)
    }
    
    # use PING solution as starting guess for other methods, because it is a 
    # direct-solve method (no optimization needed) that is approximately unbiased
    start_coeffs = c('b0'=b0, 'b1'=b1, 'b2'=b2, 'b3'=b3)
    y_hat = y_hat_func(c(b0,b1,b2,b3), X_mat)
    var_hat = calc_spe2(y, y_hat, n, p)
    
    ## GRMLN ######################################
    grmln = try(suppressWarnings(optim(par=c(start_coeffs, var_hat), fn=min_grmln, X=X_mat, p=p,
                                       method='BFGS', control=list(maxit=200))), silent=TRUE)
    if (class(grmln) != 'try-error') {
      if (grmln$convergence == 0) {
        my_list$GRMLN[i,] = grmln$par[1:p]
        if (scatter) {
          y_hat = y_hat_func(grmln$par[1:p], X_mat)
          act_vs_fit(y, y_hat, 'GRMLN')
          abline(a=0, b=1, col=cols[3], lwd=2, lty=1)
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
        y_hat = y_hat_func(beta, X_mat)
        act_vs_fit(y, y_hat, 'MUPE')
        abline(a=0, b=1, col=cols[4], lwd=2, lty=4)
      }
    }
    
    ## ZMPE #####################################################
    zmpe = try(suppressMessages(cobyla(x0=start_coeffs, fn=min_zmpe, hin=ineq_con, X=X_mat,
                                       control=list(maxeval=5000))))
    if (class(zmpe) != 'try-error') {
      if (zmpe$convergence > 0) {
        my_list$ZMPE[i,] = zmpe$par
        if (scatter) {
          y_hat = y_hat_func(zmpe$par, X_mat)
          act_vs_fit(y, y_hat, 'ZMPE')
          abline(a=0, b=1, col=cols[5], lwd=2, lty=2)
        }
      }
    }
  }
  
  if (scatter) {  # add title to scatterplot
    par(fig=c(0,1,0,1), oma=c(0,0,1.5,0), mar=c(0,0,0,0), new=TRUE)
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    mtext('Multiplicative Error Method Comparison',side=3,line=0,cex=1)
    par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(4,4,4,1))
  } else {  # plot Monte Carlo sim results
    print('100%', quote=FALSE)
    plot_multivariate(my_list)
  }
  return(my_list)
}
