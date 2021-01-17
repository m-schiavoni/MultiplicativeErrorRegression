# function to plot Monte Carlo simulation results
plot_linear = function(list_in) {
  N = list_in$N
  if (list_in$N == 1) stop('Input cannot be a single iteration.')
  twos = 2^(7:17)
  ivec = c(twos[twos<N], N)
  mylen = length(ivec)
  
  # Bootstrap Colors: yellow, purple, orange, blue, teal, green, red, cyan, pink, indigo
  cols = c('#ffc107','#6f42c1','#fd7e14','#007bff','#20c997','#28a745','#dc3545','#17a2b8','#e83e8c','#6610f2')
  
  b0 = data.frame(LogErr=vector('numeric',mylen), PING=vector('numeric',mylen), GRMLN=vector('numeric',mylen), 
                  MUPE=vector('numeric',mylen), ZMPE=vector('numeric',mylen), row.names=ivec)
  b1 = b0
  for (i in 1:mylen) {
    b0$LogErr[i] = median(list_in$LogErr$b0[1:ivec[i]], na.rm=TRUE)
    b0$PING[i] = median(list_in$PING$b0[1:ivec[i]], na.rm=TRUE)
    b0$GRMLN[i] = median(list_in$GRMLN$b0[1:ivec[i]], na.rm=TRUE)
    b0$MUPE[i] = median(list_in$MUPE$b0[1:ivec[i]], na.rm=TRUE)
    b0$ZMPE[i] = median(list_in$ZMPE$b0[1:ivec[i]], na.rm=TRUE)
    b1$LogErr[i] = median(list_in$LogErr$b1[1:ivec[i]], na.rm=TRUE)
    b1$PING[i] = median(list_in$PING$b1[1:ivec[i]], na.rm=TRUE)
    b1$GRMLN[i] = median(list_in$GRMLN$b1[1:ivec[i]], na.rm=TRUE)
    b1$MUPE[i] = median(list_in$MUPE$b1[1:ivec[i]], na.rm=TRUE)
    b1$ZMPE[i] = median(list_in$ZMPE$b1[1:ivec[i]], na.rm=TRUE)
  }
  
  par(mfrow=c(2,2), oma=c(0,0,3,0), mar=c(4.5,4,0.5,1))
  
  # b1 convergence plot
  ylims = c(min(0.97*list_in$b1, min(b1)), max(1.03*list_in$b1, max(b1)))
  matplot(b1, type='n', ylim=ylims, xaxt='n', xlab='', ylab='', cex.axis=0.9, cex.lab=0.9, new=TRUE)
  rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col='grey90',border=NA); grid(col='white',lty=1); box()
  abline(h=list_in$b1, col='red2')
  matlines(b1, xaxt='n', yaxt='n', col=cols[1:5], lwd=2, lty=c(1,3,1,4,2), new=TRUE)
  axis(side=1, at=1:mylen, labels=ivec, las=3, cex.axis=0.9)
  mtext(side=2, line=2.5, 'Median of a Estimates', cex=0.9)
  mtext(side=1, line=3, '# Monte Carlo Iterations', cex=0.9)
  
  # b1 boxplots
  df = data.frame(LogErr=list_in$LogErr$b1, PING=list_in$PING$b1, GRMLN=list_in$GRMLN$b1, 
                  MUPE=list_in$MUPE$b1, ZMPE=list_in$ZMPE$b1)
  boxplot(df, col=cols, ylab='', las=3, cex.axis=0.9, cex.lab=0.9, outline=FALSE)
  rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col='grey95',border=NA); box()
  abline(h=list_in$b1, col='red2')
  boxplot(df, col=cols, xlab='', ylab='', xaxt='n', yaxt='n', outline=FALSE, add=TRUE)
  mtext(side=2, line=2.5, 'a Estimates', cex=0.9)
  
  # b0 convergence plot
  ylims = c(min(0.97*list_in$b0, min(b0)), max(1.03*list_in$b0, max(b0)))
  matplot(b0, type='n', ylim=ylims, xaxt='n', xlab='', ylab='', cex.axis=0.9, cex.lab=0.9, new=TRUE)
  rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col='grey90',border=NA); grid(col='white',lty=1); box()
  abline(h=list_in$b0, col='red2')
  matlines(b0, xaxt='n', yaxt='n', col=cols[1:5], lwd=2, lty=c(1,3,1,4,2), new=TRUE)
  axis(side=1, at=1:mylen, labels=ivec, las=3, cex.axis=0.9)
  mtext(side=2, line=2.5, 'Median of c Estimates', cex=0.9)
  mtext(side=1, line=3, '# Monte Carlo Iterations', cex=0.9)
  
  # b0 boxplots
  df = data.frame(LogErr=list_in$LogErr$b0, PING=list_in$PING$b0, GRMLN=list_in$GRMLN$b0, 
                  MUPE=list_in$MUPE$b0, ZMPE=list_in$ZMPE$b0)
  boxplot(df, col=cols, ylab='', las=3, cex.axis=0.9, cex.lab=0.9, outline=FALSE)
  rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col='grey95',border=NA); box()
  abline(h=list_in$b0, col='red2')
  boxplot(df, col=cols, xlab='', ylab='', xaxt='n', yaxt='n', outline=FALSE, add=TRUE)
  mtext(side=2, line=2.5, 'c Estimates', cex=0.9)
  
  # title and legend
  par(fig=c(0,1,0,1), oma=c(0,0,1,0), mar=c(0,0,0,0), new=TRUE)
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  mtext(side=3, line=0, paste0('Monte Carlo Simulation Plots (y = a*x + c) [n=', list_in$n, 
                               ';  ', 'cv=', list_in$cv, ';  ', 'err=', list_in$err, ']'), cex=1.0)
  myleg = c('LogErr', 'PING', 'GRMLN', 
            'MUPE', 'ZMPE')
  legend('top', legend=myleg, col=cols[1:5], lwd=3, lty=c(1,3,1,4,2),
         xpd=TRUE, horiz=TRUE, cex=1.0, seg.len=2.5, bty='n', text.width=1.2*strwidth(myleg))
  par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(4.5,4.5,2.5,1))  # reset figure parameters
}
