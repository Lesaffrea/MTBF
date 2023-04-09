
# R.J. Marriott. 29 June 2016. (Version 2.0)

##################

# Calculate probability of failure:

Calc.Unreliability.w2p &lt;- function(beta,eta,time){
  Unreliability &lt;- (1-exp(-((time/eta)^beta)))
  return(Unreliability)
}

##################

# Calculate warranty time for target reliability:

Calc.Warranty.w2p &lt;- function(beta,eta,Rval){
  time=eta*((-log(Rval))^(1/beta))
  return(time)
}

##################

# NormalDist.rrx method:

xmax &lt;- function(x){
  # this function gets upper limit to x axis for the plot
  # to replicate the Weibull++ plot format.
  # x is the max of x
  x.max &lt;- ceiling(x/(10^(nchar(as.character(x))-1))) *
                  (10^(nchar(as.character(x))-1))
  return(x.max)
}

####

# Extract Expectation and variance from Weibull distribution:

Weibull.2p.Expectation &lt;- function(eta,beta){
  # In R, eta=scale and beta=shape parameters
  # of built-in functions.
  Expectn &lt;- eta * gamma(1 + 1/beta)
  return(Expectn)
}

Weibull.2p.Var &lt;- function(eta,beta){
  # In R, eta=scale and beta=shape parameters
  # of built-in functions.
  b &lt;- eta
  a &lt;- beta
  Var.W2p &lt;- b^2* (gamma(1 + 2/a) - (gamma(1 + 1/a))^2)
  return(Var.W2p)
}

##################

# For constructing Weibull plots:

F0inv    &lt;- function (p) log(qweibull(p, 1, 1))
# This is a fancy way to give you the y axis scale by setting
# shape and scale parameters to 1, p= median rank.

Weibull.2p.plot &lt;- function(x,y){
  # x = time; y = median ranks.
  ticks    &lt;- c(seq(0.01,0.09,0.01),(1:9)/10,seq(0.91,0.99,0.01))
  xticks &lt;- round(exp(seq(0.1,log(xmax(max(x))),
                    length.out=5)),0)
  y.trans &lt;- F0inv(y)
  plot(x,y.trans,xlim=c(exp(0.1),xmax(max(x))),
       ylim=F0inv(c(0.01,0.99)),log="x",axes=F)
  axis(1,at=xticks)
  axis(2,at=F0inv(ticks),labels=ticks)
  abline(h=F0inv(ticks),col="lightgray")
}

###########

# Produce failure rate plot: 2-parameter Weibull.

failure.rate.w2p &lt;- function(beta,eta,time){
  r &lt;- (beta/eta) *
        (time/eta)^(beta-1)
  return(r)
}

hazard.plot.w2p &lt;- function(beta,eta,time,line.colour,nincr=500){
  max.time &lt;- max(time,na.rm=F)
  t &lt;- seq(0,max.time,length.out=nincr)
  r &lt;- numeric(length(t))
  for(i in 1:length(t)){
    r[i] &lt;- failure.rate.w2p(beta,eta,t[i])
  }
  plot(t,r,type='l',bty='l',
       col=line.colour,lwd=2,
       main="",xlab="Time",
       ylab="Failure rate",
       las=1,adj=0.5,
       cex.axis=0.85,cex.lab=1.2)
}

##########

# Produce Reliability plot.

Reliability.w2p &lt;- function(beta,eta,time){
  R &lt;- exp(-(time/eta)^beta)
  return(R)
}

Reliability.plot.w2p &lt;- function(beta,eta,time,line.colour,nincr=500){
  max.time &lt;- max(time,na.rm=F)
  t &lt;- seq(0,max.time,length.out=nincr)
  R &lt;- numeric(length(t))
  for(i in 1:length(t)){
    R[i] &lt;- Reliability.w2p(beta,eta,t[i])
  }
  plot(t,R,type='l',bty='l',
       col=line.colour,lwd=2,
       main="",xlab="Time",
       ylab="Reliability",
       las=1,adj=0.5,
       cex.lab=1.2,cex.axis=0.85)
}

####

Plot.Observations &lt;- function(reliability.data,Ntotal=-999){
  # Specify Ntotal if plotting a subset of the data.
  dat &lt;- reliability.data

  fail &lt;- sort(dat[dat$event==1,"time"])
  # i.e., if there are censored observations:
  if(sum(dat$event)&lt;nrow(dat)){
    cens &lt;- sort(dat[dat$event==0,"time"])
    dat2 &lt;- data.frame(fail=fail,cens=NA)
    fail.df &lt;- data.frame(fail=NA,cens=cens)
    dat2 &lt;- rbind(dat2,fail.df)
    time.index &lt;- apply(dat2,1,sum,na.rm=T)
    dat3 &lt;- dat2[order(time.index),]
    par(yaxt="n")
    barplot(t(as.matrix(dat3[nrow(dat3):1,])),
            beside=T,horiz=T,col=c("black","red"),border=NA,
            xlab="time") # as per Fig 1.5 in Meeker &amp; Escobar.
    plot.dat &lt;- barplot(t(as.matrix(dat3[nrow(dat3):1,])),
                        beside=T,horiz=T,border=NA,plot=F)
    par(yaxt="s")
    N.dataset &lt;- ifelse(Ntotal==-999,nrow(dat),Ntotal)
    Start.y &lt;- N.dataset - nrow(dat) + 1
    Difference &lt;- abs(nrow(dat)-max(c(seq(1,nrow(dat3),by=5))))
    axis(2,at=plot.dat[1,c(seq(1,nrow(dat3),by=5))] + ceiling(plot.dat[2,Difference]),
         labels=rev(seq(Start.y,N.dataset,by=5)),
         las=1,adj=0.5,cex.axis=0.85)
  } else{
    barplot(t(as.matrix(fail[length(fail):1])),col="black",
            horiz=T,border=NA,
            xlab="time")
    plot.dat &lt;- barplot(t(as.matrix(fail[length(fail):1])),
                        horiz=T,border=NA,plot=F)
    Difference &lt;- abs(length(fail)-max(c(seq(1,length(fail),by=5))))
    axis(2,at=plot.dat[seq(1,length(fail),by=5)] + ceiling(plot.dat[2,Difference]),
         labels=rev(seq(1,length(fail),by=5)),
         las=1,adj=0.5,cex.axis=0.85)
  }
  legend("topright",
         legend=c("Failure","Suspension"),
         pch=15,
         col=c("black","red"),
         bty='n',
         horiz=F,xpd=NA,pt.cex=1.5)
  mtext("Time-ranked observation",side=2,line=2.5,adj=0.5)
}

###

Calculate.Fhat &lt;- function(reliability.data){
  # reliability.data has columns "time" and "event"={1,0}
  dat1 &lt;- reliability.data
  n &lt;- nrow(dat1)
  dat1$suspensions &lt;- 1 - dat1$event
  dat2 &lt;- aggregate(dat1[,2:3],list(time=dat1$time),sum)
  names(dat2)[2:3] &lt;- c("dj","rj")
  dat2$sum_dj &lt;- cumsum(dat2$dj)
  dat2$sum_rj &lt;- cumsum(dat2$rj)
  dat2$nj &lt;- 0 # set empty numeric vector
  m &lt;- nrow(dat2)
  attach(dat2)
  dat2$nj[1] &lt;- n
  for(i in 2:m){
    dat2$nj[i] &lt;- n - sum_dj[i-1] - sum_rj[i-1]
  }
  detach(dat2)
  dat2$pj &lt;- dat2$dj / dat2$nj
  dat2$qj &lt;- 1 - dat2$pj
  dat2$Shat &lt;- cumprod(dat2$qj)
  dat2$Fhat &lt;- 1 - dat2$Shat
  dat3 &lt;- dat2[dat2$dj&gt;0,colnames(dat2) %in% c("time","Fhat")==T]
  return(dat3)
}

Calculate.a_b &lt;- function(reliability.data){
  # reliability.data has columns "time" and "event"={1,0}
  dat1 &lt;- reliability.data
  n &lt;- nrow(dat1)
  dat1$suspensions &lt;- 1 - dat1$event
  dat2 &lt;- aggregate(dat1[,2:3],list(time=dat1$time),sum)
  names(dat2)[2:3] &lt;- c("dj","rj")
  dat2$sum_dj &lt;- cumsum(dat2$dj)
  dat2$sum_rj &lt;- cumsum(dat2$rj)
  dat2$nj &lt;- 0 # set empty numeric vector
  m &lt;- nrow(dat2)
  attach(dat2)
  dat2$nj[1] &lt;- n
  for(i in 2:m){
    dat2$nj[i] &lt;- n - sum_dj[i-1] - sum_rj[i-1]
  }
  detach(dat2)
  dat2$sigma_j &lt;- dat2$dj / (dat2$nj*(dat2$nj - dat2$dj))
  cumsum.sigma &lt;- cumsum(dat2$sigma_j)
  dat2$sigmahat &lt;- 0 # set empty numeric vector
  for(i in 2:m){
    dat2$sigmahat[i] &lt;- n * cumsum.sigma[i-1]
  }
  dat2$Khat &lt;- dat2$sigmahat / (1 + dat2$sigmahat)
  dat3 &lt;- dat2[,colnames(dat2) %in% c("time","sigmahat","Khat")==T]
  dat4 &lt;- dat3[dat3$Khat&gt; 0 &amp; dat3$Khat &lt;1, ] # Can't be zero or 1.
  a &lt;- min(dat4$Khat); b &lt;- max(dat4$Khat)
  result &lt;- list(a=a, b=b)
  return(result)
}

Calculate.e_val &lt;- function(ab.obj){
  a &lt;- ab.obj$a; b &lt;- ab.obj$b
# e_alpha is a global object created in this file (lookup table from Meeker &amp; Escobar)
  (e_val &lt;- e_alpha[e_alpha$a==e_alpha$a[which.min(abs(a - e_alpha$a))] &amp;
                      e_alpha$b==e_alpha$b[which.min(abs(b - e_alpha$b))],
                    "cl_0.95"])
  return(e_val)
}

Calc.95.simultaneous.CI &lt;- function(reliability.data,e_val){
  print(c("Adjustments to 95 % simultaneous confidence bounds to account for"))
  print(c("non-increasing values follow method of Meeker &amp; Escobar (1998)"))
  dat1 &lt;- reliability.data
  n &lt;- nrow(dat1)
  dat1$suspensions &lt;- 1 - dat1$event
  dat2 &lt;- aggregate(dat1[,2:3],list(time=dat1$time),sum)
  names(dat2)[2:3] &lt;- c("dj","rj")
  dat2$sum_dj &lt;- cumsum(dat2$dj)
  dat2$sum_rj &lt;- cumsum(dat2$rj)
  dat2$nj &lt;- 0 # set empty numeric vector
  m &lt;- nrow(dat2)
  attach(dat2)
  dat2$nj[1] &lt;- n
  for(i in 2:m){
    dat2$nj[i] &lt;- n - sum_dj[i-1] - sum_rj[i-1]
  }
  detach(dat2)
  dat2$pj &lt;- dat2$dj / dat2$nj
  dat2$qj &lt;- 1 - dat2$pj
  dat2$Shat &lt;- cumprod(dat2$qj)
  dat2$Fhat &lt;- 1 - dat2$Shat
  dat3 &lt;- dat2[dat2$dj&gt;0,]
  dat3$se.summation &lt;- dat3$pj / (dat3$nj*(1-dat3$pj))
  sum.term &lt;- cumsum(dat3$se.summation)
  attach(dat3)
  dat3$se &lt;- sqrt( (Shat)^2 * sum.term )
  detach(dat3)
  dat3$w &lt;- exp((e_val*dat3$se) / (dat3$Fhat*(1-dat3$Fhat)))
  attach(dat3)
  dat3$loUnadj &lt;- Fhat / (Fhat + (1-Fhat)*w)
  dat3$hiUnadj &lt;- Fhat / (Fhat + (1-Fhat)/w)
  detach(dat3)
  dat4 &lt;- dat3[is.nan(dat3$se)==F,]
  lo.NonDecreasing &lt;- ifelse(sum(diff(dat4$loUnadj)&lt;0)==0,"No","Yes")
  hi.NonDecreasing &lt;- ifelse(sum(diff(dat4$hiUnadj)&lt;0)==0,"No","Yes")
  if(lo.NonDecreasing=="Yes"){
    Max.lo &lt;- max(dat3$loUnadj,na.rm=T)
    dat3$lo=dat3$loUnadj
    dat3$lo[which.max(dat3$loUnadj):length(dat3$loUnadj)] &lt;- Max.lo
  }else{
    dat3$lo=dat3$loUnadj
  }
  if(hi.NonDecreasing=="Yes"){
    Min.hi &lt;- min(dat3$hiUnadj,na.rm=T)
    dat3$hi=dat3$hiUnadj
    dat3$hi[1:which.min(dat3$hiUnadj)] &lt;- Min.hi
  }else{
    dat3$hi=dat3$hiUnadj
  }
  return(dat3[dat3$Fhat&gt;0 &amp; dat3$Fhat &lt;1,
              colnames(dat3)%in%c("time","Fhat","lo","hi")])
}

####

# Reference table for e_{a,b,1-alpha/2} factors from Table 3.5
# of Meeker &amp; Escobar (1998)
e_alpha &lt;- data.frame(
  a = c(rep(c(0.005,0.01,0.05,0.1),4)),
  b = c(rep(c(0.995,0.99,0.95,0.9),each=4)),
  cl_0.95 = c(3.36,3.34,3.28,3.25,
              3.34,3.31,3.25,3.21,
              3.28,3.25,3.16,3.11,
              3.25,3.21,3.11,3.06)
)
####

# Functions to produce probability plots in Step 1 of analysis:

Normal.probability.plot &lt;- function(x,y,gridlines=F,
                                    label.individual.axes=T){
  # x = time; y = F(t).
  x &lt;- x[y &gt; 0 &amp; y &lt;1] # Can't be plotted on probability paper.
  y &lt;- y[y &gt; 0 &amp; y &lt;1]
  y.trans &lt;- qnorm(y)
  yticks &lt;- seq(round(min(y.trans),2),round(max(y.trans),2),
                length.out=5)
  xticks &lt;- round(seq(0,xmax(max(x)),length.out=5),0)
  if(label.individual.axes==T){X.label &lt;- "Time"
  } else{
    X.label &lt;- ""
  }
  if(label.individual.axes==T){Y.label &lt;- "Unreliability, F(t)=1-R(t)"
  } else{
    Y.label &lt;- ""
  }
  plot(x,y.trans,xlim=c(0,xmax(max(x))),
       ylim=c(round(min(y.trans),2),round(max(y.trans),2)),
       axes=F,pch=16,adj=0.5,
       main="Normal",xlab=X.label,ylab=Y.label)
  axis(1,at=xticks)
  axis(2,at=yticks,labels=sprintf("%.2f",round(pnorm(yticks),2)),
       las=1,adj=0.5)
  if(gridlines==T){
    abline(h=yticks,col="lightgray")
    abline(v=xticks,col="lightgray")
  }
}

Lognormal.probability.plot &lt;- function(x,y,gridlines=F,
                                       label.individual.axes=T){
  # x = time; y = F(t).
  x &lt;- x[y &gt; 0 &amp; y &lt;1] # Can't be plotted on probability paper.
  y &lt;- y[y &gt; 0 &amp; y &lt;1]
  x.trans &lt;- log(x)
  y.trans &lt;- qnorm(y)
  yticks &lt;- seq(round(min(y.trans),2),round(max(y.trans),2),
                length.out=5)
  xticks &lt;- round(seq(floor(min(x.trans)),ceiling(max(x.trans)),
                      length.out=5),0)
  if(label.individual.axes==T){X.label &lt;- "Time"
  } else{
    X.label &lt;- ""
  }
  if(label.individual.axes==T){Y.label &lt;- "Unreliability, F(t)=1-R(t)"
  } else{
    Y.label &lt;- ""
  }
  plot(x.trans,y.trans,xlim=c(floor(min(x.trans)),ceiling(max(x.trans))),
       ylim=c(round(min(y.trans),2),round(max(y.trans),2)),
       axes=F,pch=16,adj=0.5,
       main="Lognormal",xlab=X.label,ylab=Y.label)
  axis(1,at=xticks,labels=floor(exp(xticks)))
  axis(2,at=yticks,labels=sprintf("%.2f",round(pnorm(yticks),2)),
       las=1,adj=0.5)
  if(gridlines==T){
    abline(h=yticks,col="lightgray")
    abline(v=xticks,col="lightgray")
  }
}

add95CIs.Lognormal &lt;- function(CL.data){
  # CL.data is the data frame generated from using Calc.95.simultaneous.CI()
  time &lt;- CL.data$time
  lo &lt;- CL.data$lo # the lower 95% simultaneous confidence limit
  hi &lt;- CL.data$hi # the upper 95% simultaneous confidence limit
  points(log(time),qnorm(lo),pch="-",lwd=2,cex=1.2)
  points(log(time),qnorm(hi),pch="-",lwd=2,cex=1.2)
}

Weibull.backtrans.Y &lt;- function(y){1-(1/exp(exp(y)))}

Weibull.probability.plot &lt;- function(x,y,gridlines=F,
                                     label.individual.axes=T){
  # x = time; y = F(t).
  x &lt;- x[y &gt; 0 &amp; y &lt;1] # Can't be plotted on probability paper.
  y &lt;- y[y &gt; 0 &amp; y &lt;1]
  x.trans &lt;- log(x)
  y.trans &lt;- log(-log(1-y))
  yticks &lt;- seq(round(min(y.trans),2),round(max(y.trans),2),
                length.out=5)
  xticks &lt;- round(seq(floor(min(x.trans)),ceiling(max(x.trans)),
                      length.out=5),0)
  if(label.individual.axes==T){X.label &lt;- "Time"
  } else{
    X.label &lt;- ""
  }
  if(label.individual.axes==T){Y.label &lt;- "Unreliability, F(t)=1-R(t)"
  } else{
    Y.label &lt;- ""
  }
  plot(x.trans,y.trans,xlim=c(floor(min(x.trans)),ceiling(max(x.trans))),
       ylim=c(round(min(y.trans),2),round(max(y.trans),2)),
       axes=F,pch=16,adj=0.5,
       main="Weibull",xlab=X.label, ylab=Y.label)
  axis(1,at=xticks, labels=round(exp(xticks),0))
  axis(2,at=yticks,labels=sprintf("%.2f",round(Weibull.backtrans.Y(yticks),2)),
       las=1,adj=0.5)
  if(gridlines==T){
    abline(h=yticks,col="lightgray")
    abline(v=xticks,col="lightgray")
  }
}

add95CIs.Weibull &lt;- function(CL.data){
  # CL.data is the data frame generated from using Calc.95.simultaneous.CI()
  time &lt;- CL.data$time
  lo &lt;- CL.data$lo # the lower 95% simultaneous confidence limit
  hi &lt;- CL.data$hi # the upper 95% simultaneous confidence limit
  points(log(time),log(-log(1-lo)),pch="-",lwd=2,cex=1.2)
  points(log(time),log(-log(1-hi)),pch="-",lwd=2,cex=1.2)
}

Exponential.backtrans.Y &lt;- function(y){1-(1/exp(y))}

Exponential.probability.plot &lt;- function(x,y,gridlines=F,
                                         label.individual.axes=T){
  # x = time; y = F(t).
  x &lt;- x[y &gt; 0 &amp; y &lt;1] # Can't be plotted on probability paper.
  y &lt;- y[y &gt; 0 &amp; y &lt;1]
  y.trans &lt;- -log(1-y)
  yticks &lt;- seq(round(min(y.trans),2),
                round(max(y.trans),2),
                length.out=5)
  xticks &lt;- round(seq(0,xmax(max(x)),length.out=5),0)
  if(label.individual.axes==T){X.label &lt;- "Time"
  } else{
    X.label &lt;- ""
  }
  if(label.individual.axes==T){Y.label &lt;- "Unreliability, F(t)=1-R(t)"
  } else{
    Y.label &lt;- ""
  }
  plot(x,y.trans,xlim=c(0,xmax(max(x))),
       ylim=c(round(min(y.trans),2),
              round(max(y.trans),2)),
       axes=F,pch=16,adj=0.5,
       main="Exponential",xlab = X.label, ylab = Y.label)
  axis(1,at=xticks)
  axis(2,at=yticks,labels=sprintf("%.2f",round(Exponential.backtrans.Y(yticks),2)),
       las=1,adj=0.5)
  if(gridlines==T){
    abline(h=yticks,col="lightgray")
    abline(v=xticks,col="lightgray")
  }
}

add95CIs.Exponential &lt;- function(CL.data){
  # CL.data is the data frame generated from using Calc.95.simultaneous.CI()
  time &lt;- CL.data$time
  lo &lt;- CL.data$lo # the lower 95% simultaneous confidence limit
  hi &lt;- CL.data$hi # the upper 95% simultaneous confidence limit
  points(time,-log(1-lo),pch="-",lwd=2,cex=1.2)
  points(time,-log(1-hi),pch="-",lwd=2,cex=1.2)
}

Probability.Plots &lt;- function(reliability.data,gridlines=F,
                              label.individual.axes=F,dist="All"){
  dat &lt;- reliability.data
  Fhat &lt;- Calculate.Fhat(dat) # $time, $Fhat
  e_val &lt;- Calculate.e_val(Calculate.a_b(dat))
  simult.CLs &lt;- Calc.95.simultaneous.CI(dat,e_val=e_val) # $time, $lo, $hi
  if(dist=="All"){
  par(mfrow=c(2,2), mar=c(3, 3, 1, 0.25) + 0.1,cex.axis=0.75,
      oma=c(0,0,0,0),mgp=c(3,1,0),xpd=T,mai=c(0.75,0.75,0.2,0.2))
  # Plot for Weibull distbn.
  Weibull.probability.plot(Fhat$time,Fhat$Fhat,gridlines,
                           label.individual.axes)
  add95CIs.Weibull(simult.CLs)
  # Plot for Lognormal distbn.
  Lognormal.probability.plot(Fhat$time,Fhat$Fhat,gridlines,
                             label.individual.axes)
  add95CIs.Lognormal(simult.CLs)
  # Plot for Normal distbn.
  Normal.probability.plot(Fhat$time,Fhat$Fhat,gridlines,
                          label.individual.axes)
  points(simult.CLs$time,qnorm(simult.CLs$lo),pch="-",lwd=2,cex=1.2)
  points(simult.CLs$time,qnorm(simult.CLs$hi),pch="-",lwd=2,cex=1.2)
  # Plot for Exponential cdf.
  Exponential.probability.plot(Fhat$time,Fhat$Fhat,gridlines,
                               label.individual.axes)
  add95CIs.Exponential(simult.CLs)
  mtext("Time",side=1,line=-1.5,adj=0.5,outer=T,font=2)
  mtext("Unreliability, F(t)=1-R(t)",side=2,line=-1.5,adj=0.5,
        outer=T,font=2)
  }
  else if(dist=="Weibull"){
    par(mar= c(5, 5, 4, 1) + 0.1,font.lab=2,cex.axis=0.8,cex.lab=1.1,
        cex.main=1.3)
    Weibull.probability.plot(Fhat$time,Fhat$Fhat,gridlines,T)
    add95CIs.Weibull(simult.CLs)
  }
  else if(dist=="Normal"){
    par(mar= c(5, 5, 4, 1) + 0.1,font.lab=2,cex.axis=0.8,cex.lab=1.1,
        cex.main=1.3)
    Normal.probability.plot(Fhat$time,Fhat$Fhat,gridlines,T)
    points(simult.CLs$time,qnorm(simult.CLs$lo),pch="-",lwd=2,cex=1.2)
    points(simult.CLs$time,qnorm(simult.CLs$hi),pch="-",lwd=2,cex=1.2)
  }
  else if(dist=="Lognormal"){
    par(mar= c(5, 5, 4, 1) + 0.1,font.lab=2,cex.axis=0.8,cex.lab=1.1,
        cex.main=1.3)
    Lognormal.probability.plot(Fhat$time,Fhat$Fhat,gridlines,T)
    add95CIs.Lognormal(simult.CLs)
  }
  else{
    par(mar= c(5, 5, 4, 1) + 0.1,font.lab=2,cex.axis=0.8,cex.lab=1.1,
        cex.main=1.3)
    Exponential.probability.plot(Fhat$time,Fhat$Fhat,gridlines,T)
    add95CIs.Exponential(simult.CLs)
  }
  par(mfrow=c(1,1),mar= c(5, 4, 4, 2) + 0.1,cex.axis=1,xpd=F,
      oma=c(0,0,0,0),lheight=1,mgp=c(3,1,0),mai=c(1.02,0.82,0.82,0.42)) # return defaults.
}

###

# Bias-adjusted non-parametric bootstrapping to estimate approx 95% CIs for MTTF.

# "Bias-corrected percentile" method e.g., Section 13.7, Meeker &amp; Escobar.

# NB. Jeng and Meeker (1998 in Meeker &amp; Escobar 1998) demonstrated that
# simulation-based methods, for most situations, provide important
# improvements over normal-approximation methods.

MTTF.boot.percentile.adj &lt;- function(data,i){
  d &lt;- data[i,]
  mod &lt;- Lifedata.MLE(Surv(time,event)~1,
                      d,
                      dist="weibull")
  beta.b &lt;- 1/unname(exp(mod$coef[2]))
  eta.b &lt;- unname(exp(mod$coef[1]))
  MTTF &lt;- Weibull.2p.Expectation(eta=eta.b,
                                 beta=beta.b)
  return(MTTF)
}

###

# For calculating joint confidence region:

sev.pdf &lt;- function(z){exp(z - exp(z))}
sev.cdf &lt;- function(z){1-exp(-exp(z))}

loglik.sev &lt;- function(data,mu,sigma){
  t &lt;- data$time
  ll.vec &lt;- numeric(nrow(data))
  delta &lt;- data$event
  t.lnorm &lt;- (log(t)-mu)/sigma
  for(i in 1:nrow(data)){
    ll.vec[i] &lt;- (delta[i] * log(1 / (sigma*t[i]))) +
                 (delta[i] *
                    log(sev.pdf(t.lnorm[i]))) +
                 ((1-delta[i]) *
                    log(1 - sev.cdf(t.lnorm[i])))
  }
  loglik &lt;- sum(ll.vec)
  return(loglik)
}

contour.val &lt;- function(dataset,mu.input,sigma.input,maximum.loglik){
  pchisq(q=-2*(loglik.sev(dataset,mu.input,sigma.input) -
                 maximum.loglik),df = 2)
}

Get.contour.plot.data &lt;- function(data,fitted.2parameter.Weibull.model,steps){
  mod &lt;- fitted.2parameter.Weibull.model # 2 parameter Weibull model fitted using SPREDA
  mu.MLE &lt;- unname(mod$coef[1])
  sigma.MLE &lt;- unname(exp(mod$coef[2]))
  Max.loglik &lt;- loglik.sev(data,mu=mu.MLE,sigma=sigma.MLE)
  # Use 95% CIs for MLEs to determine input range for parameters.
  beta.95cl_hi &lt;- 1 / (summary(mod)$coefmat["sigma","95% Lower"])
  beta.95cl_lo &lt;- 1 / (summary(mod)$coefmat["sigma","95% Upper"])
  eta.95cl_lo &lt;- exp(summary(mod)$coefmat["(Intercept)","95% Lower"])
  eta.95cl_hi &lt;- exp(summary(mod)$coefmat["(Intercept)","95% Upper"])
  # Input ranges (steps)
  Betas &lt;- seq(0.9*beta.95cl_lo,1.1*beta.95cl_hi,length.out=steps)
  Etas &lt;- seq(0.9*eta.95cl_lo,1.1*eta.95cl_hi,length.out=steps)
  Sigmas &lt;- 1 / Betas
  Mus &lt;- log(Etas)
  ContourVals.df &lt;- data.frame(Mus=rep(Mus,steps),
                               Sigmas=rep(Sigmas,each=steps))
  for(i in 1:nrow(ContourVals.df)){
    ContourVals.df$z[i] &lt;- contour.val(data,ContourVals.df$Mus[i],
                                       ContourVals.df$Sigmas[i],Max.loglik)
  }
  ContourVals.df$Eta &lt;- exp(ContourVals.df$Mus)
  ContourVals.df$Beta &lt;- 1 / ContourVals.df$Sigmas
  return(ContourVals.df)
  # z vector is the probability that Chi-square stat with 2d.f. is less than
  # or equal to -2log likelihood ratio, for each input pair of Betas,Etas. Uses
  # the large sample Chi-square approximation for the distribution of the log-
  # likelihood ratio statistic.
}

Weibull.Confidence.Region &lt;- function(data,model,probability,
                                      title,show.contour.labels,
                                      steps=100,html="No"){
  # requires lattice package.
  # Implemented for the fit of the 2-parameter Weibull model only.
  beta.MLE &lt;- 1 / unname(exp(model$coef[2]))
  eta.MLE &lt;- unname(exp(model$coef[1]))
  reliability.data &lt;- data
  fitted.Weibull.model &lt;- model
  label.show &lt;- as.logical(show.contour.labels)
  Contour.df &lt;- Get.contour.plot.data(reliability.data,
                                      fitted.Weibull.model,steps)
  if(html=="Yes"){ # this right align is not consistent for html generation.
    title.settings &lt;- list(
      par.main.text = list(font = 2,
                          just = "right",
                          x = grid::unit(170, "mm")))
  } else{ # Use default settings.
    title.settings &lt;- list(par.main.text=trellis.par.get("par.main.text"))
  }
  contourplot(z ~ Eta * Beta, data=Contour.df,
                      at = probability, labels=label.show,
                      panel=function(data, ...){
                        panel.contourplot(...)
                        panel.points(x=eta.MLE,y=beta.MLE,pch=19,col='black')
                      },
                      par.settings=title.settings,
                      main=title)
}

###########

# Some functions used for internal referencing in R markdown:
# chunkref &lt;- local({
  # function(chunklabel) {
    # sprintf('[%s](#%s)', chunklabel, chunklabel )
  # }
# })

# secref &lt;- local({
  # function(seclabel) {
    # sprintf('[%s](#%s)', seclabel, seclabel )
  # }
# })

# pgref &lt;- local({
  # function(n)
    # sprintf('[Page-%i](#Page-%i)', n, n)
# })

# sec &lt;- local({
  # function(seclabel) {
    # sprintf('# &lt;a name="%s"/&gt; %s', seclabel, seclabel )
  # }
# })

# pgcount &lt;- local({
  # pg &lt;- 0
  # function(inc=T) {
    # if( inc ) { pg &lt;&lt;- pg + 1 }
    # return( pg )
  # }
# })

# pganchor &lt;- local({
  # function(doLabel=T) {
    # if( doLabel) {
      # sprintf('\n-----\nPage-%i\n&lt;a name="Page-%i"/&gt;\n', pgcount(inc=F), pgcount() )
    # } else {
      # sprintf('\n&lt;a name="Page-%i"/&gt;\n', pgcount() )
    # }
  # }
# })

# knit_hooks$set( anchor = function(before, options, envir) {
  # if ( before ) {
    # sprintf('&lt;a name="%s"/&gt;\n', options$label )
  # }
# })

# knit_hooks$set( echo.label = function(before, options, envir) {
  # if ( before ) {
    # sprintf('&gt; %s', options$label )
  # }
# })

# knit_hooks$set( pgbreak = function(before, options, envir) {
  # if ( !before ) {
    # pganchor();
  # }
# })
</pre></body></html>Ztext/plain