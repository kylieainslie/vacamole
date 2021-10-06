##################################################################################
# A script with helpful R functions for statistical data analysis
#
# Author: S Towers
#
# Created: Sept 5, 2013
#
# This script is not guaranteed to be free of bugs and/or errors.
#
# This script can be freely used and shared as long as the author and
# copyright information in this header remain intact.
##################################################################################


##################################################################################
# R has built-in random number generators for the binomial distibution,
# but they are expressed in different form than the form in the Lloyd-Smith
# PLoS paper.  This function translates the mu and alpha parameters
# into the format needed by R to generate the NB random numbers
##################################################################################
my_rnbinom=function(n,mu,alpha){
  size = 1/alpha
  prob = size/(mu+size)
  return(rnbinom(n,size=size,prob=prob))
}


########################################################################
########################################################################
# some shorthand functions so we don't have to type !(x%in%y)
# but instead can type x%!in%y
########################################################################
'%notin%' = function(x,y)!('%in%'(x,y))
'%!in%' = function(x,y)!('%in%'(x,y))


##################################################################################
##################################################################################
# https://stackoverflow.com/questions/14423325/confidence-intervals-for-predictions-from-logistic-regression
##################################################################################
get_prediction_and_confidence_intervals_from_binomial_fit=function(mod,newdata,critval=1.96){
  preds = predict(mod, newdata , type = "link", se.fit = TRUE)  
  upr = preds$fit + (critval * preds$se.fit)
  lwr = preds$fit - (critval * preds$se.fit)
  fit = preds$fit
  fit2 = mod$family$linkinv(fit)
  upr2 = mod$family$linkinv(upr)
  lwr2 = mod$family$linkinv(lwr)
  y = cbind(fit2,lwr2,upr2)
  return(y)
}

##################################################################################
##################################################################################
# this function histograms count data, and overlays the expected distribution
# given some glm(family="poisson") model fit results
##################################################################################
overlay_expected_distribution_from_poisson_glm_fit = function(count_data,glm_model_object){
  
  ################################################################################
  # for each data point, given the model prediction for that data point
  # and Poisson stochasticity, calculate the expected distribution
  ################################################################################
  xmin = min(count_data)
  xmax = max(count_data)
  vx = seq((xmin-5),(xmax+5))
  vexpected = rep(0,length(vx))
  for (i in 1:length(count_data)){
    vexpected = vexpected + dpois(vx,glm_model_object$fit[i])
  }
  
  ################################################################################
  # now histogram the data
  ################################################################################
  breaks = seq((min(vx)-0.5),(max(vx)+0.5))
  a = hist(count_data,breaks=breaks,plot=F)
  ymax = max(c(a$counts,vexpected))
  plot(a,ylim=c(0,1.15*ymax),col=2,xlab="Distribution of data",border=2,main="Histogram of data")
  
  ################################################################################
  # overlay the expected distribution
  ################################################################################
  lines(vx,vexpected,col=3,lwd=3)
  legend("topleft",legend=c("Distribution of data","Expected from model and Poisson stochasticity"),col=c(2,3),lwd=4,bty="n",cex=0.7)
  
  ################################################################################
  #https://stat.ethz.ch/pipermail/r-help/2017-April/446330.html
  ################################################################################
  vpois = numeric(0)
  for (iter in 1:100){
    vpois = c(vpois,rpois(length(count_data),glm_model_object$fit))
  }
  qqplot(vpois
         ,count_data
         ,xlab="Theoretical quantiles from best-fit model"
         ,ylab="Observed quantiles"
         ,main="Poisson QQplot")
  lines(c(0,1e6),c(0,1e6),col=3,lwd=3)
  legend("topleft",legend=c(paste("Expected when data are Poisson\n","distributed about model expectation",sep="")),col=3,lwd=3,bty="n",cex=0.7)
  
}

##################################################################################
##################################################################################
# R function to convert month, day and year to information like the
# date expressed as a fraction of the year (date), the weekday, the number of
# days since 1/1/1970 (jul), and the day_of_year
##################################################################################
convert_month_day_year_to_date_information = function(month,day,year){
  vjul = julian(month,day,year)
  date=as.Date(vjul,origin="1970-01-01")
  vweekday=as.numeric(format(date,"%w"))
  vday_of_year=as.numeric(format(date,"%j"))
  vyear=as.numeric(format(date,"%Y"))
  vmonth=as.numeric(format(date,"%m"))
  vday=as.numeric(format(date,"%d"))
  vdate = vyear+(vday_of_year-0.5)/365
  lgood = which(!is.na(vyear)&vyear%%4==0)
  vdate[lgood] = vyear[lgood]+(vday_of_year[lgood]-0.5)/366
  return(list(date=vdate,weekday=vweekday,jul=vjul,day_of_year=vday_of_year))
}


##################################################################################
##################################################################################
# function to return the number of daylight hours by day_of_year at 
# some degree of latitude
##################################################################################
length_of_day = function(day_of_year,lat){
  lat=lat*pi/180
  P = asin(0.39795*cos(0.2163108+2*atan(0.9671396*tan(0.00860*(day_of_year-186)))))
  D = 24-(24/pi)*acos((sin(0.8333*pi/180)+sin(lat)*sin(P))/(cos(lat)*cos(P)))
  return(D)
}

##################################################################################
##################################################################################
##################################################################################
create_continuous_time_series = function(t_expect,t,y){
  ynew = rep(NA,length(t_expect))
  ynew[t_expect%in%t] = y
  return(list(t=t_expect,y=ynew))
}

##################################################################################
# a function to overlay a Normal distribution with the same mean and variance as
# the data over the data.  
# For example, x might be a linear model fit residuals
##################################################################################
norm_overlay=function(x,nbins=20,xlo=min(x),xhi=max(x),main="",weights=rep(1,length(x)),breaks=NULL){
  xmin = min(c(x*weights,mean(x)-4*sd(x)))
  xmax = max(c(x*weights,mean(x)+4*sd(x)))
  xlim = c(xmin,xmax)
  if (length(breaks)==0){
    a = hist(x*weights,breaks=seq(xmin,xmax,length=nbins),xlab="Histogram of data",ylab="\043 per bin",main=main,xlim=xlim,col=3)
  }else{
    xlim = c(min(breaks),max(breaks))
    a = hist(x*weights,breaks=breaks,xlab="Histogram of data",ylab="\043 per bin",main=main,xlim=xlim,col=3)
  }
  
  # overlay a Normal distribution with mean and std deviation the same as that
  # of the data
  amids = a$mids
  apred = dnorm(amids,mean=mean(x),sd=sd(x))
  apred = length(x)*apred/(sum(apred))
  lines(amids,apred,col=2,lwd=4)
  legend("topleft",legend=c("Data","Normal distribution"),col=c(3,2),lwd=5,bty="n",cex=0.7)
}

##################################################################################
# a function to overlay a Normal distribution with the same mean and variance as
# the data over the data.  Also does the QQ-plot
# For example, x might be a linear model fit residuals
##################################################################################
norm_overlay_and_qq=function(x,nbins=20,xlo=min(x),xhi=max(x),weights=rep(1,length(x)),main="",breaks=NULL){
  
  xmin = min(c(x,mean(x)-4*sd(x)))
  xmax = max(c(x,mean(x)+4*sd(x)))
  xlim = c(xmin,xmax)
  
  if (length(breaks)==0){
    a = hist(x*weights,breaks=seq(xmin,xmax,length=nbins),xlab="Histogram of data",ylab="\043 per bin",main=main,xlim=xlim,col=3)
  }else{
    xlim = c(min(breaks),max(breaks))
    a = hist(x*weights,breaks=breaks,xlab="Histogram of data",ylab="\043 per bin",main=main,xlim=xlim,col=3)
  }
  
  # overlay a Normal distribution with mean and std deviation the same as that
  # of the data
  amids = a$mids
  apred = dnorm(amids,mean=mean(x),sd=sd(x))
  apred = length(x)*apred/(sum(apred))
  lines(amids,apred,col=2,lwd=4)
  legend("topleft",legend=c("Data","Normal distribution"),col=c(3,2),lwd=5,bty="n",cex=0.7)
  
  qqnorm((x-mean(x))*weights/sd(x),cex=1.5,ylab="Normalized data values",xlab="Expected value for Normal data",pch=20)
  qqline((x-mean(x))*weights/sd(x),col=3,lwd=2)
  
}


##################################################################################
# A function to calculate the standard error on the mean of the sample
##################################################################################
se = function(x){
  return(sd(x)/sqrt(length(x)))
}

##################################################################################
# A function to summarize the analysis of variance of a model, to get the
# percentage of the total sums of squares described by the 
# The total sum of squares a model is the explained sum of squares plus the
# residual sum of squares
# The explained sum of squares is Sum((hat(y)-mean(y))^2)
##################################################################################
anova_summary = function(myfit){
  vname = attr(summary(myfit)$terms,"term.labels")
  y = myfit$fitted.values+myfit$residuals
  totsumsq = sum((y-mean(y))^2)
  myanova = anova(myfit)
  vss = myanova[1:length(vname),2]
  pertot = vss/totsumsq
  cat("Fraction of total sum of squares described by each variable added in turn:\n")
  for (i in 1:length(vname)){
    cat("  ",vname[i],pertot[i],"\n")
  }
}

##################################################################################
##################################################################################
factorplot = function(y,ey,add=F,col=1,colb="grey"){
  x = seq(1,length(y))
  scale = 1.5
  ymin = min(y-scale*ey)
  ymax = max(y+scale*ey)
  xmin = min(x)-0.5
  xmax = max(x)+0.5
  sig = sd(y)
  mu = mean(y)
  ymin = min(ymin,(mu-scale*sig))
  ymax = max(ymax,(mu+scale*sig))
  if (add){
    points(x,y,ylim=c(ymin,ymax),pch=20,col=col,xlim=c(xmin,xmax),xaxs="i",yaxs="i",xlab="factor levels",ylab="average within level")
  }else{
    plot(x,y,ylim=c(ymin,ymax),pch=20,col=col,xlim=c(xmin,xmax),xaxs="i",yaxs="i",xlab="factor levels",ylab="average within level")
  }
  for (i in 1:length(x)){
    lines(c(x[i],x[i]),c(y[i]-ey[i],y[i]+ey[i]),col=col,lwd=4)
  }
  if (!add){
    legend("bottomleft",legend=c("mean and sd within factor level","mean across factor levels","+/- 1sd for entire sample"),lty=c(1,1,3),col=c(col,colb,colb),lwd=4)
    lines(c(xmin,xmax),c(mean(y),mean(y)),col=colb,lwd=3)
    lines(c(xmin,xmax),c(mean(y)-sig,mean(y)-sig),col=colb,lwd=3,lty=3)
    lines(c(xmin,xmax),c(mean(y)+sig,mean(y)+sig),col=colb,lwd=3,lty=3)
  }
}


##################################################################################
##################################################################################
##################################################################################
harmonic_regression=function(x,y,TestPeriod,plot=T){
  xa = sin(2*pi*x/TestPeriod)
  xb = cos(2*pi*x/TestPeriod)
  b = lm(y~xa+xb)
  V = vcov(b)
  beta1 = b$coefficients[2]
  beta2 = b$coefficients[3]
  
  beta1 = as.numeric(beta1)
  beta2 = as.numeric(beta2)
  
  A = sqrt(beta1^2+beta2^2)
  fpartial = c(0,beta1/A,beta2/A)
  eA = sqrt(t(fpartial)%*%V%*%(fpartial))
  
  phi = atan2(beta1,beta2)
  fpartial = c(0,-beta1/A^2,+beta2/A^2)
  ephi = sqrt(t(fpartial)%*%V%*%(fpartial))
  ephi = ephi*sqrt(2*pi)
  
  phi = phi*TestPeriod/(2*pi)
  ephi = ephi*TestPeriod/(2*pi)
  
  # pamplitude is the probability that the null hypothesis A=0 is true
  # pphase is the probability that the null hypothesis phi=0 is true
  pamplitude = 1-pchisq((A/eA)^2,2)
  pphase = 1-pchisq((phi/ephi)^2,1)
  if (plot){
    #par(mfrow=c(1,1))
    plot(x,y,lwd=3,type="l",main=paste("Testing period =",TestPeriod))
    lines(x,b$fitted.values,col=2,lwd=3)
  }
  return(list(amplitude=A,eamplitude=eA,phase=phi,ephase=ephi,pamplitude=pamplitude,pphase=pphase,fit=b))
}


##################################################################################
##################################################################################
##################################################################################
harmonic_regression_multi=function(x,y,TestPeriods,plot=T){
  y = y-mean(y)
  d = data.frame(y=y)
  for (i in 1:length(TestPeriods)){
    aname = paste("xa",i,sep="")
    bname = paste("xb",i,sep="")
    a=(sin(2*pi*x/TestPeriods[i]))
    b=(cos(2*pi*x/TestPeriods[i]))
    g = data.frame(a,b)
    d = cbind(d,g)
    names(d)[2+(i-1)*2] = aname
    names(d)[3+(i-1)*2] = bname
  }
  formula="y~.-1"
  b = lm(formula=formula,data=d)
  V = vcov(b)
  wa = numeric(0)
  ewa = numeric(0)
  wphi = numeric(0)
  ewphi = numeric(0)
  wpamplitude = numeric(0)
  wpphase = numeric(0)
  mjacob = matrix(0,ncol=(2*length(TestPeriods)),nrow=(length(TestPeriods)*2))
  for (i in 1:length(TestPeriods)){
    beta1 = b$coefficients[1+(i-1)*2]
    beta2 = b$coefficients[2+(i-1)*2]
    A = sqrt(beta1^2+beta2^2)
    fpartial = rep(0,length(b$coefficients))
    fpartial[1+(i-1)*2] = beta1/A
    fpartial[2+(i-1)*2] = beta2/A
    mjacob[1+(i-1)*2,i] = beta1/A
    mjacob[2+(i-1)*2,i] = beta2/A
    eA = sqrt(t(fpartial)%*%V%*%(fpartial))
    wa = append(wa,A)
    ewa = append(ewa,eA)
    
    phi = atan2(beta1,beta2)
    fpartial = rep(0,length(b$coefficients))
    fpartial[1+(i-1)*2] = -beta1/A^2
    fpartial[2+(i-1)*2] = +beta2/A^2
    mjacob[1+(i-1)*2,i+length(TestPeriods)] = -beta1/A^2
    mjacob[2+(i-1)*2,i+length(TestPeriods)] = +beta2/A^2
    ephi = sqrt(t(fpartial)%*%V%*%(fpartial))
    ephi = ephi*sqrt(2*pi)
    phi = phi*TestPeriods[i]/(2*pi)
    ephi = ephi*TestPeriods[i]/(2*pi)
    wphi = append(wphi,phi)
    ewphi = append(ewphi,ephi)
    
    # pamplitude is the probability that the null hypothesis A=0 is true
    # pphase is the probability that the null hypothesis phi=0 is true
    pamplitude = 1-pchisq((A/eA)^2,2)
    pphase = 1-pchisq((phi/ephi)^2,1)
    wpamplitude = append(wpamplitude,pamplitude)
    wpphase = append(wpphase,pphase)
  }
  ###########################################
  # calculate the covariance matrix
  # of the parameter estimates.  The
  # amplitudes are in the upper block
  # and the phases are in the lower block
  ###########################################
  D = t(mjacob)%*%V%*%mjacob
  Dcor = cov2cor(D)
  vsigma = sqrt(diag(D))
  vsigma[(length(TestPeriods)+1):(2*length(TestPeriods))] =  vsigma[(length(TestPeriods)+1):(2*length(TestPeriods))]*TestPeriods/sqrt(2*pi)
  Dnew = D
  for (i in 1:(2*length(TestPeriods))){
    for (j in 1:(2*length(TestPeriods))){
      Dnew[i,j] = Dcor[i,j]*vsigma[i]*vsigma[j]
    }
  }
  ###########################################
  # plot the results
  ###########################################
  if (plot){
    #par(mfrow=c(1,1))
    plot(x,y,lwd=3,type="l",main=paste("Testing periods =",TestPeriods))
    lines(x,b$fitted.values,col=2,lwd=3)
  }
  return(list(amplitude=wa,eamplitude=ewa,phase=wphi,ephase=ewphi,pamplitude=wpamplitude,pphase=wpphase,vcov=Dnew,vcor=Dcor,fit=b))
}

##################################################################################
##################################################################################
##################################################################################
lomb_freq=function(t,y,TestFrequencies,main="",xlab="Frequency",ylab="1-(pvalue)",plot=T,alpha=0.05,plotp=T,add=F,col=4){
  if (min(abs(diff(TestFrequencies)))<1/(length(y)+1)){
    cat("lomb_freq: the minimum difference in frequencies is ",min(abs(diff(TestFrequencies))),"\n")
    cat("lomb_freq: the difference in frequencies has to be at least ",1/length(y),"in order to obtain independent tests of significance\n")
  }
  a = lomb(t,y,1/TestFrequencies,plot=F)
  if (plot){
    if (!add){
      plot(TestFrequencies,1-a$pamplitude,
           pch=20, col=col,
           ylab=ylab, xlab=xlab,
           main=main,
           ylim=c(0,1.000),xaxs="i",yaxs="i",cex=0.2)
      lines(1/TestFrequencies,1-a$pamplitude,col=col,lwd=3)
      if (plotp) lines(c(-1000,10000),c(a$pbonferroni,a$pbonferroni),col=2,lwd=2,lty=3)
    }else{
      lines(TestPeriods,1-a$pamplitude,col=col,lwd=3)
    }
  }
  
  return(list(freq=TestFrequencies,pamplitude=a$pamplitude,vsignificant=(a$vsignificant),pbonferroni=a$pbonferroni))
}

##################################################################################
##################################################################################
##################################################################################
lomb=function(t,y,TestPeriods,main="",xlab="Period",ylab="1-(pvalue)",plot=T,alpha=0.05,plotp=T,add=F,col=4){
  TestFrequencies = 1/TestPeriods
  if (min(abs(diff(TestFrequencies)))<1/(length(y)+1)){
    cat("lomb: the minimum difference in frequencies is ",min(abs(diff(TestFrequencies))),"\n")
    cat("lomb: the difference in frequencies has to be at least ",1/length(y),"in order to obtain independent tests of significance\n")
  }
  N = as.integer(length(TestPeriods)/2)
  LS=LombScargle(t,y,TestFrequencies,N)
  pbonferroni = (1-alpha/length(TestPeriods))
  if (plot){
    if (!add){
      plot(1/TestFrequencies,1-LS$Probability,
           pch=20, col=col,
           ylab=ylab, xlab=xlab,
           main=main,
           ylim=c(0,1.000),xaxs="i",yaxs="i",cex=0.2)
      lines(1/TestFrequencies,1-LS$Probability,col=col,lwd=3)
      if (plotp) lines(c(-1000,10000),c(pbonferroni,pbonferroni),col=2,lwd=2,lty=3)
    }else{
      lines(TestPeriods,1-LS$Probability,col=col,lwd=3)
    }
  }
  vsignificant =((1-LS$Probability)>pbonferroni)
  
  # pamplitude is the probability that the null hypothesis amplitude=0 is true
  return(list(T=1/TestFrequencies,pamplitude=(LS$Probability),vsignificant=(vsignificant),pbonferroni=pbonferroni))
}

##################################################################################
##################################################################################
# Lomb Scargle periodogram
# from
# http://research.stowers-institute.org/efg/2005/LombScargle/R/cosine.R
#
# example of its use in ~/epidemiology/fits/data/geneva/geneva.R
##################################################################################
LombScargle <- function(t, h, TestFrequencies, Nindependent)
{
  stopifnot( length(t) == length(h) )
  
  hResidual    <- h - mean(h)
  
  SpectralPowerDensity <- rep(0, length(TestFrequencies))
  
  for (i in 1:length(TestFrequencies))
  {
    Omega       <- 2*pi*TestFrequencies[i]
    TwoOmegaT   <- 2*Omega*t
    Tau         <- atan2( sum(sin(TwoOmegaT)) , sum(cos(TwoOmegaT))) /
      (2*Omega)
    OmegaTMinusTau <- Omega * (t - Tau)
    SpectralPowerDensity[i] <- (sum (hResidual * cos(OmegaTMinusTau))^2) /
      sum( cos(OmegaTMinusTau)^2 ) +
      (sum (hResidual * sin(OmegaTMinusTau))^2) /
      sum( sin(OmegaTMinusTau)^2 )
  }
  
  # The "normalized" spectral density refers to the variance term in the
  # denominator.  With this term the SpectralPowerDensity has an
  # exponential probability distribution with unit mean.
  SpectralPowerDensity <- SpectralPowerDensity / ( 2 * var(h) )
  
  Probability <- 1-((1-exp(-SpectralPowerDensity)))
  
  PeakIndex    <- match(max(SpectralPowerDensity), SpectralPowerDensity)
  LikelyPeriod <- 1 / TestFrequencies[PeakIndex]
  LikelyPvalue <- Probability[PeakIndex]
  
  return( list( Probability=Probability,
                SpectralPowerDensity=SpectralPowerDensity,
                PeakIndex=PeakIndex,
                LikelyPeriod=LikelyPeriod,
                LikelyPvalue=LikelyPvalue
  ) )
}




#qqplot(qgamma(ppoints(length(x)),shape=10,scale=1),y,xlab="")