#===========================================================================
#   1.0 Set Modeling Environment 
#===========================================================================
# Clear data 
rm(list=ls(all=TRUE))
# load packages
library(coda)
library(R2jags)

#---------------------------------------------------------------------------
#   1.1 Set Data input output file and directory 
#---------------------------------------------------------------------------
# Set directory where data file is located
data_dir <- 'c:/Projects/BEG/UCI/'
# set data name 

#file_name <- 'Kasilof_Sockeye.csv'
#results_file <- 'Kasilof_JAGS.txt'
file_name <- 'Kenai_Sockeye.csv'
results_file <- 'Kenai_JAGS.txt'
results_file2 <- 'Kenai_JAGS79.txt'

dat1 <- read.csv(paste(data_dir,file_name, sep=''), na.strings ='', header = TRUE)
# Run this for years from 1979 - 2012
#dat1 <- dat1[dat1$Year >=1979,]

t_col <- function(color, percent = 50, name = NULL) {
#	  color = color name
#	percent = % transparency
#	   name = an optional name for the color
## Get RGB values for named color
  rgb.val <- col2rgb(color)

## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100,
               names = name)
## Save the color
  invisible(t.col)
  }
  
# Number of years 
nyrs <- length(dat1$Year)
# Spawner
S <- (dat1$Spawner)
# Recruit
R <- (dat1$Return)
# Mulitplier (This will make S small, and improve beta estimate) 
d <- floor(log10(mean(S)))

#---------------------------------------------------------------------------
#  1.2: Create JAGS data file (datnew)				               
#---------------------------------------------------------------------------
datnew<-list(nyrs=nyrs,S=S,R=R,d = d,ar1=1,rw=0)

#===========================================================================	
#  2.0: JAGS SR Model code					   
#===========================================================================
parameters.CR <- c('lnalpha','beta','phi','e0','sigma','lnalphai') 
  jag.model.CR <- function(){
    for(y in 1:nyrs){
      s[y] <- S[y]/(10^d)
      fit[y] <- log(S[y]) + lnalpha - beta * s[y] + rw*v[y]
      e[y] <- log(R[y]) - fit[y]
	  w[y] ~dnorm(0,tauw)
	  v[y] ~dnorm(0,tauv)
     }
# ar1 = 0 and rw=0 in standard analysis   
# ar1 = 1 and rw=0 when AR1 error moddel is considered. 
# ar1 = 0 and rw=1 when time-varying alpha is considered 
    mu[1] <-  fit[1] + ar1*phi * e0;
    cw[1] <- w[1]
    for(y in 2:nyrs){	   
      cw[y] <- cw[y-1]+ w[y]
      mu[y] <- fit[y] + rw*cw[y] + ar1*phi*e[y-1]
    }
    
#   Define Priors
    lnalpha ~ dunif(0,10)
    beta ~ dunif(0,10)
    sigma ~ dunif(0,10)
	sigmaw ~ dunif(0,10)
	sigmav ~ dunif(0,10)
    phi ~ dunif(-1,1)
    e0 ~ dnorm(0,0.001) 
    Tau <- 1/(sigma*sigma)
	tauw <- 1/(sigmaw*sigmaw)	
	tauv <- 1/(sigmav*sigmav)	
# Extract time-varyihg alapha
    for(y in 1:nyrs){
     lnalphai[y] <- lnalpha+cw[y] 
    }
# Likelihood 
    for(y in 1:nyrs){     
      R[y] ~ dlnorm(mu[y],Tau)
    }  
  }
  ar1 <- TRUE
  kf <- TRUE
    parameters <- c('lnalpha','beta','sigma') 
    if(ar1==TRUE){ ar1.p <- c('phi','e0')} else {ar1.p <- NULL}
    if(kf==TRUE){kf.p <- c('sigmaw','sigmav','lnalphai')} else {kf.p <- NULL}
    parameters.CR <- c(parameters, ar1.p,kf.p)
sim2 <- jags(data=datnew, parameters.to.save=parameters.CR, model.file= jag.model.CR,n.chains=1, 
	n.iter=10000,n.burnin=2000,n.thin=10,DIC=TRUE, working.directory=data_dir)	
  mcmc <- as.mcmc(sim2)
 
  plot(mcmc)
  
  


datnew<-list(nyrs=nyrs,S=S,R=R,d = d)
parameters.K <- c('lnalpha','beta','sigma','sigmav','sigmaw') 
 jag.model.K <- function(){
    for(y in 1:nyrs){
      s[y] <- S[y]/(10^d)
      fit[y] <- log(S[y]) + lnalpha[y] - beta * s[y] + v[y]
	  v[y] ~ dnorm(0,tauv)
     }
# ar1 = 0 and rw=0 in standard analysis   
# ar1 = 1 and rw=0 when AR1 error moddel is considered. 
# ar1 = 0 and rw=1 when time-varying alpha is considered 
    for(y in 2:nyrs){	   
       lnalpha[y] <- lnalpha[y-1] + w[y]
	   w[y] ~dnorm(0,tauw)
    }
    
#   Define Priors
    lnalpha[1] ~ dunif(0,10)
    beta ~ dunif(0,10)
    sigma ~ dunif(0,10)
    sigmav ~ dunif(0,10)
    sigmaw ~ dunif(0,10)
    tauv <- 1/(sigmav*sigmav)
	tauw <- 1/(sigmaw*sigmaw) 
    Tau <- 1/(sigma*sigma)
# Random walk for alpha

# Extract time-varyihg alapha

# Likelihood 
    for(y in 1:nyrs){     
      R[y] ~ dlnorm(fit[y],Tau)
    }  
  }


sim2 <- jags(data=datnew, parameters.to.save=parameters.K, model.file= jag.model.K,n.chains=1, 
	n.iter=100000,n.burnin=20000,n.thin=10,DIC=TRUE, working.directory=data_dir)	


sim

#---------------------------------------------------------------
#  AR1 Ricker 
#---------------------------------------------------------------
parameters.AR1 <- c('lnalpha','beta','phi','e0','sigma') 
jag.model.AR1 <- function(){
  for(y in 1:nyrs){
   s[y] <- S[y]/(10^d)
   fit[y] = log(S[y]) + lnalpha - beta * s[y]
   e[y] = log(R[y]) - fit[y]
      }
   mu[1] = fit[1] + phi * e0;	  
  for(y in 2:nyrs){	   
   mu[y] = fit[y] + phi * e[y-1]
   }
#     Define Priors
   lnalpha ~ dunif(0,10)
   beta ~ dunif(0,10)
   sigma ~ dunif(0,10)
   phi ~ dunif(-1,1)
   e0 ~ dnorm(0,0.001) 
   Tau <- 1/(sigma*sigma)
   tau.cor <- Tau*(1-phi*phi)
# Likelihood 
   for(y in 1:nyrs){     
     R[y] ~ dlnorm(mu[y],tau.cor)
	}  
}


parameters.AR12 <- c('lnalpha','beta','phi','sigma') 
jag.model.AR12 <- function(){
   s <- S/(10^d)  
  fit[1] = log(S[1]) + lnalpha - beta * s[1]
   mu[1] = fit[1]+ phi*e0
   e[1] = log(R[1]) - fit[1] - phi*e0
 for(y in 2:nyrs){   
   fit[y] = log(S[y]) + lnalpha - beta * s[y]
   mu[y] = fit[y] + phi*e[y-1]  
   e[y] <- log(R[y])-fit[y] - phi*e[y-1]
	}  
 for(y in 1:nyrs){   	
   R[y] ~ dlnorm(mu[y],tau.cor)	
 }  
#     Define Priors
   lnalpha ~ dunif(0,10)
   beta ~ dunif(0,10)
   sigma ~ dunif(0,10)
   phi ~ dunif(-1,1)
   e0 ~ dnorm(0,0.001) 
   Tau <- 1/(sigma*sigma)
   tau.cor <- Tau*(1-phi*phi)
# Likelihood 
  
}

#---------------------------------------------------------------
#  Classic Ricker  CRe
#---------------------------------------------------------------
parameters.CRe <- c('lnalpha','beta','sigma','sigma.x','sigma.y','Se','Re') 
jag.model.CRe <- function(){
  for(y in 1:nyrs){
   s[y] <- Se[y]/(10^d)
   lnRm[y] = log(Se[y]) + lnalpha - beta * s[y]
      }
#     Define Priors
   lnalpha ~ dunif(0,10)
   beta ~ dunif(0,10)
   sigma ~ dunif(0,10)
   phi ~ dunif(-1,1)
   Tau <- 1/(sigma*sigma)
   sigma.x ~ dunif(0,10)
   sigma.y  ~ dunif(0,10)
   Tau.x <- 1/(sigma.x*sigma.x)
   Tau.y <- 1/(sigma.y*sigma.y)  
   for(y in 1:nyrs){
    mu.x[y] ~ dnorm(12,0.01)
    Se[y] ~ dlnorm(mu.x[y],Tau.x)
   }
# Likelihood 
   for(y in 1:nyrs){  
     S[y] ~ dlnorm(log(Se[y]),Tau.x)   
     Re[y] ~ dlnorm(lnRm[y],Tau)
	 R[y] ~ dlnorm(log(Re[y]),Tau.y)
	}  
}
#---------------------------------------------------------------
#  Classic Ricker 
#---------------------------------------------------------------
parameters.CR <- c('lnalpha','beta','sigma') 
jag.model.CR <- function(){
  for(y in 1:nyrs){
   s[y] <- S[y]/(10^d)
   lnRm[y] = log(S[y]) + lnalpha - beta * s[y]
      }
#     Define Priors
   lnalpha ~ dunif(0,10)
   beta ~ dunif(0,10)
   sigma ~ dunif(0,10)
   phi ~ dunif(-1,1)
   Tau <- 1/(sigma*sigma)
# Likelihood 
   for(y in 1:nyrs){  
     R[y] ~ dlnorm(lnRm[y],Tau)
	}  
}
foo <- jags(data=datnew, parameters.to.save=parameters.AR1, model.file= jag.model.AR1,n.chains=1, 
	n.iter=100000,n.burnin=20000,n.thin=10,DIC=TRUE, working.directory=data_dir)	

	
#---------------------------------------------------------------
#  Ricker Brood Interaction 
#---------------------------------------------------------------
parameters.BI <- c('lnalpha','beta1','beta2','lnS0','sigma') 
jag.model.BI<- function(){
  for(y in 1:nyrs){
     s[y] <- S[y]/(10^d)
   lnRm1[y] <- log(S[y]) + lnalpha - beta1*s[y]
      }  
   lnRm[1] <- lnRm1[1] + beta2*exp(lnS0)/(10^d)   
  for(y in 2:nyrs){	   
   lnRm[y] <- lnRm1[y] + beta2*s[y-1]
   }
   
#    Define Priors
   lnalpha ~ dunif(0,10)
   beta1 ~ dunif(0,10)
   sigma ~ dunif(0,10)
   beta2 ~ dunif(-10,10) 
   lnS0 ~ dunif(0,16)
   Tau <- 1/(sigma*sigma)
# Likelihood 
   for(y in 1:nyrs){     
     R[y] ~ dlnorm(lnRm[y],Tau)
	}  
}

#---------------------------------------------------------------
#  Ricker Brood Interaction 
#---------------------------------------------------------------
parameters.BI2 <- c('lnalpha','beta3','lnS0','sigma') 
jag.model.BI2<- function(){
  for(y in 1:nyrs){
     s[y] <- S[y]/(10^d)
      }  
   lnRm[1] <- log(S[1]) + lnalpha - beta3*(s[1])*exp(lnS0)/(10^d)   
  for(y in 2:nyrs){	   
   lnRm[y] <- log(S[y]) + lnalpha - beta3*s[y]*s[y-1]
   }   
# Define Priors
   lnalpha ~ dunif(0,10)
   sigma ~ dunif(0,100)
   beta3 ~ dunif(-10,10)   
   lnS0 ~ dunif(0,16)
   Tau <- 1/(sigma*sigma)
# Likelihood 
   for(y in 1:nyrs){     
     R[y] ~ dlnorm(lnRm[y],Tau)
	}  
}

#---------------------------------------------------------------
#  Ricker Brood Interaction ER
#---------------------------------------------------------------
parameters.BI2e <- c('lnalpha','beta3','lnS0','sigma','sigma.x','sigma.y','Se','Re') 
jag.model.BI2e <- function(){
  for(y in 1:nyrs){
     s[y] <- Se[y]/(10^d)
      }  
   lnRm[1] <- log(Se[1]) + lnalpha - beta3*(s[1])*exp(lnS0)/(10^d)   
  for(y in 2:nyrs){	   
   lnRm[y] <- log(Se[y]) + lnalpha - beta3*s[y]*s[y-1]
   }   
# Define Priors
   lnalpha ~ dunif(0,10)
   sigma ~ dunif(0,100)
   beta3 ~ dunif(-10,10)   
   lnS0 ~ dunif(0,16)
   Tau <- 1/(sigma*sigma)
   sigma.x ~ dunif(0,10)
   sigma.y  ~ dunif(0,10)
   Tau.x <- 1/(sigma.x*sigma.x)
   Tau.y <- 1/(sigma.y*sigma.y)  
   for(y in 1:nyrs){
    mu.x[y] ~ dnorm(12,0.01)
    Se[y] ~ dlnorm(mu.x[y],Tau.x)
   }
   # Likelihood 
   for(y in 1:nyrs){     
     S[y] ~ dlnorm(log(Se[y]),Tau.x)   
     Re[y] ~ dlnorm(lnRm[y],Tau)
	 R[y] ~ dlnorm(log(Re[y]),Tau.y)
	}  
}
#sim <- jags(data=datnew, parameters.to.save=parameters.BI2e, model.file= jag.model.BI2e,n.chains=1, 
#	n.iter=100000,n.burnin=20000,n.thin=10,DIC=TRUE, working.directory=data_dir)	


#---------------------------------------------------------------
#  Ricker Brood Interaction 
#---------------------------------------------------------------
parameters.BI3 <- c('lnalpha','beta3','lnS0','sigma') 
jag.model.BI3<- function(){
  for(y in 1:nyrs){
     s[y] <- S[y]/(10^d)
      }  
  for(y in 1:(nyrs-1)){	   
   lnRm[y] <- log(S[y]) + lnalpha - beta3*s[y]*s[y+1]
   }
   lnRm[nyrs] <- log(S[nyrs]) + lnalpha - beta3*(s[nyrs])*exp(lnS0)/(10^d)      
#    Define Priors
   lnalpha ~ dunif(0,10)
   sigma ~ dunif(0,100)
   beta3 ~ dunif(-10,10)   
   lnS0 ~ dunif(0,16)
   Tau <- 1/(sigma*sigma)
# Likelihood 
   for(y in 1:nyrs){     
     R[y] ~ dlnorm(lnRm[y],Tau)
	}  
}

#---------------------------------------------------------------
#  Beverton Holt  
#---------------------------------------------------------------
parameters.BH <- c('lnalpha','beta','sigma') 
jag.model.BH <- function(){
  for(y in 1:nyrs){
   s[y] <- S[y]/(10^d)
   lnRm[y] <- lnalpha + log(S[y]) -log(1+beta*s[y])
      }
#     Define Priors
   lnalpha ~ dunif(0,10)
   beta ~ dunif(0,10)
   sigma ~ dunif(0,10)
   Tau <- 1/(sigma*sigma)
# Likelihood 
   for(y in 1:nyrs){     
     R[y] ~ dlnorm(lnRm[y],Tau)
	}  
}
#---------------------------------------------------------------
#  Deriso-Shunute  
#---------------------------------------------------------------
parameters.DS <- c('lnalpha','beta','c','sigma')
jag.model.DS <- function(){
  for(y in 1:nyrs){
     s[y] <- S[y]/(10^d)
   lnS[y] <- log(S[y])
   lnR[y] <- log(R[y])
   lnRm[y] = lnS[y] + lnalpha - log(1 + beta*c*s[y])/c 
      }
#     Define Priors
   lnalpha ~ dunif(0,10)
   beta ~ dunif(0,10)
   sigma ~ dunif(0,10)
   c ~ dunif(0,1)
   Tau <- 1/(sigma*sigma)
# Likelihood 
   for(y in 1:nyrs){     
     R[y] ~ dlnorm(lnRm[y],Tau)
	}  
}

#-------------------------------------------------------------------------------
# SR functions for post pcoessing 
#-------------------------------------------------------------------------------
# Classic Ricker
SR.CR <- function(lnalpha,beta,S,d){
     s <- S/(10^d)
   lnR <- log(S) + lnalpha - beta*s
   R <- exp(lnR)
   return(R)
}
# Ricker 
SR.BI <- function(lnalpha,beta1,beta2,S,d){
     s <- S/(10^d)
   lnR <- log(S) + lnalpha - (beta1-beta2)*s
   R <- exp(lnR)
   return(R)
}
# Ricker BI Carlson Model
SR.BI2 <- function(lnalpha,beta3,S,d){
     s <- S/(10^d)
   lnR <- log(S) + lnalpha - beta3*s^2
   R <- exp(lnR)
   return(R)
}
#  Beverton-Holt
SR.BH <- function(lnalpha,beta,S,d){
     s <- S/(10^d)
   lnR <- lnalpha +log(S) - log(1+beta*s)
   R <- exp(lnR)
   return(R)
}
#  Deriso-Shunute   
SR.DS <- function(lnalpha,beta,c,S,d){
   s <- S/(10^d)
   lnR <- log(S) + lnalpha - log(1 + beta*c*s)/c 
   R <- exp(lnR)
   return(R)
}

# Store Models  
nmodels <- 6
models <- list()
models$model1 = jag.model.CR
models$model2 = jag.model.AR1
models$model3 = jag.model.BI
models$model4 = jag.model.BH
models$model5 = jag.model.DS
models$model6 = jag.model.BI2

# Stor$e Model Parameters
parlist <- list()
parlist$par1 = parameters.CR
parlist$par2 = parameters.AR1
parlist$par3 = parameters.BI
parlist$par4 = parameters.BH
parlist$par5 = parameters.DS
parlist$par6 = parameters.BI2


################################################################
#  3.0: Run JAGS                                              
################################################################
#simlist <- list()
#for (i in 1:nmodels){
#sim <- jags(data=datnew, parameters.to.save=parlist[[i]], model.file= models[[i]],n.chains=1, 
#	n.iter=100000,n.burnin=20000,n.thin=10,DIC=TRUE, working.directory=data_dir)	
#simlist[[i]] <- sim
#}
#results_file <- 'Kasilof_JAGAR1.txt'
#dput(simlist,paste0(data_dir,results_file))

#sim2 <- jags(data=datnew, parameters.to.save=parameters.AR12, model.file= jag.model.AR12,n.chains=1, 
#	n.iter=100000,n.burnin=20000,n.thin=10,DIC=TRUE, working.directory=data_dir)	

################################################################
#  4.0: Read JAGS Results                                              
################################################################
simlist <- dget(paste0(data_dir,results_file))
simlist2 <- dget(paste0(data_dir,results_file2))

# Kasilof Range 
#S <- seq(0,1600000,length.out = 501)
# Escapement range
#eg <- c(140000,320000)


# Keni Range 
S <- seq(0,6000000,length.out = 1001)
# Escapement range
eg <- c(750000,1300000)

k <- 1000
main <- c('Classic Ricker','AR1 Ricker','Beverton-Holt','Deriso-Shnute','Ricker BYI','Ricker.BY2')
y <- c(1,2,4,5,3,6)

mdparlist <- list()
pred.Ym <-matrix(0,nrow=length(S),ncol=nmodels)
for (i in y){
sim <- simlist[[i]]  #read sim outputs
mcmc <- as.mcmc(sim)
post.samp <- as.matrix(mcmc)
mdparlist[[i]] <- data.frame(t(apply(post.samp, 2, FUN = median)))
pars <- mdparlist[[i]]
# Ricker CR and Ricker AR1
if(i ==1| i==2){
pred.Ym[,i] <- SR.CR(pars$lnalpha,pars$beta,S,d) 
}
# Brood interaction additive 
if (i == 3){
pred.Ym[,i] <- SR.BI(pars$lnalpha,pars$beta1,pars$beta2,S,d)
}
# Beverton Holt
if (i == 4){
pred.Ym[,i] <- SR.BH(pars$lnalpha,pars$beta,S,d)
}
# Derisso Shunte
if (i == 5){
pred.Ym[,i] <- SR.DS(pars$lnalpha,pars$beta,pars$c,S,d)
}
# Multiplicative Brood inteaction 
if (i == 6){
pred.Ym[,i] <- SR.BI2(pars$lnalpha,pars$beta3,S,d)
}
}

RS <- data.frame(S,pred.Ym)

SMSY <- matrix(0,nmodels,7)
for(i in 1:nmodels){
#MSY
SMSY[i,1] <- max(RS[,i+1]-RS[,1])
#SMSY
SMSY[i,2] <- RS[(RS[,i+1]-RS[,1])==SMSY[i,1],1]
foo <- RS[sign(RS[,i+1]-RS[,1]-0.9*SMSY[i,1])==1,1]
SMSY[i,3] <- min(foo)
SMSY[i,4] <- max(foo)
#UMSY
SMSY[i,5] <- SMSY[i,1]/(SMSY[i,1]+SMSY[i,2])
#SMAXR
SMSY[i,6] <- RS[RS[,i+1]==max(RS[,i+1]),1]
#SEQ
foo <- RS[sign(RS[,i+1]-RS[,1])==1,1]
SMSY[i,7] <- max(foo)
}
colnames(SMSY) <- c('MSY','Smsy','Smsy_L90','Smsy_U90','UMSY','Smax','SEQ')
SMSY



windowsFonts(
  A=windowsFont("Arial Black"),
  B=windowsFont("Bookman Old Style"),
  C=windowsFont("Comic Sans MS"),
  D=windowsFont("Times New Roman")
)
  
windows(record = TRUE) 
par(mfrow=c(nmodels,1),mar = c(2,1,1,1),oma = c(3,3,3,3),yaxs='i',xaxs='i',bty='l', family = 'D') 
gcols <- gray.colors(5, start = 0.1, end = 0.8, gamma = 2.2, alpha = NULL)
#---------------------------------------------------------------------------
#  Spawner - Recruit 
#---------------------------------------------------------------------------
for(i in 1:nmodels){
# SR points 
plot(Return~Spawner,data=dat1/k,pch=19,col=1,xlim=c(0,max(S)/k),ylim=c(0,max(dat1$Return)/k),main=main[i],xlab='',ylab='')
abline(0,1,lty=1)
# SR model lines 
lines(RS[,1]/k,RS[,y[i]+1]/k,lwd=2,lty = 1)
# Yield line 
lines(RS[,1]/k,(RS[,y[i]+1]-RS[,1])/k,lwd=1,lty = 2)
# SMSY
abline(v=SMSY[y[i],2]/k,lwd=2,lty = 1)
#points(Return~Spawner,data=dat1[dat1$Year>2009,]/k,pch=19,col=1,cex=1.5)
#text(Return~Spawner, data=dat1/k,labels=Year*k, cex= 0.7, pos=3)

}
########  Add Texts  ###########################################################
mult <- ifelse(k==1000000,paste0('(x million)'),ifelse(k>1,paste0('(x',k,')'),''))
mtext(paste('Recruits',mult), side = 2, line = 1, outer = TRUE)
mtext(paste("Escapement",mult), side = 1, line = 1, outer = TRUE)



#---------------------------------------------------------------------------
#  Model Specific Yields Curve 
#---------------------------------------------------------------------------
par(mfrow=c(1,1),mar = c(2,1,1,1),oma = c(3,3,3,3),yaxs='i',xaxs='i',bty='l', family = 'D') 
# Select Model 
i <- 1
dat1$Yield <- dat1$Return - dat1$Spawner
# Yield Curve 
plot(Yield~Spawner,data=dat1/k,pch=19,col=1,xlim=c(0,max(RS[,1])/k),ylim=c(min(0,min(dat1$Yield))/k,1.2*max(dat1$Yield)/k),main=main[i],xlab='',ylab='')
# Yield line 
lines(RS[,1]/k,(RS[,i+1]-RS[,1])/k,lwd=1.5,lty = 1)
lines(RS2[,1]/k,(RS2[,i+1]-RS[,1])/k,lwd=1.5,lty = 2)
abline(h = 0,lty=2)
mycol <- t_col('grey',50)
rect(eg[1]/k,min(0,min(dat1$Yield))/k,eg[2]/k,1.2*max(dat1$Yield)/k,col=mycol,border=0)
# SMSY
#points(Return~Spawner,data=dat1[dat1$Year>2009,]/k,pch=19,col=1,cex=1.5)
text((Return-Spawner)~Spawner, data=dat1/k,labels=Year*k, cex= 0.7, pos=3)
# Caluculate mean, range, n data 
dat2 <- dat1[with(dat1,Spawner>=eg[1] & Spawner<=eg[2]),]
eg.n <- length(dat2$Spawner)
eg.m <- mean(dat2$Yield)
eg.min <- min(dat2$Yield)
eg.max <- max(dat2$Yield)
tex <- c(paste0('Mean = ',round(eg.m,-3)/k,'k'),paste0('Range = ',round(eg.min,-3)/k,'k - ',format(round(eg.max,-3)/k,big.mark=","),'k'),paste('n =',eg.n))
#legend('topright',legend=tex,xjust=0,yjust=-1,box.col='white',bg='white',cex=1.3)
########  Add Texts  ###########################################################
mult <- ifelse(k==1000000,paste0('(x million)'),ifelse(k>1,paste0('(x',k,')'),''))
mtext(paste('Yield',mult), side = 2, line = 1, outer = TRUE)
mtext(paste("Escapement",mult), side = 1, line = 1, outer = TRUE)


#---------------------------------------------------------------------------
#  Model Specific SR Curve 
#---------------------------------------------------------------------------
par(mfrow=c(1,1),mar = c(2,1,1,1),oma = c(3,3,3,3),yaxs='i',xaxs='i',bty='l', family = 'D') 
# Select Model 
i <- 1
# Yield Curve 
plot(Return~Spawner,data=dat1/k,pch=19,col=1,xlim=c(0,max(RS[,1])/k),ylim=c(0,1.2*max(dat1$R)/k),main=main[i],xlab='',ylab='')
# Yield line 
lines(RS[,1]/k,(RS[,i+1])/k,lwd=1.5,lty = 1)
lines(RS2[,1]/k,(RS2[,i+1])/k,lwd=1.5,lty = 2)
abline(h = 0,lty=2)
mycol <- t_col('grey',50)
rect(eg[1]/k,min(0,min(dat1$R))/k,eg[2]/k,1.2*max(dat1$R)/k,col=mycol,border=0)
# SMSY
#points(Return~Spawner,data=dat1[dat1$Year>2009,]/k,pch=19,col=1,cex=1.5)
text((Return)~Spawner, data=dat1/k,labels=Year*k, cex= 0.7, pos=3)
abline(0,1)
# Caluculate mean, range, n data 
dat2 <- dat1[with(dat1,Spawner>=eg[1] & Spawner<=eg[2]),]
eg.n <- length(dat2$Spawner)
eg.m <- mean(dat2$Yield)
eg.min <- min(dat2$Yield)
eg.max <- max(dat2$Yield)
tex <- c(paste0('Mean = ',round(eg.m,-3)/k,'k'),paste0('Range = ',round(eg.min,-3)/k,'k - ',format(round(eg.max,-3)/k,big.mark=","),'k'),paste('n =',eg.n))
#legend('topright',legend=tex,xjust=0,yjust=-1,box.col='white',bg='white',cex=1.3)
########  Add Texts  ###########################################################
mult <- ifelse(k==1000000,paste0('(x million)'),ifelse(k>1,paste0('(x',k,')'),''))
mtext(paste('Recruit',mult), side = 2, line = 1, outer = TRUE)
mtext(paste("Escapement",mult), side = 1, line = 1, outer = TRUE)




#---------------------------------------------------------------------------
#  Smsy Profile Analyses 
#---------------------------------------------------------------------------
#windows(record = TRUE) 
#par(mfrow=c(nmodels,1),mar = c(2,2,2,2),oma = c(3,3,3,3),yaxs='i',xaxs='i',bty='l',family = 'D') 

percent.list <- c(0.9,0.85,0.8)
prof.list <- list()
for(i in 1:nmodels){
# Read data from JAGS output
sim <- simlist[[i]]
# Extract mcmc and create matrix
mcmc <- as.mcmc(sim)
post.samp <- as.matrix(mcmc)
# n is the number of simulation 
n <- dim(post.samp)[1]
# create temporary matrix
pred.Y <-  matrix(0,nrow = n,ncol=length(S))
print(i)
if(i ==1| i==2){
for(j in 1:n){pred.Y[j,] <- SR.CR(post.samp[j,'lnalpha'],post.samp[j,'beta'],S,d)} 
}
if (i == 3){
for(j in 1:n){pred.Y[j,] <- SR.BI(post.samp[j,'lnalpha'],post.samp[j,'beta1'],post.samp[j,'beta2'],S,d)} 
}
if (i == 4){
for(j in 1:n){pred.Y[j,] <- SR.BH(post.samp[j,'lnalpha'],post.samp[j,'beta'],S,d)} 
}
if (i == 5){
for(j in 1:n){pred.Y[j,] <- SR.DS(post.samp[j,'lnalpha'],post.samp[j,'beta'],post.samp[j,'c'],S,d)} 
}
if (i == 6){
for(j in 1:n){pred.Y[j,] <- SR.BI2(post.samp[j,'lnalpha'],post.samp[j,'beta3'],S,d)} 
}
# create temporary matrix
foo <-  matrix(0,nrow = n,ncol=length(S))
prof <- matrix(0,nrow = length(S),ncol=length(percent.list))
for(k in 1:length(percent.list)){
for(l in 1:n){foo[l,] <- ifelse((pred.Y[l,]-S)-percent.list[k]*max(pred.Y[l,]-S)>0,1,0)} 
# Calculate mean 
prof[,k] <- colMeans(foo)
}
print(head(prof))
prof.list[[i]] <- prof
}


# Profile function
plot.profile <- function(prof.data,eg,lineset,k,lw){
plot(prof.data[,1]/k,prof.data[,2],type='l',ylim=c(0,1),xlab=NULL,ylab=NULL,lwd=lw,lty=lineset[1])
ncols <- dim(prof.data)[2]
if(ncols>2){
for(i in 3:ncols){
lines(prof.data[,1]/k,prof.data[,i],lwd=lw,lty=lineset[i-1])
 }
 }
mycol <- t_col('grey',50)
rect(eg[1]/k,0,eg[2]/k,1,col=mycol,border=0)
}



lw<-1.5
k <- 1000
#windows(record = TRUE) 
#par(mfrow=c(nmodels,1),mar = c(2,2,2,2),oma = c(3,3,3,3),yaxs='i',xaxs='i',bty='l',family = 'D') 
for(i in 1:nmodels){
plot.profile(cbind(S,prof.list[[y[i]]]),eg,c(1,3,4),k,lw)
abline(h=0.9,lwd=lw)
abline(h=0.5,lwd=lw,lty=2)
legend('topright',col=1,lwd=lw,lty=c(1,3,4),legend =  c('90%','85%','80%'),
xjust=0,yjust=-1,box.col='white',bg='white')
title(main = main[i], sub = NULL, xlab = NULL, ylab = NULL,
      line = NA, outer = FALSE)
rug(x=seq(0,max(S),100000)/k,side =1,tick=-0.02)	  
}
########  Add Texts  ###########################################################
mult <- ifelse(k==1000000,paste0('(x million)'),ifelse(k>1,paste0('(x',k,')'),''))
mtext(paste('Probability'), side = 2, line = 1, outer = TRUE)
mtext(paste("Escapement",mult), side = 1, line = 1, outer = TRUE)





	  
plot.profile(mfoo.bi2,c(700000,1300000),c(1,3,4),1000,lw)
abline(h=0.9,lwd=lw)
abline(h=0.5,lwd=lw,lty=2)
legend('toprigh',col=1,lwd=lw,lty=c(1,3,4),legend =  c('90%','85%','80%'),
xjust=0,yjust=-1,box.col='white',bg='white')
title(main = 'Kenai Smsy Profile Ricker BYI2', sub = NULL, xlab = NULL, ylab = NULL,
      line = NA, outer = FALSE)
rug(x=seq(0,3000000,100000)/k,side =1,tick=-0.02)

plot.profile(mfoo.bh,c(700000,1300000),c(1,3,4),1000,lw)
abline(h=0.9,lwd=lw)
abline(h=0.5,lwd=lw,lty=2)
legend('toprigh',col=1,lwd=lw,lty=c(1,3,4),legend =  c('90%','85%','80%'),
xjust=0,yjust=-1,box.col='white',bg='white')
title(main = 'Kenai Smsy Profile BH', sub = NULL, xlab = NULL, ylab = NULL,
      line = NA, outer = FALSE)	 
rug(x=seq(0,3000000,100000)/k,side =1,tick=-0.02)
	  
plot.profile(mfoo.ds,c(700000,1300000),c(1,3,4),1000,lw)
abline(h=0.9,lwd=lw)
abline(h=0.5,lwd=lw,lty=2)
legend('toprigh',col=1,lwd=lw,lty=c(1,3,4),legend =  c('90%','85%','80%'),
xjust=0,yjust=-1,box.col='white',bg='white')
title(main = 'Kenai Smsy Profile DS', sub = NULL, xlab = NULL, ylab = NULL,
      line = NA, outer = FALSE)	 
rug(x=seq(0,3000000,100000)/k,side =1,tick=-0.02)
	  	  




#---------------------------------------------------------------------------
#  Time Series 
#---------------------------------------------------------------------------
k <- 1000000
ltype = c(1,2,4,5,6)
plot(Return~Years,data=dat1/k,pch=19,col=1,xlim=c(0,max(S)/k),ylim=c(0,max(dat1$Return)/k),xlab='',ylab='')
abline(0,1,lty=3)
for (i in 1:5){
lines(RS[,1]/k,(RS[,i+1]-RS[,1])/k,lwd=2,lty = ltype[i],col=i)
abline(v=Smsy[i]/k,lwd=2,col=i)
}
text((Return-Spawners)~Spawners, data=dat1/k,labels=Year*k, cex= 0.8, pos=3)
legend('toprigh',col=c(1:5),lwd=2,lty=c(1,2,4,5,6),legend = c('Ricker','Ricker AR1','Ricker BI','Beverton-Holt','Derisso-Shunute'),
bty='o',xjust=0,yjust=-1,bg='white')

########  Add Texts  ###########################################################
mult <- ifelse(k==1000000,paste0('(x million)'),ifelse(k>1,paste0('(x',k,')'),''))
mtext(paste('Return',mult), side = 2, line = 1, outer = TRUE)
mtext(paste("Escapement",mult), side = 1, line = 1, outer = TRUE)


#horse tail 

plot(S/k,SR.BI2(post.samp[1,'lnalpha'],post.samp[1,'beta3'],S,d)/k,type='l',ylim=c(0,1.2*max(dat1$Return)/k))
for(i in 1:8000){
lines(S/k,SR.BI2(post.samp[i,'lnalpha'],post.samp[i,'beta3'],S,d)/k,type='l')
}


