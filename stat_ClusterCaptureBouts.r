### \author Kostas Lagogiannis 2019-04-17 
## Figure 4 of manuscript - Clustering of capture bouts 
## This is a mixture of two Gaussians Model for Clustering Capture speed (fast/ slow), based on Speed & Distance from Prey
## The covariances also of the model also show the relationship between the capture bout speed and distance to prey
## Points are assigned to the fast cluster if their posterior prob of occupying that cluster is above  0.7 (minClusterLikelyhood) (see addClusterColumns function)


library(rjags)
library(runjags)


library(ggplot2) ##install.packages("ggplot2")
library(ggExtra)##  install.packages("ggExtra") ##devtools::install_github("daattali/ggExtra").
library(ggpubr) ##install.packages("ggpubr")


source("common_lib.R")

strmodel_capspeedVsDistance <- "
var x_rand[2,2];

model {

##Draw capt speed from 2d gaussian
for (i in 1:N)
{
  ##Draw from gaussian model  as determined by mod flag
  c[i,1:2] ~ dmnorm(mu[mID[i]+1,],prec[mID[i]+1, , ]) ## data in column 1 and 2
  mID[i] ~ dbern(0.5) ##Se Gaussian class membership randomly
  
}

## ????XXXX Fit Bernouli distribution on Number of Hunt |Events that have a high-speed strike ??
## We used  a normal for Probability of Strike Swim 
pS  ~ dnorm(sum(mID)/N,1000)T(0,1)
mStrikeCount ~ dbin(pS,N )

##Covariance matrix and its inverse -> the precision matrix
## for each Gaussian in the mixture (1 and 2)
for  (g in 1:2)
{
  prec[g,1:2,1:2] <- inverse(cov[g,,])
  
  cov[g,1,1] <- sigma[g,1]*sigma[g,1]
  cov[g,1,2] <- sigma[g,1]*sigma[g,2]*rho[g]
  cov[g,2,1] <- cov[g,1,2]
  cov[g,2,2] <- sigma[g,2]*sigma[g,2]
  
  ## Priors 
  sigma[g,1] ~ dunif(0,1) ##dist prey - Keep it broad within the expected limits 
  
  rho[g] ~ dunif(-1,1) ##The covar coefficient
}
  ## Low Speed Captcha cluster
  mu[1,1] ~ dnorm(0.5,0.01)T(0.0,) ##Distance prey
  mu[1,2] ~ dnorm(5,1)T(0,) ##cap speed
  sigma[1,2] ~ dunif(0,2) ##the low cap speed sigma 

  ## High speed Capture Cluster
  mu[2,1] ~ dnorm(0.5,0.01)T(0.0,) ##Distance prey ##precision=1/sigma^2
  mu[2,2] ~ dnorm(35,1)T(mu[1,2],) ##cap speed
  sigma[2,2] ~ dunif(0,10) ##the high cap speed sigma 

## Synthesize data from the distribution
x_rand[1,] ~ dmnorm(mu[1,],prec[1,,])
x_rand[2,] ~ dmnorm(mu[2,],prec[2,,])

} "


##  Init  datastruct that we pass to model ##
steps <- 5500
thin <- 2
chains = 3
str_vars <- c("mu","rho","sigma","cov","x_rand","mID","mStrikeCount","pS")
##load  behavioural data of each group 
ldata_LF <- readRDS(file=paste0("dat/huntEpisodeDataMergedWithLarvalSuccess_LF.rds") )
ldata_NF <- readRDS(file=paste0("dat/huntEpisodeDataMergedWithLarvalSuccess_NF.rds") )
ldata_DF <-  readRDS(file=paste0("dat/huntEpisodeDataMergedWithLarvalSuccess_DF.rds") )

##Convert To convenient format for selecting Columns
distVsSpeed_LF <- data.frame(unlist(ldata_LF$c))
distVsSpeed_NF <- data.frame(unlist(ldata_NF$c))
distVsSpeed_DF <- data.frame(unlist(ldata_DF$c))

### RUN MODEL on Each Group independently###
jags_model_LF <- jags.model(textConnection(strmodel_capspeedVsDistance),
                            data = list(N=NROW(distVsSpeed_LF),
                                    c=cbind(dist=distVsSpeed_LF$DistanceToPrey,speed=distVsSpeed_LF$CaptureSpeed)), 
                            n.adapt = 500, n.chains = chains, quiet = F)
update(jags_model_LF, 500)
draw_LF=jags.samples(jags_model_LF,steps,thin=thin,variable.names=str_vars)

## Not Fed
jags_model_NF <- jags.model(textConnection(strmodel_capspeedVsDistance), 
                            data = list(N=NROW(distVsSpeed_NF),
                                        c=cbind(dist=distVsSpeed_NF$DistanceToPrey,speed=distVsSpeed_NF$CaptureSpeed)), 
                            n.adapt = 500, n.chains = 3, quiet = F)
update(jags_model_NF)
draw_NF=jags.samples(jags_model_NF,steps,thin=2,variable.names=str_vars)

##  DRY  Fed
jags_model_DF <- jags.model(textConnection(strmodel_capspeedVsDistance), 
                            data = list(N=NROW(distVsSpeed_DF),
                                        c=cbind(dist=distVsSpeed_DF$DistanceToPrey,speed=distVsSpeed_DF$CaptureSpeed)), 
                            n.adapt = 500, n.chains = 3, quiet = F)
update(jags_model_DF, 500)
draw_DF=jags.samples(jags_model_DF,steps,thin=2,variable.names=str_vars)

## Plot The gaussian cluster models and data coloured according to fast/slow cluster membershipt - As in Fig 4 of manuscript
  plotClusteredData(distVsSpeed_NF,draw_NF)
  plotClusteredData(distVsSpeed_LF,draw_LF)
  plotClusteredData(distVsSpeed_DF,draw_DF)

## Extra validation plot
## Check Clustering Model and how they split the Distribution of  Capture Speeds of each group##
  par(mar = c(3.9,4.3,1,1))
  layout(matrix(c(1,2,3),1,3, byrow = FALSE))
  npchain<-3
  plotCaptureSpeedFit(distVsSpeed_NF,draw_NF,1,npchain)
  #title(main="Model capture Speed")
  plotCaptureSpeedFit(distVsSpeed_LF,draw_LF,2,npchain)
  plotCaptureSpeedFit(distVsSpeed_DF,draw_DF,3,npchain)
  
############### END ###


