#####################################################
## \author Konstantinos Lagogiannis 2019
## 
## \short Figure 3 Estimating Hunt efficiency
##        Sample script that models labelled hunt events to identify statistical differences in hunt success and hunt rate and consumption
##############################################

source("common_lib.R")
library("MASS") ##For contours
library("rjags")

## The mixture Of Poisson Drawing from Gammaa, gives a negative binomial
## ID identifies the group ID of each larva
modelNBinom="model {
         
         for(i in 1:3) {
            ### Hunt Rate Model
             r[i] ~ dgamma(1,1) ##
             q[i] ~ dunif(0.0,1)# 

             ##Success/Fail Prob Estimation 
             f[i] ~ dbeta(1,1) ##Prob of capture Fail
             p[i] ~ dbeta(1,1) ##Prob of capture success
             t[i] ~ dbeta(1,1) ##Prob Of Enganging With Prey Given HuntMode Is On
         }

        
         for(i in 1:NTOT){
             ##Here we can extend this model to infer the probability to reaching final capture once a hunt event is initiated 
             Events[i] ~  dnegbin(q[ID[i]],r[ID[i]] ) ##The data on total number of hunt event for a larva (Can also be used for hunt events that reached capture in TrackPrey)
             TrackPrey[i] ~ dbinom(t[ID[i]],Events[i]) ##Number of hunt events that resulted in attempt at prey
             Success[i] ~ dbinom(p[ID[i]],TrackPrey[i]) ##Subset of Successful capture hunt events
             Fail[i] ~ dbinom(f[ID[i]],TrackPrey[i]) ##Reciprocal, Failed Capture events
         }
}"



message(paste(" Loading Hunt Event List to Analyse... "))
##Load From Central Function
datHuntLabelledEventsSB <- readRDS(file=paste0("dat/datHuntLabelledEvents.rds"))

##Select subset of  evoked hunt events only
datHuntLabelledEventsSB_evoked <-  datHuntLabelledEventsSB[datHuntLabelledEventsSB$groupID %in% c("LL","DL","NL"),]
datFishSuccessRate <- getHuntSuccessPerFish(datHuntLabelledEventsSB_evoked)

tblResSB <- table(convertToScoreLabel(datHuntLabelledEventsSB_evoked$huntScore),datHuntLabelledEventsSB_evoked$groupID)

datFishSuccessRate$groupID <- factor(datFishSuccessRate$groupID)
strGroups <-levels(datFishSuccessRate$groupID)

## Get the number larvae from each of the evoked test conditions
NRecCount_DL <- table(datFishSuccessRate$groupID)["DL"]
NRecCount_NL <- table(datFishSuccessRate$groupID)["NL"]
NRecCount_LL <- table(datFishSuccessRate$groupID)["LL"]

##Make data structure for jags model
datatest=list(Success=datFishSuccessRate$Success,
              Fail=datFishSuccessRate$Fails,
              TrackPrey=datFishSuccessRate$Fails+datFishSuccessRate$Success, ##Number of Prey Engangements
              Events=datFishSuccessRate$CaptureEvents, ######Includes No_Target - ie cases where stimuli Trigger Fish HuntMode But no Prey Tracking Seems to take place
              ID=as.numeric(datFishSuccessRate$groupID),
              NTOT=nrow(datFishSuccessRate));

varnames1=c("q","p","t","f","r")
burn_in=1000; ##Initial N of steps to establish model
steps=20000; ##Number of samples to obtain
thin=10;
chains=3 ## Parallel Sample Chains

##Setup the Jags Model and run sampler
m=jags.model(textConnection(modelNBinom),data=datatest,n.chains=chains);
update(m,burn_in)
draw=jags.samples(m,steps,thin=thin,variable.names=varnames1)


### Draw Distribution oF Hunt Rates - 
## for the exp draw (z= p/(1-p)) ## But it is the same for Rate Of Gamma Too / Or inverse for scale
plotsamples <- 1500
schain <-1:3
Range_ylim <- c(1,25)
cex <- 1.2

### It looks like Tha Gamma scale is the inverse:  gamma rate!
##13/11/19 Spotted difference in GammaRate q/(1-q) was inverse - now matches HuntRateInPreyRange model analysis
HEventHuntGammaRate_LL <-tail(draw$q[2,,schain],plotsamples)/(1-tail(draw$q[2,,schain],plotsamples));
HEventHuntGammaRate_DL <- tail(draw$q[1,,schain],plotsamples)/(1-tail(draw$q[1,,schain],plotsamples));      
HEventHuntGammaRate_NL <- tail(draw$q[3,,schain],plotsamples)/(1-tail(draw$q[3,,schain],plotsamples));      
HEventHuntGammaShape_LL <- tail(draw$r[2,,schain],plotsamples)
HEventHuntGammaShape_DL <- tail(draw$r[1,,schain],plotsamples)
HEventHuntGammaShape_NL <- tail(draw$r[3,,schain],plotsamples)

HEventSuccess_LL <-tail(draw$p[2,,schain],plotsamples)
HEventSuccess_DL <-tail(draw$p[1,,schain],plotsamples)
HEventSuccess_NL <-tail(draw$p[3,,schain],plotsamples)

###Mean Rates As Exp OF Gamma
MeanHuntRate_LL <- HEventHuntGammaShape_LL*1/HEventHuntGammaRate_LL
MeanHuntRate_DL <- HEventHuntGammaShape_DL*1/HEventHuntGammaRate_DL
MeanHuntRate_NL <- HEventHuntGammaShape_NL*1/HEventHuntGammaRate_NL

HConsumptionRate_NL <- HEventSuccess_NL*MeanHuntRate_NL
HConsumptionRate_LL <- HEventSuccess_LL*MeanHuntRate_LL
HConsumptionRate_DL <- HEventSuccess_DL*MeanHuntRate_DL

#'' Measure correlation between success probability and efficiency - 
chain <-1 ##chose  chain for plotting samples
HCovRateAndEfficiency_NL <- cor( tail(HEventSuccess_NL[,chain],plotsamples),tail(MeanHuntRate_NL[,chain],plotsamples) )
HCovRateAndEfficiency_LL <- cor( tail(HEventSuccess_LL[,chain],plotsamples),tail(MeanHuntRate_LL[,chain],plotsamples) )
HCovRateAndEfficiency_DL <- cor( tail(HEventSuccess_DL[,chain],plotsamples),tail(MeanHuntRate_DL[,chain],plotsamples) )


## Fig 3 The Efficiency Hunt Rate plot ##
plotWidthIn <- 8
  outer = FALSE
  line <- 2.6 ## SubFig Label Params
  lineGroupLabel <- line - 32 ##pie chart group label
  cex = 1.4
  adj  = 0.5
  padj <- -0
  las <- 1

  ##Margin: (Bottom,Left,Top,Right )
  par(mar = c(4.2,4.8,1.1,1))
  #layout(matrix(c(1,,2), 1, 2, byrow = TRUE))
  
####### Efficiency Inference Plot ## Taken From stat_SyccessVsFailModel.r ####
  nlevels <- 5
  zLL <- kde2d(c(HEventSuccess_LL[,schain]), c(MeanHuntRate_LL[,schain]),n=80)
  zNL <-  kde2d(c(HEventSuccess_NL[,schain]), c(MeanHuntRate_NL[,schain]),n=80)
  zDL <-  kde2d(c(HEventSuccess_DL[,schain]), c(MeanHuntRate_DL[,schain]),n=80)
  
  plot(HEventSuccess_DL, MeanHuntRate_DL,col=colourHPoint[3],ylim=Range_ylim,xlim=c(0,0.5),pch=pchL[3],
       main=NA, #"Bayesian Estimation for Hunt Rate and Efficiency",
       xlab=NA,#"Probability of Success q",
       ylab=NA,cex.main =cex,cex.axis=1.5 )#(expression(paste("Hunt Rate ",lambda ) ) )  ) #paste("Hunt Rate", )
  points(HEventSuccess_LL, MeanHuntRate_LL,col=colourHPoint[2],ylim=Range_ylim,xlim=c(0.1,0.5),pch=pchL[2])
  points(HEventSuccess_NL, MeanHuntRate_NL,col=colourHPoint[1],ylim=Range_ylim,xlim=c(0.1,0.5),pch=pchL[1])
  mtext(side = 1, cex=cex, line = line,expression(paste("Probability of success (",q,")" ) )  ) 
  mtext(side = 2, cex=cex, line = line, expression(paste("Estimated capture attempts / 10min" ) )  ) #(",lambda,")"
  #mtext("B",at="topleft",outer=F,side=2,col="black",font=2,las=las,line=4,padj=-11,adj=0,cex=cex,cex.main=4)
  
  contour(zDL, drawlabels=FALSE, nlevels=nlevels,add=TRUE)
  contour(zLL, drawlabels=FALSE, nlevels=nlevels,add=TRUE)
  
  contour(zNL, drawlabels=FALSE, nlevels=nlevels,add=TRUE)
  #legend("topright", legend=paste(strGroups," n=",c(NRecCount_DL,NRecCount_LL,NRecCount_NL)),fill=colourL)
  legend("topright",cex=cex+0.2,
         legend=c(  expression (),
                    bquote(NF[""] ~ '#' ~ .(NRecCount_NL)  ),
                    bquote(LF[""] ~ '#' ~ .(NRecCount_LL)  ),
                    bquote(DF[""] ~ '#' ~ .(NRecCount_DL)  )  ), #paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
         pch=pchL, col=colourLegL)

 
####### Consumption Estimates ####
  par(mar = c(4.2,4.7,1.1,1))
  #### Consumption
  plot(density(HConsumptionRate_NL),xlim=c(0,8),ylim=c(0,1.5),col=colourLegL[1],lwd=4,lty=1
       ,cex.main =cex,cex.axis=1.5, xlab=NA,ylab=NA,main=NA) #"Mean consumption per larva"
  lines(density(HConsumptionRate_LL),col=colourLegL[2],lwd=4,lty=2)
  lines(density(HConsumptionRate_DL),col=colourLegL[3],lwd=4,lty=3)
  legend("topright",cex=cex+0.2,
         legend=c(  expression (),
                    bquote(NF[""] ~ '#' ~ .(NRecCount_NL)  ),
                    bquote(LF[""] ~ '#' ~ .(NRecCount_LL)  ),
                    bquote(DF[""] ~ '#' ~ .(NRecCount_DL)  )  ), ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
         col=colourLegL,lty=c(1,2,3),lwd=4)
  mtext(side = 1, cex=cex, line = line, expression(paste("Estimated consumption (Prey/10min)  ") ))
  mtext(side = 2, cex=cex, line = line, expression("Density function") )



  #### Plot Hunt Power ####
  par(mar = c(4.2,4.7,1.1,1))
  plotHuntPowerDataCDF(datHuntLabelledEventsSB)
  

  #### Plot Hunt Efficiency ####
  par(mar = c(4.2,4.7,1.1,1))
  plotHuntEfficiencyDataCDF(datHuntLabelledEventsSB)


######################

## Z-Score ###
muEff_LF <- mean(datFishSuccessRate[datFishSuccessRate$groupID == "LL","Efficiency" ], na.rm = TRUE)
muEff_NF <- mean(datFishSuccessRate[datFishSuccessRate$groupID == "NL","Efficiency" ], na.rm = TRUE)
muEff_DF <- mean(datFishSuccessRate[datFishSuccessRate$groupID == "DL","Efficiency" ], na.rm = TRUE)
  
stdNFLF <- sd(datFishSuccessRate[datFishSuccessRate$groupID %in% c("NL","LL"),"Efficiency" ], na.rm = TRUE)
stdDFLF <- sd(datFishSuccessRate[datFishSuccessRate$groupID %in% c("DL","LL"),"Efficiency" ], na.rm = TRUE)

zScoreNFLF <- (muEff_LF-muEff_NF)/stdNFLF
zScoreDFLF <- (muEff_LF-muEff_DF)/stdDFLF

message(paste("Z-Score NF-LF Efficiency",zScoreNFLF) )
message(paste("Z-Score DF-LF Efficiency",zScoreDFLF) )