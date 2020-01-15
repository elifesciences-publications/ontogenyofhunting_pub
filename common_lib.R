##

### \notes main Tracking handles the initial import of tracked videos into data frames, and the detection of hunting episodes - 
## Hunting episdes are placed in as "setn15-HuntEvents-SB-Updated-Merged2 - and labelling updates this file.
##
### Hunt event Motion Analysis 
##  For close analysis of hunting the DataLabelling directory has a script to draw hunt events of as pecific assigned label or unlabelled ones
##  The output files of fish and prey (food)  are placed in HuntEventS_Retracked and from there they need to be moved to the subdirectory where LF, DF, NF 
##  retracked (Success / Fail separated) events are placed/ 
##  To import them run the HuntEpisodeAnalysis/runimportHuntEpisodeTrackFiles.r - This will merge food and larva motion records and produce
##  an imported hunt event register. From there analysis and a plot for each hunt episode is produced by running runHuntEpisodeAnalysis.r
#
## A script that Randomly and blindly allows to manually label hunt events exists in the DataLablling dir, main_LabellingBlind - this updates the 
## detected huntevent register ("setn15-HuntEvents-SB-Updated-Merged2)
#####################
#library(MASS)

DIM_PXRADIUS <- 790 #Is the Radius Of the dish In the Video
DIM_MMPERPX <- 35/DIM_PXRADIUS ##35mm Opening of The viewport Assumed
DIM_DISTTOMOUTH_PX <- 14 ## Estimated Distance from CEntroid To Mouth based on head template size used in tracker
DIM_DISTTOMOUTH_MM <- DIM_DISTTOMOUTH_PX*DIM_MMPERPX ## Estimated Distance from CEntroid To Mouth based on head template size used in tracker
G_APPROXFPS              <- 410
G_THRESHUNTANGLE         <- 20 #Define Min Angle Both Eyes need for a hunting event to be assumed
G_THRESHUNTVERGENCEANGLE <- 45 ## When Eyes pointing Inwards Their Vergence (L-R)needs to exceed this value for Hunting To be considered
G_THRESHCLIPEYEDATA      <- 40 ##Limit To Which Eye Angle Data is filtered to lie within
G_MINGAPBETWEENEPISODES  <- 300
G_MINEPISODEDURATION     <- 100
G_MIN_BOUTSPEED          <- 0.2 ##mm/frame - Need to be above to be considered A Motion Bout
G_THRES_CAPTURE_SPEED    <-  16 ###15##mm/sec ##Theshold on Body Speed above which a hunt event is marked to have a capture strike
G_THRES_MOTION_BOUT_SPEED <- 2.9 ##Got from Clustering #4 ##mm/sec
PREY_COUNT_FRAMEWINDOW   <- 1600 ##Number oF Frames Over which to count Prey Stats at Beginning And End Of Experiments
G_MIN_TURNBOUT_ANGLE     <- 10 ##
G_THRES_TAILFQ_BOUT      <- 9.5 ##Hz
nFrWidth                 <- 20 ## Sliding Window Filter Width - Reduced From 50 to 20 to improve Meanf sliding window speed estimation lags
nEyeFilterWidth          <- nFrWidth*2
MIN_BOUT_DURATION        <- 10 ##Used in HuntEpisodeAnalysis_lib
MIN_BOUT_PAUSE           <- 25
G_MIN_BOUTSCORE          <- 1 ##Number of coincident event needed to detect bout in frame - (a combo of Tail Fq, Centroid Speed,Turn speed)

## Plot Options ##
FONTSZ_AXISLAB <- 1.2
FONTSZ_AXIS    <- 1.2 ##Axis Tik Number Size

line = 2.8 ## SubFig Label Params
lineAxis = 2.7
lineTitle = 2.7
lineXAxis = 3.0
cex = 1.4
adj  = 1.0
padj <- -8.0
las <- 1

##For the 3 Groups 
###  NF, LF, DF , Black Colouring 

colourLegE <- c("#FB9A99C8","#B2DF8AC8","#A6CEE3C8","#FDBF6F23") 
colourLegL <- c("#E31A1CD3", "#33A02CD3", "#1F78B4D3", "#FF7F00D3") 
colourL    <-c("#03B303AF","#E60303AF","#0303E6AF")
colourHPoint <- c("#E31A1C37","#33A02C37","#1F78B437","#FB9A9937","#B2DF8A37","#A6CEE337","#FDBF6F37")
colourHLine <- colourLegL

colourH <- c("#0303E663","#03B30363","#E6030363")
colourD <- c("#0303E623","#03B30323","#E6030323")
colourL <- c("#0303E6AF","#03B303AF","#E60303AF")


pchL <- c(1,2,0,17,15,16,4) ## The style of bullet used for each group DL, LL, NL
lineTypeL <- c(2,1,3,4) ## The style of bullet used for each group DL, LL, NL

## Condition Labels
strDataLabels <- expression("NF-s","LF-s","DF-s","NF-e","LF-e","DF-e" )
strGroupID <- c("DL","LL","NL") 
## GLOBAL VARS ###


########### Capture Type Clustering Functions ####
##################################################

## Plots the Data Density and the 2 Gaussians fititng high and low speed capture swims
plotCaptureSpeedFit <- function(datSpeed,drawMCMC,colourIdx,nchain = 1)
{
  xquant <- seq(0,70,1)
  XLIM <- c(0,60)
  YLIM <- c(0,0.15)
  pdistBW <- 5 ## mm/sec
  strKern <- "gaussian"
  #ntail <- NROW(drawMCMC$mu[1,2,,nchain])*0.10
  ntail <- min(50,NROW(drawMCMC$mu[1,1,,1])*0.10)
  
  plot(density(datSpeed$CaptureSpeed,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=XLIM,ylim=YLIM,cex=cex,cex.axis=cex 
       ,main=NA,xlab = NA,ylab=NA)
  for (i in 1:(ntail-1) )
  {
    lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[1,2,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[1,2,ntail-i,nchain],1)),type='l',col=colourHLine[colourIdx],lty=1 )
    lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[2,2,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[2,2,ntail-i,nchain],1)),type='l',col=colourHLine[colourIdx],lty=2 )
  }
  
  dens<- density(datSpeed$CaptureSpeed,bw=pdistBW,kernel=strKern)
  lines(dens,col="black",lwd=4,xlim=XLIM )
  legend("topright",title="",cex=cex,
         legend=c( paste0("",dens$n, "# Data density "), #(Bw:",prettyNum(digits=2, pdistBW ),")" ) ,
                   paste("Model low speed " ),
                   paste("Model high speed " )),
         col=c("black",colourLegL[colourIdx],colourLegL[colourIdx]),lwd=c(3,1,1),lty=c(1,1,2) ) 
  
  mtext(side = 1,cex=cex, line = 3.2, expression("Capture speed (mm/sec) " ))
  mtext(side = 2,cex=cex, line = 2.5, expression("Density function " ))
  
}



###Used for drawing contour in ggplot -
## Draw the model fit above the cluster points
getFastClusterGrid <- function(drawMCMC)
{
  ### Add the cluster contours ###
  xran <- seq(0,0.8,0.05) ##Distance Grid
  yran <- seq(0,70,1) ##Speed Grid
  
  ##Example Code for PLotting Inferred Multivariate Cluster
  nsteps <- NROW(drawMCMC$mID[,,1][1,])
  mat_cov_fast <- rowMeans(drawMCMC$cov[2,,,(nsteps-1000):nsteps,1],dim=2) ##Average over samples
  
  mat_mu <- rowMeans(drawMCMC$mu[,,(nsteps-1000):nsteps,1],dim=2)
  #valGrid <- matrix( expand.grid(distance=xran,speed=yran),nrow=NROW(xran),ncol=NROW(yran) )
  valGrid <- expand.grid(distance=xran,speed=yran)
  #matrix(valGrid$distance,ncol=NROW(xran))
  cluster_z <- mvtnorm::dmvnorm(valGrid,mean=mat_mu[2,],sigma=mat_cov_fast )
  cluster_fast <- cbind(valGrid,cluster_z, factor(rep(16,NROW(valGrid) ),levels=c(1,16),labels=c("slow","fast")  ) )
  
  names(cluster_fast) <- c("DistanceToPrey", "CaptureSpeed", "Density","Cluster")
  
  return(cluster_fast)
}

###Used for drawing contour in ggplot
getSlowClusterGrid <- function(drawMCMC)
{
  ### Add the cluster contours ###
  xran <- seq(0,0.8,0.05) ##Distance Grid
  yran <- seq(0,70,1) ##Speed Grid
  
  ##Example Code for PLotting Inferred Multivariate Cluster
  nsteps <- NROW(drawMCMC$mID[,,1][1,])
  mat_cov_slow <- rowMeans(drawMCMC$cov[1,,,(nsteps-1000):nsteps,1],dim=2) ##Average over samples
  
  mat_mu <- rowMeans(drawMCMC$mu[,,(nsteps-1000):nsteps,1],dim=2)
  #valGrid <- matrix( expand.grid(distance=xran,speed=yran),nrow=NROW(xran),ncol=NROW(yran) )
  valGrid <- expand.grid(distance=xran,speed=yran)
  #matrix(valGrid$distance,ncol=NROW(xran))
  cluster_z <- mvtnorm::dmvnorm(valGrid,mean=mat_mu[1,],sigma=mat_cov_slow )
  cluster_slow <- cbind(valGrid,cluster_z, factor(rep(1,NROW(valGrid) ),levels=c(1,16),labels=c("slow","fast")  ) )
  
  names(cluster_slow) <- c("DistanceToPrey", "CaptureSpeed", "Density","Cluster")
  
  return(cluster_slow)
}




## Uses the sampled occupancy for cluster membership, assigns cluster depending on whether the 
# average time spend in fast cluster (mID=1) is above minClusterLikelyhood
addClusterColumns <- function(datDistVsSpeed,draw_F)
{
  ## Set Cluster
  minClusterLikelyhood <- 0.7
  nsamples <- 200
  ch <- 1
  tsamples <- NROW(draw_F$mID[1,,ch]) ##Get Total Number of Samples in Chain
  
  
  ###Add Cluster Column / Set According to mean occupancy being above 0.7
  datDistVsSpeed$fastClustScore <- apply(draw_F$mID[,(tsamples-nsamples):tsamples,ch],1,mean) 
  datDistVsSpeed$pchL <- 1
  datDistVsSpeed$pchL[datDistVsSpeed$fastClustScore > minClusterLikelyhood] <- 16 ##Set marker Point Type used for display
  datDistVsSpeed$Cluster <- factor(labels=c("slow","fast"),datDistVsSpeed$pchL)  ##Define Cluster type factor
  return (datDistVsSpeed)
}
### Capture Type Clusting plot as in Figure 4 function to show the inferred Gaussian models for capture type clustering
## the data points coloured according to cluster membership along with their densities on the plot margins 
## \note Uses ggplot
plotClusteredData <- function(distVsSpeed_F,draw_F)
{
  distVsSpeed_F <- addClusterColumns(distVsSpeed_F,draw_F)
  
  p_NF = ggplot( distVsSpeed_F, aes(DistanceToPrey, CaptureSpeed,color =Cluster,fill=Cluster))  +
    ggtitle(NULL) +
    theme(axis.title =  element_text(family="Helvetica",face="bold", size=16),
          axis.text = element_text(family="Helvetica",face="bold", size=16),
          plot.margin = unit(c(1,1,1,1), "mm"), legend.position = "none") +
    fill_palette("jco")
  
  p_NF = p_NF + geom_point( size = 3, alpha = 0.6,aes(color =distVsSpeed_F$Cluster) ) +  xlim(0, 0.8) +  ylim(0, 80) +
    scale_color_manual( values = c("#00AFBB", "#E7B800", "#FC4E07") )
  
  contour_fast <- getFastClusterGrid(draw_F) ## Draw the mvtnorm model fit contour
  contour_slow <- getSlowClusterGrid(draw_F)
  p_NF = p_NF +
    geom_contour(contour_fast, mapping = aes(x = DistanceToPrey, y = CaptureSpeed, z = Density) ,linetype=2 ) +
    geom_contour(contour_slow, mapping = aes(x = DistanceToPrey, y = CaptureSpeed, z = Density) ,linetype=2 ) +
    scale_x_continuous(name="Distance to prey (mm)", limits=c(0, 0.8)) +
    scale_y_continuous(name="Capture Speed (mm/sec)", limits=c(0, 80)) 
  #theme_linedraw()
  
  ggMarginal(p_NF, x="DistanceToPrey",y="CaptureSpeed", type = "density",groupColour = TRUE,groupFill=TRUE,show.legend=FALSE) 
}  

################################


## Analysis of Hunt Success 
### Makes Data frame With Number of Success Vs Failures From Labelled DatHuntEvent
getHuntSuccessPerFish <- function(datHuntLabelledEvents)
{
  #datHuntLabelledEvents <- datHuntLabelledEvents[datHuntLabelledEvents$eventID != 0,] ##Exclude the Artificial Event 0 Used such that all ExpID are In the DatHuntEvent
  tblResSB <- table(convertToScoreLabel(datHuntLabelledEvents$huntScore),datHuntLabelledEvents$groupID)
  tblFishScores <- table(datHuntLabelledEvents$expID, convertToScoreLabel(datHuntLabelledEvents$huntScore) )
  
  tblFishScoresLabelled<- tblFishScores[tblFishScores[,1] < 2, ] ##Pick Only THose ExpId (Fish) Whose Labelling Has (almost!) Finished
  ##Choose The Columns With the Scores Of Interest Success 3, Success-SpitBackOut 12 etc
  ##No_Targer is Column 5
  
  
  ## Find Tbl Indexes Indicating Success 
  tblIdxSuccess <- which (grepl("Success",row.names(tblResSB) ) ) 
  tblIdxFail <- which (grepl("Fail",row.names(tblResSB) ) ) 
  
  tblIdxNotHuntMode <- which (grepl("Out_Of_Range",row.names(tblResSB) ) | 
                                grepl("UnLabelled",row.names(tblResSB) ) | 
                                grepl("Duplicate",row.names(tblResSB) ) | 
                                grepl("Near-Hunt State",row.names(tblResSB) ) )  
  
  tblIdxEscape <- which (grepl("Escape",row.names(tblResSB) ) ) 
  tblIdxFail <- which (grepl("Fail",row.names(tblResSB) ) ) 
  
  
  
  datFishSuccessRate <- data.frame( cbind("Success" = rowSums(tblFishScoresLabelled[,tblIdxSuccess]),#tblFishScoresLabelled[,"Success"]+tblFishScoresLabelled[,"Success-SpitBackOut"]+tblFishScoresLabelled[,"Success-OnStrike"]+tblFishScoresLabelled[,"Success-OnStrike-SpitBackOut"]+tblFishScoresLabelled[,"Success-OnApproach"] +tblFishScoresLabelled[,"Success-OnApproach-AfterStrike"],
                                          "Fails_NS"= tblFishScoresLabelled[,"Fail-No Strike"],
                                          "Fails_WS"=tblFishScoresLabelled[,"Fail-With Strike"],
                                          "Fails"= rowSums(tblFishScoresLabelled[,tblIdxFail]),  #tblFishScoresLabelled[,"Fail"]+tblFishScoresLabelled[,"Fail-No Strike"]+tblFishScoresLabelled[,"Fail-With Strike"],
                                          "HuntEvents"=rowSums(tblFishScoresLabelled[, !(1:NCOL(tblFishScoresLabelled) %in% tblIdxNotHuntMode) ] ),  #rowSums(tblFishScoresLabelled[,c("Success","Success-SpitBackOut","Success-OnStrike","Success-OnStrike-SpitBackOut","Success-OnApproach","Success-OnApproach-AfterStrike","Fail","Fail-No Strike","Fail-With Strike","No_Target")]) , ##Ad The No Target To indicate Triggering Of Hunt Mode (Col 5)
                                          "CaptureEvents"=rowSums(tblFishScoresLabelled[,c(tblIdxSuccess,tblIdxFail)] ),  #rowSums(tblFishScoresLabelled[,c("Success","Success-SpitBackOut","Success-OnStrike","Success-OnStrike-SpitBackOut","Success-OnApproach","Success-OnApproach-AfterStrike","Fail","Fail-No Strike","Fail-With Strike","No_Target")]) , ##Ad The No Target To indicate Triggering Of Hunt Mode (Col 5)
                                          "expID"=NA,
                                          "groupID"=NA,
                                          "dataSetID"=NA) ) #
  
  vScoreIdx        <- ((datFishSuccessRate[,"Success"]*datFishSuccessRate[,"Success"])/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails"]))
  vEfficiency      <- ((datFishSuccessRate[,"Success"])/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails"]))  
  datFishSuccessRate <- cbind(datFishSuccessRate,HuntPower=vScoreIdx,Efficiency=vEfficiency)
  ##Add Group Label To the resulting Data Frame
  for (e in row.names(tblFishScoresLabelled) )
  {
    datFishSuccessRate[e,"expID"] <- unique( datHuntLabelledEvents[!is.na(datHuntLabelledEvents$expID) & datHuntLabelledEvents$expID == e,"expID"] )
    datFishSuccessRate[e,"groupID"] <- unique( datHuntLabelledEvents[!is.na(datHuntLabelledEvents$expID) & datHuntLabelledEvents$expID == e,"groupID"] )
    datFishSuccessRate[e,"dataSetID"] <- unique( datHuntLabelledEvents[!is.na(datHuntLabelledEvents$expID) & datHuntLabelledEvents$expID == e,"dataSetID"] )
  }
  
  return (datFishSuccessRate)
  
}


### Hunt outcome Labelling/Scorring convertion Functions
vHuntEventLabels <- c("UnLabelled","NA","Success","Fail","No_Target","Not_HuntMode/Delete","Escape","Out_Of_Range","Duplicate/Overlapping","Fail-No Strike","Fail-With Strike",
                      "Success-SpitBackOut",
                      "Debri-Triggered","Near-Hunt State","Success-OnStrike","Success-OnStrike-SpitBackOut",
                      "Success-OnApproach", ##For Ones that either do not strike but simply push on to the prey, 
                      "Success-OnApproach-AfterStrike" #or those that strike but only capture prey after and not during the strike (Denotes lesser ability to judge distance)
)

convertToScoreLabel <- function (huntScore) { 
  return (factor(x=huntScore,levels=seq(0,NROW(vHuntEventLabels)-1),labels=vHuntEventLabels ) )
}

####

## The figure with CDF of hunt efficiency - plotted as power, so as to show differences in actual consumption 
## instead of just ratio of success vs fail
plotHuntPowerDataCDF <- function(datHuntEventAllGroupToLabel)
{
  datFishSuccessRate <- getHuntSuccessPerFish(datHuntEventAllGroupToLabel)
  ## Hunt POwer ##
  vScoreIdx        <- datFishSuccessRate[,"HuntPower"] # ((datFishSuccessRate[,"Success"]*datFishSuccessRate[,"Success"])/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails"]))
  vEfficiencyRatio <- (datFishSuccessRate[,"Success"]/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails"]))
  vEfficiencyRatio_Strike <- (datFishSuccessRate[,"Success"]/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails_WS"]))
  vEfficiencyRatio_NStrike <- (datFishSuccessRate[,"Success"]/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails_NS"]))
  #vScoreIdx[is.nan(vScoreIdx) ] <- 0
  tblEventsTracked <- table(datHuntEventAllGroupToLabel$expID, datHuntEventAllGroupToLabel$markTracked,useNA="always" )
  datFishSuccessRateMerged <- cbind(datFishSuccessRate,
                                    markUnTrackable=data.frame(tblEventsTracked[row.names(datFishSuccessRate),1]),
                                    markTracked=data.frame(tblEventsTracked[row.names(datFishSuccessRate),2]),
                                    notTracked=data.frame(tblEventsTracked[row.names(datFishSuccessRate),3]))
  
  datFishSuccessRateMerged <- cbind(datFishSuccessRateMerged,vScoreIdx,vEfficiencyRatio,vEfficiencyRatio_Strike,vEfficiencyRatio_NStrike)
  
  
  
  ## Subset only the active/Larvae - ones that have hunted 
  datFishSuccessRateActive <- datFishSuccessRateMerged[!is.nan(datFishSuccessRateMerged$vScoreIdx),]
  
  ## Plot Density of Hunting POWER S^2/(S+F)
  densDLScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$vScoreIdx)
  densNLScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$vScoreIdx)
  densLLScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$vScoreIdx)
  
  cdfDLScore <- ecdf(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$vScoreIdx)
  cdfNLScore <- ecdf(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$vScoreIdx)
  cdfLLScore <- ecdf(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$vScoreIdx)
  
  
  #par(mar = c(3.9,4.3,1,1)) 
  plot(cdfNLScore,lty=2,lwd=3,col=colourLegL [1],xlim=c(0,12),pch=pchL[1],ylim=c(0.03,1.01),
       main=NA,ylab=NA,  xlab=NA,cex.main =cex,cex.axis=cex,cex=cex )
  
  plot(cdfLLScore,add=T,lty=1,lwd=3,col=colourLegL[2],pch=pchL[2],ylim=c(0,1.01),cex=cex)
  plot(cdfDLScore,add=T,lty=1,lwd=3,col=colourLegL[3],pch=pchL[3],ylim=c(0,1.01),cex=cex)
  #axis(side=1)
  mtext(side = 1,cex=cex, line = line, expression( "Hunt power " ~ N[S]^2/(N[S]+N[F]) ,paste("") )   )
  mtext(side = 2,cex=cex, line = line, expression("Cumulative function " ))
  
  legend("bottomright",legend=paste(c("NL #","LL #","DL #"),c(densNLScore$n,densLLScore$n,densDLScore$n) ),
         col = colourLegL,pch=pchL,cex=cex+0.2)
}




## The figure with CDF of hunt efficiency - plotted as power, so as to show differences in actual consumption 
## instead of just ratio of success vs fail
plotHuntEfficiencyDataCDF <- function(datHuntEventAllGroupToLabel,showLegend=FALSE)
{
  datFishSuccessRate <- getHuntSuccessPerFish(datHuntEventAllGroupToLabel)
  vScoreIdx        <- ((datFishSuccessRate[,"Success"]*datFishSuccessRate[,"Success"])/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails"]))
  vEfficiencyRatio <- (datFishSuccessRate[,"Success"]/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails"]))
  vEfficiencyRatio_Strike <- (datFishSuccessRate[,"Success"]/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails_WS"]))
  vEfficiencyRatio_NStrike <- (datFishSuccessRate[,"Success"]/(datFishSuccessRate[,"Success"]+datFishSuccessRate[,"Fails_NS"]))
  #vScoreIdx[is.nan(vScoreIdx) ] <- 0
  tblEventsTracked <- table(datHuntEventAllGroupToLabel$expID, datHuntEventAllGroupToLabel$markTracked,useNA="always" )
  datFishSuccessRateMerged <- cbind(datFishSuccessRate,
                                    markUnTrackable=data.frame(tblEventsTracked[row.names(datFishSuccessRate),1]),
                                    markTracked=data.frame(tblEventsTracked[row.names(datFishSuccessRate),2]),
                                    notTracked=data.frame(tblEventsTracked[row.names(datFishSuccessRate),3]))
  
  datFishSuccessRateMerged <- cbind(datFishSuccessRateMerged,vScoreIdx,vEfficiencyRatio,vEfficiencyRatio_Strike,vEfficiencyRatio_NStrike)
  
  
  
  ## Subset only the active/Larvae - ones that have hunted 
  datFishSuccessRateActive <- datFishSuccessRateMerged[!is.nan(datFishSuccessRateMerged$vEfficiencyRatio),]
  
  ## Plot Density of Hunting POWER S^2/(S+F)
  densDLScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$vEfficiencyRatio)
  densNLScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$vEfficiencyRatio)
  densLLScore <- density(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$vEfficiencyRatio)
  
  cdfDLScore <- ecdf(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "DL",]$vEfficiencyRatio)
  cdfNLScore <- ecdf(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "NL",]$vEfficiencyRatio)
  cdfLLScore <- ecdf(datFishSuccessRateActive[datFishSuccessRateActive$groupID == "LL",]$vEfficiencyRatio)
  
  
  #par(mar = c(3.9,4.3,1,1)) 
  plot(cdfNLScore,lty=2,lwd=3,col=colourLegL [1],xlim=c(0,1),pch=pchL[1],ylim=c(0.03,1.01),
       main=NA,ylab=NA,  xlab=NA,cex.main =cex,cex.axis=cex,cex=cex )
  
  plot(cdfLLScore,add=T,lty=1,lwd=3,col=colourLegL[2],pch=pchL[2],ylim=c(0,1.01),cex=cex)
  plot(cdfDLScore,add=T,lty=1,lwd=3,col=colourLegL[3],pch=pchL[3],ylim=c(0,1.01),cex=cex)
  #axis(side=1)
  mtext(side = 1,cex=cex, line = line, expression( "Capture Efficiency  "  ,paste("") )   )
  mtext(side = 2,cex=cex, line = line, expression("Cumulative function " ))
  if (showLegend) 
  {
    legend("bottomright",legend=paste(c("NL #","LL #","DL #"),c(densNLScore$n,densLLScore$n,densDLScore$n) ),
           col = colourLegL,pch=pchL,cex=cex+0.2)
  }
}


