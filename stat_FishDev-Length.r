## kostasl 11 Jun 2019
## stat Model for comparing larval development across the 3 rearing groups
## we assube 3 gaussians, one for each group, and  infer mean and std dev for each

library(rjags)
source("common_lib.R")


## Simple Gaussian model to estimate the mean larval length of each rearing group
modelNorm="model {
         
         for(i in 1:3) {
             mu[i]   ~ dunif(3,6 ) ## Mean length in mm
             prec[i] ~ dunif(0,50) ##dgamma(1,1) ##
             sigma[i] <- sqrt(1/prec[i]) 
        }

         for(i in 1:NTOT){
             LarvaLength[i] ~  dnorm(mu[groupID[i] ],prec[groupID[i]] )

         }
}"


################### LOAD CSV OF PX Larva lengths ################
message(paste(" Loading Measured fish length in pixels data ... "))
datFlatPxLength <- read.csv(file= paste0("dat/FishSTDFengths_ALFGroups.csv"))


strGroupID = levels(datFlatPxLength$groupID)

Jagsdata=list(groupID=datFlatPxLength$groupID,
              LarvaLength=datFlatPxLength$LengthPx*DIM_MMPERPX, 
              NTOT=nrow(datFlatPxLength));

varnames1=c("mu","sigma")
burn_in=1000;
steps=15000;
thin=10;
chains=3
ntail = 100 ##Number of Samples To Plot

strModelName = "modelFishSize.tmp"
fileConn=file(strModelName)
writeLines(modelNorm,fileConn);
close(fileConn)

sizemodel=jags.model(file=strModelName,data=Jagsdata,n.chains=chains);
update(sizemodel,burn_in)
draw=jags.samples(sizemodel,steps,thin=thin,variable.names=varnames1)

##
dmodelSizeNF <- density(draw$mu[which(strGroupID == "NF"),,],bw=0.1)
dmodelSizeLF <- density(draw$mu[which(strGroupID == "LF"),,],bw=0.1)
dmodelSizeDF <- density(draw$mu[which(strGroupID == "DF"),,],bw=0.1)

dSizeNF <- density(datFlatPxLength[datFlatPxLength$groupID == "NF",]$LengthPx*DIM_MMPERPX,bw=0.1)
dSizeLF <- density(datFlatPxLength[datFlatPxLength$groupID == "LF",]$LengthPx*DIM_MMPERPX,bw=0.1)
dSizeDF <- density(datFlatPxLength[datFlatPxLength$groupID == "DF",]$LengthPx *DIM_MMPERPX,bw=0.1)



xquant <- seq(0,6,diff(dSizeNF$x)[1])
XLIM <- c(3.5,5)
YLIM <- c(0,5)
pdistBW <- 2 ## mm/sec
strKern <- "gaussian"
ntail <- 20 #nrow(draw$mu[which(strGroupID == "NF"),,])*0.5
norm <- max(dSizeNF$y)

##T- Test Between NF and LF Larvae Shows Significance 
t.test(datFlatPxLength[datFlatPxLength$groupID == "NF",]$LengthPx*DIM_MMPERPX,
       datFlatPxLength[datFlatPxLength$groupID == "LF",]$LengthPx*DIM_MMPERPX)

##T- Test Between DF and LF Larvae Shows Significance 
t.test(datFlatPxLength[datFlatPxLength$groupID == "DF",]$LengthPx*DIM_MMPERPX,
       datFlatPxLength[datFlatPxLength$groupID == "LF",]$LengthPx*DIM_MMPERPX)


#######PLOT RESULTS
## FIGURE CONTROL FOR LARVAL SIZE
#### PLOT Fitted Larva Length  ###
pdf(file= paste(strPlotExportPath,"/stat/stat_LarvaLFengthsGaussian.pdf" ,sep=""),width = 14,height = 3.5)
  
  layout(matrix(c(1,2,3,4),1,4, byrow = TRUE))
  ##Margin: (Bottom,Left,Top,Right )
  par(mar = c(3.9,4.9,1,1))
 
  line = 8.3 ## SubFig Label Params
  lineAxis = 3.0
  cex = 1.4
  adj  = -3
  padj <- -7.2
  las <- 1
  
  plot(dSizeNF,col=colourLegL[1] , xlim=XLIM,ylim=YLIM ,lwd=2 ,type='l',main=NA,xlab=NA,ylab=NA,lty=2,
       cex=cex,cex.axis=cex,cex.lab=cex)
  ##Show Model Fit
  for (i in 1:(ntail-1) )
  {
    lines(xquant,
         dnorm(xquant,mean=tail(draw$mu[which(strGroupID == "NF"),ntail-i,],1),
                                  sd=tail(draw$sigma[which(strGroupID == "NF"),ntail-i,],1))
         ,type='l', col=colourLegL[1],lwd=1,lty=1)
  }
  lines(dSizeNF,col="black",lwd=4,lty=2)
  mtext(side = 1,cex=cex, line = lineAxis, expression("Larval length (mm) " ) )
  mtext(side = 2,cex=cex, line = lineAxis, expression("Density function " ))
  mtext("A",at="topleft",outer=F,side=2,col="black",las=las,font=2,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex) #las=las,line=line,padj=padj,adj=adj
  legend("topright",title="NF",
         legend=c( paste0("",dSizeNF$n, "# Data "), #(Bw:",prettyNum(digits=2, pdistBW ),")" ) ,
                   paste("Model  " )),
         col=c("black",colourLegL[1]),lwd=c(3,1),lty=c(2,1),seg.len = 3,cex=cex ) 
  
  ## 2nd panel
  plot(dSizeLF,col=colourLegL[2] , xlim=XLIM,ylim=YLIM ,lwd=2 ,type='l',main=NA,xlab=NA,ylab=NA,lty=2,cex=cex,cex.axis=cex,cex.lab=cex )
  ##Show Model Fit
  for (i in 1:(ntail-1) )
  {
    lines(xquant,
          dnorm(xquant,mean=tail(draw$mu[which(strGroupID == "LF"),ntail-i,],1),
                sd=tail(draw$sigma[which(strGroupID == "LF"),ntail-i,],1)),
          type='l', col=colourLegL[2],lwd=1,lty=1)
  }
  lines(dSizeLF,col="black" , xlim=XLIM,ylim=YLIM ,lwd=4 ,type='l',xlab=NA,ylab=NA,lty=2 )
  mtext(side = 1,cex=cex, line = lineAxis, expression("Larval length (mm) " ))
  mtext(side = 2,cex=cex, line = lineAxis, expression("Density function " ))
  mtext("B",at="topleft",outer=F,side=2,col="black",las=las,font=2,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex) #las=las,line=line,padj=padj,adj=adj
  legend("topright",title="LF",
         legend=c( paste0("",dSizeLF$n, "# Data "), #(Bw:",prettyNum(digits=2, pdistBW ),")" ) ,
                   paste("Model  " )),
         col=c("black",colourLegL[2]),lwd=c(3,1),lty=c(2,1),seg.len = 3,cex=cex) 
  
  
  plot(dSizeDF,col=colourLegL[3] , xlim=XLIM,ylim=YLIM ,lwd=2 ,type='l',main=NA,xlab=NA,ylab=NA,lty=2,cex=cex,cex.axis=cex,cex.lab=cex )
  ##Show Model Fit
  for (i in 1:(ntail-1) )
  {
    lines(xquant,
          dnorm(xquant,mean=tail(draw$mu[which(strGroupID == "DF"),ntail-i,],1),
                sd=tail(draw$sigma[which(strGroupID == "DF"),ntail-i,],1))
          ,type='l',col=colourLegL[3],lwd=1,lty=1)
  }
  lines(dSizeDF,col="black",xlim=XLIM,lwd=4,lty=2)
  
  mtext(side = 1,cex=cex, line = lineAxis, expression("Larval length (mm) " ))
  mtext(side = 2,cex=cex, line = lineAxis, expression("Density function " ))
  mtext("C",at="topleft",outer=F,side=2,col="black",las=las,font=2,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex) #las=las,line=line,padj=padj,adj=adj
  legend("topright",title="DF",
         legend=c( paste0("",dSizeDF$n, "# Data "), #(Bw:",prettyNum(digits=2, pdistBW ),")" ) ,
                   paste("Model  " )),
         col=c("black",colourLegL[3]),lwd=c(3,1),lty=c(2,1),seg.len = 3 ,cex=cex) 
  
  ##MPlot parameter Density
  plot(dmodelSizeNF,col=colourLegL[1] , xlim=XLIM,ylim=YLIM,lwd=4 ,type='l',main=NA,xlab=NA,ylab=NA,lty=1,cex=cex,cex.axis=cex,cex.lab=cex )
  lines(dmodelSizeLF,col=colourLegL[2] , xlim=XLIM ,lwd=4 ,type='l',main=NA,xlab=NA,ylab=NA,lty=2 )
  lines(dmodelSizeDF,col=colourLegL[3] , xlim=XLIM ,lwd=4 ,type='l',main=NA,xlab=NA,ylab=NA,lty=3 )
  legend("topright", legend=c("NF","LF","DF"), 
         col=colourLegL,lty=c(1,2,3),lwd=3,seg.len=4,cex=cex)
  mtext(side = 1,cex=cex, line = lineAxis, expression("Estimated mean length (mm) " ))
  mtext(side = 2,cex=cex, line = lineAxis, expression("Density function " ))
  mtext("D",at="topleft",outer=F,side=2,col="black",las=las,font=2,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex) #las=las,line=line,padj=padj,adj=adj

dev.off()

