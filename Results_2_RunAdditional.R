#------------------------------------------------------------------------------#
# Results_1_Run.R                                                              #
# Analyse and present intervention results by round                            #
# Last updated 2025-07-25 by KCH                                               #
#------------------------------------------------------------------------------#

# Clear environment
rm(list = ls())

# Load required libraries
library(deSolve)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(stringr)
library(patchwork)
library(ggbreak)
library(gt)
library(gtsummary)
library(ggpubr)
library(webshot2)
library(writexl)

# Set working directory (to this file's location)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# Assign intervention colours
colnoint     <- '#5c5c5c'   #'#283845' #'#899E8B'
colsympxpert <- '#BC3C29FF' #'#8E5572' #'#1F4E79' '#4A4E69'
colcxrxpert  <- '#4A4E69'   #'#4A4E69' #'#8E5572' '#C08497'
colxpert     <- '#E18727FF' #'#537D8D' #'#899E8B' '#899E8B'
colcxr       <- '#6F99ADFF' #'#899E8B' #'#2C6E49' '#88665D' '#537D8D'


# Assign disease state colours
colsTB <- colorRampPalette(c("white",colxpert))(5)[[5]]  #'#8C2F39' 
colaTB <- colorRampPalette(c("white",colxpert))(5)[[3]]  #'#F29559'
colnTB <- colorRampPalette(c("white",colxpert))(5)[[2]]  #'#FFD166'
colnot <- colorRampPalette(c("white",colxpert))(5)[[5]]  #'#8AA29E'

# Set scenario parameters
pop      <- 100000  # Population
plotmos  <- 332     # Number of months of follow-up
datapts  <- 51      # Number of data points in model output
outnum   <- 15      # Number of outcomes 
algnum   <- 5       # Number of screening algorithms
covnum   <- 10      # Number of screening coverage levels
durnum   <- 5       # Number of screening years
prevlist <- c(250,500,1000)                             # List of prevalences
covlist  <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)  # List of coverages
durlist  <- c(5)                                        # List of durations

# Source functions and files
source("./Results_0_Functions.R")   

# Calculate outputs
for (prev in prevlist){

  # Create folder to store results
  dir.create(file.path(paste("./Results_",prev,sep="")), showWarnings=FALSE)

  # Set up output structure
  conpostb1    <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  conpostb2    <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  conpostb3    <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  conpostb4    <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  conpostb5    <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  conposnot1   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  conposnot2   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  conposnot3   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  conposnot4   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  conposnot5   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  fptp1   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  fptp2   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  fptp3   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  fptp4   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  fptp5   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  
  yr5results <- array(NA,dim=c(covnum,durnum,outnum,(algnum-1),3))
  
  # Calculate outputs of interest
  for (cov in covlist){
    for (int in durlist){
     
      if (cov==0.1) covn <-  1
      if (cov==0.2) covn <-  2
      if (cov==0.3) covn <-  3
      if (cov==0.4) covn <-  4
      if (cov==0.5) covn <-  5
      if (cov==0.6) covn <-  6
      if (cov==0.7) covn <-  7
      if (cov==0.8) covn <-  8
      if (cov==0.9) covn <-  9
      if (cov==1.0) covn <- 10
      
      # Load model output
      load(file=paste("./Output_",prev,"/run_sets_",int,"_",cov,".Rdata",sep=""))
      fulldata <- array(NA,dim=c(length(run_sets),algnum,plotmos,datapts))
      
      # Restructure data
      for (i in 1:length(run_sets)){
        print(i)
        run_sets_i <- run_sets[[i]]
        
        load(file=paste("./Output_",prev,"_Xpert","/run_sets_",int,"_",cov,".Rdata",sep=""))
        run_sets_i_Xpert <- run_sets[[i]][[1]]
        
        set <- list(as.data.frame(run_sets_i[[1]]),
                    as.data.frame(run_sets_i[[2]]), 
                    as.data.frame(run_sets_i[[3]]),
                    as.data.frame(run_sets_i_Xpert),
                    as.data.frame(run_sets_i[[5]]))
        
        load(file=paste("./Output_",prev,"/run_sets_",int,"_",cov,".Rdata",sep=""))
        
        for(n in 1:length(set)){
          fulldata[i,n,,] <- data.matrix(set[[n]][1:plotmos,1:datapts])
          
          # Number of individuals with positive confirmatory results
          conpostb1[i,n]    <- sum(set[[n]][24,c("ConPossTB",
                                                    "ConPosaTB",
                                                    "ConPosnTB")])
          conposnot1[i,n]   <- sum(set[[n]][24,c("ConPosS",
                                                    "ConPosI",
                                                    "ConPosRn",
                                                    "ConPosRt")])
          fptp1[i,n]        <- conposnot1[i,n] / conpostb1[i,n] 
          conpostb2[i,n]    <- sum(set[[n]][36,c("ConPossTB",
                                                    "ConPosaTB",
                                                    "ConPosnTB")]) - 
                               sum(set[[n]][24,c("ConPossTB",
                                                    "ConPosaTB",
                                                    "ConPosnTB")])
          conposnot2[i,n]   <- sum(set[[n]][36,c("ConPosS",
                                                    "ConPosI",
                                                    "ConPosRn",
                                                    "ConPosRt")]) - 
                               sum(set[[n]][24,c("ConPosS",
                                                    "ConPosI",
                                                    "ConPosRn",
                                                    "ConPosRt")])
          fptp2[i,n]        <- conposnot2[i,n] / conpostb2[i,n]
          conpostb3[i,n]    <- sum(set[[n]][48,c("ConPossTB",
                                                    "ConPosaTB",
                                                    "ConPosnTB")]) - 
                               sum(set[[n]][36,c("ConPossTB",
                                                    "ConPosaTB",
                                                    "ConPosnTB")])
          conposnot3[i,n]   <- sum(set[[n]][48,c("ConPosS",
                                                    "ConPosI",
                                                    "ConPosRn",
                                                    "ConPosRt")]) - 
                               sum(set[[n]][36,c("ConPosS",
                                                      "ConPosI",
                                                      "ConPosRn",
                                                      "ConPosRt")])
          fptp3[i,n]        <- conposnot3[i,n] / conpostb3[i,n] 
          conpostb4[i,n]    <- sum(set[[n]][60,c("ConPossTB",
                                                    "ConPosaTB",
                                                    "ConPosnTB")]) - 
                               sum(set[[n]][48,c("ConPossTB",
                                                      "ConPosaTB",
                                                      "ConPosnTB")])
          conposnot4[i,n]   <- sum(set[[n]][60,c("ConPosS",
                                                    "ConPosI",
                                                    "ConPosRn",
                                                    "ConPosRt")]) - 
                               sum(set[[n]][48,c("ConPosS",
                                                      "ConPosI",
                                                      "ConPosRn",
                                                      "ConPosRt")])
          fptp4[i,n]        <- conposnot4[i,n] / conpostb4[i,n]
          conpostb5[i,n]    <- sum(set[[n]][72,c("ConPossTB",
                                                    "ConPosaTB",
                                                    "ConPosnTB")]) - 
            sum(set[[n]][60,c("ConPossTB",
                                                      "ConPosaTB",
                                                      "ConPosnTB")])
          conposnot5[i,n]   <- sum(set[[n]][72,c("ConPosS",
                                                    "ConPosI",
                                                    "ConPosRn",
                                                    "ConPosRt")]) - 
            sum(set[[n]][60,c("ConPosS",
                                                      "ConPosI",
                                                      "ConPosRn",
                                                      "ConPosRt")])
          fptp5[i,n]        <- conposnot5[i,n] / conpostb5[i,n]
          
        }
      }
      
      ## SCREENING PERFORMANCE ##
      
      # Number of individuals with positive results: TB
      yr5results[covn,int, 1,,] <- rbind(quantile(conpostb1$sympxpert,probs=c(0.025,0.5,0.975)),
                                         quantile(conpostb1$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                         quantile(conpostb1$xpert    ,probs=c(0.025,0.5,0.975)),
                                         quantile(conpostb1$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of individuals with positive results: notTB
      yr5results[covn,int, 2,,] <- rbind(quantile(conposnot1$sympxpert,probs=c(0.025,0.5,0.975)),
                                         quantile(conposnot1$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                         quantile(conposnot1$xpert    ,probs=c(0.025,0.5,0.975)),
                                         quantile(conposnot1$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of individuals with positive results: TB
      yr5results[covn,int, 3,,] <- rbind(quantile(conpostb2$sympxpert,probs=c(0.025,0.5,0.975)),
                                         quantile(conpostb2$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                         quantile(conpostb2$xpert    ,probs=c(0.025,0.5,0.975)),
                                         quantile(conpostb2$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of individuals with positive results: notTB
      yr5results[covn,int, 4,,] <- rbind(quantile(conposnot2$sympxpert,probs=c(0.025,0.5,0.975)),
                                         quantile(conposnot2$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                         quantile(conposnot2$xpert    ,probs=c(0.025,0.5,0.975)),
                                         quantile(conposnot2$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of individuals with positive results: TB
      yr5results[covn,int, 5,,] <- rbind(quantile(conpostb3$sympxpert,probs=c(0.025,0.5,0.975)),
                                         quantile(conpostb3$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                         quantile(conpostb3$xpert    ,probs=c(0.025,0.5,0.975)),
                                         quantile(conpostb3$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of individuals with positive results: notTB
      yr5results[covn,int, 6,,] <- rbind(quantile(conposnot3$sympxpert,probs=c(0.025,0.5,0.975)),
                                         quantile(conposnot3$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                         quantile(conposnot3$xpert    ,probs=c(0.025,0.5,0.975)),
                                         quantile(conposnot3$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of individuals with positive results: TB
      yr5results[covn,int, 7,,] <- rbind(quantile(conpostb4$sympxpert,probs=c(0.025,0.5,0.975)),
                                         quantile(conpostb4$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                         quantile(conpostb4$xpert    ,probs=c(0.025,0.5,0.975)),
                                         quantile(conpostb4$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of individuals with positive results: notTB
      yr5results[covn,int, 8,,] <- rbind(quantile(conposnot4$sympxpert,probs=c(0.025,0.5,0.975)),
                                         quantile(conposnot4$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                         quantile(conposnot4$xpert    ,probs=c(0.025,0.5,0.975)),
                                         quantile(conposnot4$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of individuals with positive results: TB
      yr5results[covn,int, 9,,] <- rbind(quantile(conpostb5$sympxpert,probs=c(0.025,0.5,0.975)),
                                         quantile(conpostb5$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                         quantile(conpostb5$xpert    ,probs=c(0.025,0.5,0.975)),
                                         quantile(conpostb5$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of individuals with positive results: notTB
      yr5results[covn,int,10,,] <- rbind(quantile(conposnot5$sympxpert,probs=c(0.025,0.5,0.975)),
                                         quantile(conposnot5$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                         quantile(conposnot5$xpert    ,probs=c(0.025,0.5,0.975)),
                                         quantile(conposnot5$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of false positives per true positive
      yr5results[covn,int,11,,] <- rbind(quantile(fptp1$sympxpert,probs=c(0.025,0.5,0.975)),
                                         quantile(fptp1$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                         quantile(fptp1$xpert    ,probs=c(0.025,0.5,0.975)),
                                         quantile(fptp1$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of false positives per true positive
      yr5results[covn,int,12,,] <- rbind(quantile(fptp2$sympxpert,probs=c(0.025,0.5,0.975)),
                                         quantile(fptp2$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                         quantile(fptp2$xpert    ,probs=c(0.025,0.5,0.975)),
                                         quantile(fptp2$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of false positives per true positive
      yr5results[covn,int,13,,] <- rbind(quantile(fptp3$sympxpert,probs=c(0.025,0.5,0.975)),
                                         quantile(fptp3$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                         quantile(fptp3$xpert    ,probs=c(0.025,0.5,0.975)),
                                         quantile(fptp3$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of false positives per true positive
      yr5results[covn,int,14,,] <- rbind(quantile(fptp4$sympxpert,probs=c(0.025,0.5,0.975)),
                                         quantile(fptp4$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                         quantile(fptp4$xpert    ,probs=c(0.025,0.5,0.975)),
                                         quantile(fptp4$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of false positives per true positive
      yr5results[covn,int,15,,] <- rbind(quantile(fptp5$sympxpert,probs=c(0.025,0.5,0.975)),
                                         quantile(fptp5$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                         quantile(fptp5$xpert    ,probs=c(0.025,0.5,0.975)),
                                         quantile(fptp5$cxr      ,probs=c(0.025,0.5,0.975)))
      
      save(yr5results,file=paste("./Results_",prev,"/res/yr5results.Rdata",sep=""))
    }
  }  

  for (covn in c(1:covnum)){
    for (int in c(1:durnum)){
      for (out in c(1:outnum)){
        temp <- array(0,dim=c(4,3))
        temp <- yr5results[covn,int,out,,]
        write.csv(temp,file=paste("./Results_",prev,"/res/yr5res_",out,"_",covn,"_",int,".csv",sep=""))
      }
    }
  }
}

fptplist <- list()
ylim     <- list()
fptpdata <- c()
tpdata <- c()
fpdata <- c()
for (prev in prevlist){
  # True positives by round
  load(paste("./Results_",prev,"/res/yr5results.Rdata",sep=""))
  
  plotlist <- list()
  ylimit <- max(yr5results[,5,c(1,3,5,7,9),,3])
  
  for (test in c(1:(algnum-1))){
    outdata <- c()
    if (test==1){
      testcol  <- colsympxpert
      testname <- "Cough+Xpert"
    }
    if (test==2){
      testcol  <- colcxrxpert
      testname <- "CXR+Xpert"
    }
    if (test==3){
      testcol  <- colxpert
      testname <- "Xpert"
    }
    if (test==4){
      testcol  <- colcxr
      testname <- "CXR"
    }
    testcolrange <- colorRampPalette(c("white",testcol))
    for (covn in c(1:covnum)){
      for (int in c(1:durnum)){
        
      if (covn== 1) cov <- '10%'
      if (covn== 2) cov <- '20%'
      if (covn== 3) cov <- '30%'
      if (covn== 4) cov <- '40%'
      if (covn== 5) cov <- '50%'
      if (covn== 6) cov <- '60%'
      if (covn== 7) cov <- '70%'
      if (covn== 8) cov <- '80%'
      if (covn== 9) cov <- '90%'
      if (covn==10) cov <- '100%'
    
      temp <- yr5results[covn,5,int*2-1,test,]
      tpdata <- rbind(tpdata,c(temp,testname,cov,int,prev))
      outdata <- rbind(outdata,c(temp,covn,cov,int))
      }
    }
    outdata <- as.data.frame(outdata)
    colnames(outdata) <- c("lower","median","upper","covn","cov","int")
    outdata$int <- as.character(outdata$int)  
    outdata$median <- as.numeric(outdata$median)  
    outdata$lower  <- as.numeric(outdata$lower)  
    outdata$upper  <- as.numeric(outdata$upper)  
    outdataplot <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
      geom_bar(position=position_dodge(),stat="identity",colour="black") +
      geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
      scale_x_discrete("Population coverage") +
      scale_y_continuous("True positives per round",expand=c(0,0),labels=scales::comma) + #,lim=c(0,ylimit)) +
      scale_fill_manual("Rounds", values=testcolrange(durnum+1)[2:(durnum+1)]) +
      ggtitle(testname) +
      theme_classic() +
      theme(legend.position="none",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5)#,axis.title.y=element_blank()
      ) + #
      guides(fill=guide_legend(title.position="top",title.hjust=0.5))
    plotlist[[test]] <- outdataplot
  }
  figtp <- (plotlist[[1]] / plotlist[[3]] / plotlist[[4]])
  ggsave(filename=paste("./Results_",prev,"/plot_tpbyround.tiff",sep=""),plot=figtp,width=20,height=20, unit="cm")

  # False positives by round
  load(paste("./Results_",prev,"/res/yr5results.Rdata",sep=""))
  
  plotlist <- list()
  ylimit <- max(yr5results[,5,c(2,4,6,8,10),,3])
  
  for (test in c(1:(algnum-1))){
    outdata <- c()
    if (test==1){
      testcol  <- colsympxpert
      testname <- "Cough+Xpert"
    }
    if (test==2){
      testcol  <- colcxrxpert
      testname <- "CXR+Xpert"
    }
    if (test==3){
      testcol  <- colxpert
      testname <- "Xpert"
    }
    if (test==4){
      testcol  <- colcxr
      testname <- "CXR"
    }
    testcolrange <- colorRampPalette(c("white",testcol))
    for (covn in c(1:covnum)){
      for (int in c(1:durnum)){
        
        if (covn== 1) cov <- '10%'
        if (covn== 2) cov <- '20%'
        if (covn== 3) cov <- '30%'
        if (covn== 4) cov <- '40%'
        if (covn== 5) cov <- '50%'
        if (covn== 6) cov <- '60%'
        if (covn== 7) cov <- '70%'
        if (covn== 8) cov <- '80%'
        if (covn== 9) cov <- '90%'
        if (covn==10) cov <- '100%'
        
        temp <- yr5results[covn,5,int*2,test,]
        fpdata <- rbind(fpdata,c(temp,testname,cov,int,prev))
        outdata <- rbind(outdata,c(temp,covn,cov,int))
      }
    }
    outdata <- as.data.frame(outdata)
    colnames(outdata) <- c("lower","median","upper","covn","cov","int")
    outdata$int <- as.character(outdata$int)  
    outdata$median <- as.numeric(outdata$median)  
    outdata$lower  <- as.numeric(outdata$lower)  
    outdata$upper  <- as.numeric(outdata$upper)  
    outdataplot <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
      geom_bar(position=position_dodge(),stat="identity",colour="black") +
      geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
      scale_x_discrete("Population coverage") +
      scale_y_continuous("False positives per round",expand=c(0,0),labels=scales::comma) + #,lim=c(0,ylimit)) +
      scale_fill_manual("Rounds", values=testcolrange(durnum+1)[2:(durnum+1)]) +
      ggtitle(testname) +
      theme_classic() +
      theme(legend.position="none",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5)#,axis.title.y=element_blank()
      ) + #
      guides(fill=guide_legend(title.position="top",title.hjust=0.5))
    plotlist[[test]] <- outdataplot
  }
  figfp <- (plotlist[[1]] / plotlist[[3]] / plotlist[[4]])
  ggsave(filename=paste("./Results_",prev,"/plot_fpbyround.tiff",sep=""),plot=figfp,width=20,height=20, unit="cm")

  #FP:TP by round
  load(paste("./Results_",prev,"/res/yr5results.Rdata",sep=""))
  
  plotlist <- list()
  for (test in c(1:(algnum-1))){
    outdata <- c()
    if (test==1){
      testcol  <- colsympxpert
      testname <- "Cough+Xpert"
    }
    if (test==2){
      testcol  <- colcxrxpert
      testname <- "CXR+Xpert"
    }
    if (test==3){
      testcol  <- colxpert
      testname <- "Xpert"
    }
    if (test==4){
      testcol  <- colcxr
      testname <- "CXR"
    }
    testcolrange <- colorRampPalette(c("white",testcol))
    if (prev == 250){
      if (test == 4)  ylim[[test]] <- max(yr5results[,5,c(11:15),test,2])*1.05
      else            ylim[[test]] <- max(yr5results[,5,c(11:15),test,3])*1.05
    }
    ylimit <- ylim[[test]]
    for (covn in c(1:covnum)){
      for (int in c(1:durnum)){
        
        if (covn== 1) cov <- '10%'
        if (covn== 2) cov <- '20%'
        if (covn== 3) cov <- '30%'
        if (covn== 4) cov <- '40%'
        if (covn== 5) cov <- '50%'
        if (covn== 6) cov <- '60%'
        if (covn== 7) cov <- '70%'
        if (covn== 8) cov <- '80%'
        if (covn== 9) cov <- '90%'
        if (covn==10) cov <- '100%'
        
        temp <- yr5results[covn,5,int+10,test,]
        fptpdata <- rbind(fptpdata,c(temp,testname,cov,int,prev))
        if (prev == 250 && test == 4 && cov == "100%" && int == 5){
          temp[3] <- temp[2]
        }
        outdata <- rbind(outdata,c(temp,covn,cov,int))
      }
    }
    outdata <- as.data.frame(outdata)
    colnames(outdata) <- c("lower","median","upper","covn","cov","int")
    outdata$int <- as.character(outdata$int)  
    outdata$median <- as.numeric(outdata$median)  
    outdata$lower  <- as.numeric(outdata$lower)  
    outdata$upper  <- as.numeric(outdata$upper)  
    arrow <- outdata %>%
      filter(int == 5) %>%
      filter(cov == "100%")
      
    if (prev == 250 && test == 1){
      outdataplot <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
        geom_bar(position=position_dodge(),stat="identity",colour="black") +
        geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
        scale_x_discrete("") +
        scale_y_continuous("",expand=c(0,0),labels=scales::comma,lim=c(0,ylimit)) +
        scale_fill_manual("Round", values=testcolrange(durnum+1)[2:(durnum+1)]) +
        ggtitle(paste("Baseline prevalence of ",prev," per 100,000",sep="")) +
        theme_classic(base_size=18) +
        theme(legend.position="none",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank()
        ) + #
        guides(fill=guide_legend(title.position="top",title.hjust=0.5))
    }else if (prev == 250 && test == 4){
      outdataplot <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
        geom_bar(position=position_dodge(),stat="identity",colour="black") +
        geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
        geom_segment(arrow, mapping=aes(x='100%',xend='100%',y=lower,yend=median*1.04),arrow=arrow(),position=position_nudge(.36)) +
        scale_x_discrete("") +
        scale_y_continuous("",expand=c(0,0),labels=scales::comma,lim=c(0,ylimit)) +
        scale_fill_manual("Round", values=testcolrange(durnum+1)[2:(durnum+1)]) +
        ggtitle("") +
        theme_classic(base_size=18) +
        theme(legend.position="none",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5),axis.text.y=element_blank()
        ) + #
        guides(fill=guide_legend(title.position="top",title.hjust=0.5))
    }else if (prev == 250 && !test %in% c(1,4)){
      outdataplot <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
        geom_bar(position=position_dodge(),stat="identity",colour="black") +
        geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
        scale_x_discrete("") +
        scale_y_continuous("",expand=c(0,0),labels=scales::comma,lim=c(0,ylimit)) +
        scale_fill_manual("Round", values=testcolrange(durnum+1)[2:(durnum+1)]) +
        ggtitle("") +
        theme_classic(base_size=18) +
        theme(legend.position="none",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank()
        ) + #
        guides(fill=guide_legend(title.position="top",title.hjust=0.5))
    }else if (prev == 500 && test == 1){
      outdataplot <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
        geom_bar(position=position_dodge(),stat="identity",colour="black") +
        geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
        scale_x_discrete("") +
        scale_y_continuous("",expand=c(0,0),labels=scales::comma,lim=c(0,ylimit)) +
        scale_fill_manual("Round", values=testcolrange(durnum+1)[2:(durnum+1)]) +
        ggtitle(paste("Baseline prevalence of ",prev," per 100,000",sep="")) +
        theme_classic(base_size=18) +
        theme(legend.position="none",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank()
        ) + #
        guides(fill=guide_legend(title.position="top",title.hjust=0.5))
    }else if (prev == 500 && test == 4){
      outdataplot <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
        geom_bar(position=position_dodge(),stat="identity",colour="black") +
        geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
        scale_x_discrete("Population coverage") +
        scale_y_continuous("",expand=c(0,0),labels=scales::comma,lim=c(0,ylimit)) +
        scale_fill_manual("Round", values=testcolrange(durnum+1)[2:(durnum+1)]) +
        ggtitle("") +
        theme_classic(base_size=18) +
        theme(legend.position="none",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5),axis.title.y=element_blank(),axis.text.y=element_blank()
        ) + #
        guides(fill=guide_legend(title.position="top",title.hjust=0.5))
    }else if (prev == 1000 && test == 1){
      outdataplot <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
        geom_bar(position=position_dodge(),stat="identity",colour="black") +
        geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
        scale_x_discrete("") +
        scale_y_continuous(paste(" \n",testname,sep=""),expand=c(0,0),labels=scales::comma,lim=c(0,ylimit)) + 
        scale_fill_manual("Round", values=testcolrange(durnum+1)[2:(durnum+1)]) +
        ggtitle(paste("Baseline prevalence of ",prev," per 100,000",sep="")) +
        theme_classic(base_size=18) +
        theme(legend.position="none",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5),axis.title.x=element_blank(),axis.text.x=element_blank()
        ) + #
        guides(fill=guide_legend(title.position="top",title.hjust=0.5))
    }else if (prev == 1000 && test == 3){
      outdataplot <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
        geom_bar(position=position_dodge(),stat="identity",colour="black") +
        geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
        scale_x_discrete("") +
        scale_y_continuous(paste("False positives per true positive\n",testname,sep=""),expand=c(0,0),labels=scales::comma,lim=c(0,ylimit)) + 
        scale_fill_manual("Round", values=testcolrange(durnum+1)[2:(durnum+1)]) +
        ggtitle("") +
        theme_classic(base_size=18) +
        theme(legend.position="none",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5),axis.title.x=element_blank(),axis.text.x=element_blank()
        ) + #
        guides(fill=guide_legend(title.position="top",title.hjust=0.5))
    }else if (prev == 1000 && test == 4){
      outdataplot <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
        geom_bar(position=position_dodge(),stat="identity",colour="black") +
        geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
        scale_x_discrete("") +
        scale_y_continuous(paste(" \n",testname,sep=""),expand=c(0,0),labels=scales::comma,lim=c(0,ylimit)) +
        scale_fill_manual("Round", values=testcolrange(durnum+1)[2:(durnum+1)]) +
        ggtitle("") +
        theme_classic(base_size=18) +
        theme(legend.position="none",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5)
        ) + #
        guides(fill=guide_legend(title.position="top",title.hjust=0.5))
    }else{
      outdataplot <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
        geom_bar(position=position_dodge(),stat="identity",colour="black") +
        geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
        scale_x_discrete("") +
        scale_y_continuous("",expand=c(0,0),labels=scales::comma,lim=c(0,ylimit)) +
        scale_fill_manual("Round", values=testcolrange(durnum+1)[2:(durnum+1)]) +
        ggtitle("") +
        theme_classic(base_size=18) +
        theme(legend.position="none",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank()
        ) + #
        guides(fill=guide_legend(title.position="top",title.hjust=0.5))
    }
    plotlist[[test]] <- outdataplot
    if(prev == 250){
      fptplist[[test]] <- outdataplot
    }else if (prev==500){
      fptplist[[test+4]] <- outdataplot
    }else if (prev==1000){
      fptplist[[test+8]] <- outdataplot
    }
  }
  figfptp <- (plotlist[[1]] / plotlist[[3]] / plotlist[[4]])
  ggsave(filename=paste("./Results_",prev,"/plot_fptpbyround2.tiff",sep=""),plot=figfptp,width=20,height=20, unit="cm")
  
#  fig <- figtp | figfp | figfptp
#  ggsave(filename=paste("./Results_",prev,"/plot_resbyround.tiff",sep=""),plot=fig,width=60,height=20, unit="cm")
  
}

colnames(fptpdata) <- c("low","med","high","approach","cov","int","prev")
fptpdata <- as.data.frame(fptpdata)
tab <- fptpdata %>%
  mutate(combo = paste(round(as.numeric(med ),digits=1), " \n(", 
                       round(as.numeric(low ),digits=1), "-",
                       round(as.numeric(high),digits=1), ")", sep="")) %>%
  select(approach,prev,cov,int,combo)
for (prevlev in prevlist){
  fptp <- tab %>%
    filter(prev  == prevlev) %>%
    pivot_wider(names_from=cov,values_from=combo) %>%
    arrange(factor(approach,levels=c("Cough+Xpert","Xpert","CXR"))) %>%
    select(!c(prev))
  colnames(fptp) <- c("Algorithm","Rounds","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")
  write_xlsx(fptp,paste("./Results_",prevlev,"/tablefptp.xlsx",sep=""))
}

colnames(tpdata) <- c("low","med","high","approach","cov","int","prev")
tpdata <- as.data.frame(tpdata)
tab <- tpdata %>%
  mutate(combo = paste(round(as.numeric(med ),digits=0), " \n(", 
                       round(as.numeric(low ),digits=0), "-",
                       round(as.numeric(high),digits=0), ")", sep="")) %>%
  select(approach,prev,cov,int,combo)
for (prevlev in prevlist){
  tp <- tab %>%
    filter(prev  == prevlev) %>%
    pivot_wider(names_from=cov,values_from=combo) %>%
    arrange(factor(approach,levels=c("Cough+Xpert","Xpert","CXR"))) %>%
    select(!c(prev))
  colnames(tp) <- c("Algorithm","Rounds","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")
  write_xlsx(tp,paste("./Results_",prevlev,"/tabletp.xlsx",sep=""))
}

colnames(fpdata) <- c("low","med","high","approach","cov","int","prev")
fpdata <- as.data.frame(fpdata)
tab <- fpdata %>%
  mutate(combo = paste(round(as.numeric(med ),digits=0), " \n(", 
                       round(as.numeric(low ),digits=0), "-",
                       round(as.numeric(high),digits=0), ")", sep="")) %>%
  select(approach,prev,cov,int,combo)
for (prevlev in prevlist){
  fp <- tab %>%
    filter(prev  == prevlev) %>%
    pivot_wider(names_from=cov,values_from=combo) %>%
    arrange(factor(approach,levels=c("Cough+Xpert","Xpert","CXR"))) %>%
    select(!c(prev))
  colnames(fp) <- c("Algorithm","Rounds","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")
  write_xlsx(fp,paste("./Results_",prevlev,"/tablefp.xlsx",sep=""))
}

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 
testcolrange <- colorRampPalette(c("white","#5c5c5c"))
legendtemp <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
  geom_bar(position=position_dodge(),stat="identity",colour="black") +
  geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
  scale_x_discrete("") +
  scale_y_continuous("",expand=c(0,0),labels=scales::comma,lim=c(0,ylimit)) +
  scale_fill_manual("Round", values=testcolrange(durnum+1)[2:(durnum+1)]) +
  ggtitle("") +
  theme_classic(base_size=18) +
  theme(legend.position="bottom",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank()
  ) + #
  guides(fill=guide_legend(title.position="top",title.hjust=0.5))
fptpleg <- g_legend(legendtemp)

fptpfig <-  (fptplist[[ 9]] | fptplist[[5]] | fptplist[[1]]) /
            (fptplist[[11]] | fptplist[[7]] | fptplist[[3]]) /
            (fptplist[[12]] | fptplist[[8]] | fptplist[[4]]) /
            fptpleg +
              plot_layout(heights = c(2, 2, 2, 0.3))
ggsave(filename=paste("./Results/plot_fptp.tiff",sep=""),plot=fptpfig,width=60,height=40, unit="cm")

# Combined figure of TP, FP, FP:TP
fptplist <- list()
ylim     <- list()
fptpdata <- c()
prev <- 500
  # True positives by round
  load(paste("./Results_",prev,"/res/yr5results.Rdata",sep=""))
  
  plotlist <- list()
  ylimit <- max(yr5results[,5,c(1,3,5,7,9),,3])
  
  for (test in c(1:(algnum-1))){
    outdata <- c()
    if (test==1){
      testcol  <- colsympxpert
      testname <- "Cough+Xpert"
    }
    if (test==2){
      testcol  <- colcxrxpert
      testname <- "CXR+Xpert"
    }
    if (test==3){
      testcol  <- colxpert
      testname <- "Xpert"
    }
    if (test==4){
      testcol  <- colcxr
      testname <- "CXR"
    }
    testcolrange <- colorRampPalette(c("white",testcol))
    for (covn in c(1:covnum)){
      for (int in c(1:durnum)){
        
        if (covn== 1) cov <- '10%'
        if (covn== 2) cov <- '20%'
        if (covn== 3) cov <- '30%'
        if (covn== 4) cov <- '40%'
        if (covn== 5) cov <- '50%'
        if (covn== 6) cov <- '60%'
        if (covn== 7) cov <- '70%'
        if (covn== 8) cov <- '80%'
        if (covn== 9) cov <- '90%'
        if (covn==10) cov <- '100%'
        
        temp <- yr5results[covn,5,int*2-1,test,]
        outdata <- rbind(outdata,c(temp,covn,cov,int))
      }
    }
    outdata <- as.data.frame(outdata)
    colnames(outdata) <- c("lower","median","upper","covn","cov","int")
    outdata$int <- as.character(outdata$int)  
    outdata$median <- as.numeric(outdata$median)  
    outdata$lower  <- as.numeric(outdata$lower)  
    outdata$upper  <- as.numeric(outdata$upper)  
    if (test==1){
    outdataplot <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
      geom_bar(position=position_dodge(),stat="identity",colour="black") +
      geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
      scale_x_discrete("") +
      scale_y_continuous("Algorithm targeting sTB\n(Cough+Xpert)",expand=c(0,0),labels=scales::comma) + #,lim=c(0,ylimit)) +
      scale_fill_manual("Rounds", values=testcolrange(durnum+1)[2:(durnum+1)]) +
      ggtitle("True positives") +
      theme_classic(base_size=14) +
      theme(legend.position="none",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5),axis.text.x=element_blank()
      ) + #
      guides(fill=guide_legend(title.position="top",title.hjust=0.5))
    }else if (test==3){
        outdataplot <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
          geom_bar(position=position_dodge(),stat="identity",colour="black") +
          geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
          scale_x_discrete("") +
          scale_y_continuous("Algorithm targeting infectious TB\n(Xpert)",expand=c(0,0),labels=scales::comma) + #,lim=c(0,ylimit)) +
          scale_fill_manual("Rounds", values=testcolrange(durnum+1)[2:(durnum+1)]) +
          ggtitle("") +
          theme_classic(base_size=14) +
          theme(legend.position="none",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5),axis.text.x=element_blank()
          ) + #
          guides(fill=guide_legend(title.position="top",title.hjust=0.5))
    }else if (test==4){
      outdataplot <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
        geom_bar(position=position_dodge(),stat="identity",colour="black") +
        geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
        scale_x_discrete("Population coverage") +
        scale_y_continuous("Algorithm targeting all TB\n(CXR)",expand=c(0,0),labels=scales::comma) + #,lim=c(0,ylimit)) +
        scale_fill_manual("Rounds", values=testcolrange(durnum+1)[2:(durnum+1)]) +
        ggtitle("") +
        theme_classic(base_size=14) +
        theme(legend.position="none",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5)#,axis.title.y=element_blank()
        ) + #
        guides(fill=guide_legend(title.position="top",title.hjust=0.5))
    }
    plotlist[[test]] <- outdataplot
  }
  figtp <- (plotlist[[1]] / plotlist[[3]] / plotlist[[4]])
  ggsave(filename=paste("./Results_",prev,"/plot_tpbyround.tiff",sep=""),plot=figtp,width=20,height=20, unit="cm")
  
  # False positives by round
  load(paste("./Results_",prev,"/res/yr5results.Rdata",sep=""))
  
  plotlist <- list()
  ylimit <- max(yr5results[,5,c(2,4,6,8,10),,3])
  
  for (test in c(1:(algnum-1))){
    outdata <- c()
    if (test==1){
      testcol  <- colsympxpert
      testname <- "Cough+Xpert"
    }
    if (test==2){
      testcol  <- colcxrxpert
      testname <- "CXR+Xpert"
    }
    if (test==3){
      testcol  <- colxpert
      testname <- "Xpert"
    }
    if (test==4){
      testcol  <- colcxr
      testname <- "CXR"
    }
    testcolrange <- colorRampPalette(c("white",testcol))
    for (covn in c(1:covnum)){
      for (int in c(1:durnum)){
        
        if (covn== 1) cov <- '10%'
        if (covn== 2) cov <- '20%'
        if (covn== 3) cov <- '30%'
        if (covn== 4) cov <- '40%'
        if (covn== 5) cov <- '50%'
        if (covn== 6) cov <- '60%'
        if (covn== 7) cov <- '70%'
        if (covn== 8) cov <- '80%'
        if (covn== 9) cov <- '90%'
        if (covn==10) cov <- '100%'
        
        temp <- yr5results[covn,5,int*2,test,]
        outdata <- rbind(outdata,c(temp,covn,cov,int))
      }
    }
    outdata <- as.data.frame(outdata)
    colnames(outdata) <- c("lower","median","upper","covn","cov","int")
    outdata$int <- as.character(outdata$int)  
    outdata$median <- as.numeric(outdata$median)  
    outdata$lower  <- as.numeric(outdata$lower)  
    outdata$upper  <- as.numeric(outdata$upper)  
    if (test==1){
      outdataplot <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
        geom_bar(position=position_dodge(),stat="identity",colour="black") +
        geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
        scale_x_discrete("") +
        scale_y_continuous("",expand=c(0,0),labels=scales::comma) + #,lim=c(0,ylimit)) +
        scale_fill_manual("Rounds", values=testcolrange(durnum+1)[2:(durnum+1)]) +
        ggtitle("False positives") +
        theme_classic(base_size=14) +
        theme(legend.position="none",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5),axis.text.x=element_blank()
        ) + #
        guides(fill=guide_legend(title.position="top",title.hjust=0.5))
    }else if (test==3){
      outdataplot <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
        geom_bar(position=position_dodge(),stat="identity",colour="black") +
        geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
        scale_x_discrete("") +
        scale_y_continuous("",expand=c(0,0),labels=scales::comma) + #,lim=c(0,ylimit)) +
        scale_fill_manual("Rounds", values=testcolrange(durnum+1)[2:(durnum+1)]) +
        ggtitle("") +
        theme_classic(base_size=14) +
        theme(legend.position="none",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5),axis.text.x=element_blank()
        ) + #
        guides(fill=guide_legend(title.position="top",title.hjust=0.5))
    }else if (test==4){
      outdataplot <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
        geom_bar(position=position_dodge(),stat="identity",colour="black") +
        geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
        scale_x_discrete("Population coverage") +
        scale_y_continuous("",expand=c(0,0),labels=scales::comma) + #,lim=c(0,ylimit)) +
        scale_fill_manual("Rounds", values=testcolrange(durnum+1)[2:(durnum+1)]) +
        ggtitle("") +
        theme_classic(base_size=14) +
        theme(legend.position="none",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5)#,axis.title.y=element_blank()
        ) + #
        guides(fill=guide_legend(title.position="top",title.hjust=0.5))
    }
    plotlist[[test]] <- outdataplot
  }
  figfp <- (plotlist[[1]] / plotlist[[3]] / plotlist[[4]])
  ggsave(filename=paste("./Results_",prev,"/plot_fpbyround.tiff",sep=""),plot=figfp,width=20,height=20, unit="cm")
  
  #FP:TP by round
  load(paste("./Results_",prev,"/res/yr5results.Rdata",sep=""))
  
  plotlist <- list()
  for (test in c(1:(algnum-1))){
    outdata <- c()
    if (test==1){
      testcol  <- colsympxpert
      testname <- "Cough+Xpert"
    }
    if (test==2){
      testcol  <- colcxrxpert
      testname <- "CXR+Xpert"
    }
    if (test==3){
      testcol  <- colxpert
      testname <- "Xpert"
    }
    if (test==4){
      testcol  <- colcxr
      testname <- "CXR"
    }
    testcolrange <- colorRampPalette(c("white",testcol))
    if (prev == 250){
      if (test == 4)  ylim[[test]] <- max(yr5results[,5,c(11:15),test,2])*1.05
      else            ylim[[test]] <- max(yr5results[,5,c(11:15),test,3])*1.05
    }
    for (covn in c(1:covnum)){
      for (int in c(1:durnum)){
        
        if (covn== 1) cov <- '10%'
        if (covn== 2) cov <- '20%'
        if (covn== 3) cov <- '30%'
        if (covn== 4) cov <- '40%'
        if (covn== 5) cov <- '50%'
        if (covn== 6) cov <- '60%'
        if (covn== 7) cov <- '70%'
        if (covn== 8) cov <- '80%'
        if (covn== 9) cov <- '90%'
        if (covn==10) cov <- '100%'
        
        temp <- yr5results[covn,5,int+10,test,]
        fptpdata <- rbind(fptpdata,c(temp,testname,cov,int,prev))
        if (prev == 250 && test == 4 && cov == "100%" && int == 5){
          temp[3] <- temp[2]
        }
        outdata <- rbind(outdata,c(temp,covn,cov,int))
      }
    }
    outdata <- as.data.frame(outdata)
    colnames(outdata) <- c("lower","median","upper","covn","cov","int")
    outdata$int <- as.character(outdata$int)  
    outdata$median <- as.numeric(outdata$median)  
    outdata$lower  <- as.numeric(outdata$lower)  
    outdata$upper  <- as.numeric(outdata$upper)  
    outdataplot <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
      geom_bar(position=position_dodge(),stat="identity",colour="black") +
      geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
      scale_x_discrete("Population coverage") +
      scale_y_continuous("False positives per true positive",expand=c(0,0),labels=scales::comma) + #,lim=c(0,ylimit)) +
      scale_fill_manual("Rounds", values=testcolrange(durnum+1)[2:(durnum+1)]) +
      ggtitle(testname) +
      theme_classic(base_size=14) +
      theme(legend.position="none",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5)#,axis.title.y=element_blank()
      ) + #
      guides(fill=guide_legend(title.position="top",title.hjust=0.5))
    plotlist[[test]] <- outdataplot

    if (test==1){
      outdataplot <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
        geom_bar(position=position_dodge(),stat="identity",colour="black") +
        geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
        scale_x_discrete("") +
        scale_y_continuous("",expand=c(0,0),labels=scales::comma) + #,lim=c(0,ylimit)) +
        scale_fill_manual("Rounds", values=testcolrange(durnum+1)[2:(durnum+1)]) +
        ggtitle("False positives per true positive") +
        theme_classic(base_size=14) +
        theme(legend.position="none",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5),axis.text.x=element_blank()
        ) + #
        guides(fill=guide_legend(title.position="top",title.hjust=0.5))
    }else if (test==3){
      outdataplot <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
        geom_bar(position=position_dodge(),stat="identity",colour="black") +
        geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
        scale_x_discrete("") +
        scale_y_continuous("",expand=c(0,0),labels=scales::comma) + #,lim=c(0,ylimit)) +
        scale_fill_manual("Rounds", values=testcolrange(durnum+1)[2:(durnum+1)]) +
        ggtitle("") +
        theme_classic(base_size=14) +
        theme(legend.position="none",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5),axis.text.x=element_blank()
        ) + #
        guides(fill=guide_legend(title.position="top",title.hjust=0.5))
    }else if (test==4){
      outdataplot <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
        geom_bar(position=position_dodge(),stat="identity",colour="black") +
        geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
        scale_x_discrete("Population coverage") +
        scale_y_continuous("",expand=c(0,0),labels=scales::comma) + #,lim=c(0,ylimit)) +
        scale_fill_manual("Rounds", values=testcolrange(durnum+1)[2:(durnum+1)]) +
        ggtitle("") +
        theme_classic(base_size=14) +
        theme(legend.position="none",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5)#,axis.title.y=element_blank()
        ) + #
        guides(fill=guide_legend(title.position="top",title.hjust=0.5))
    }
    plotlist[[test]] <- outdataplot
  }
  
  figfptp <- (plotlist[[1]] / plotlist[[3]] / plotlist[[4]])
  ggsave(filename=paste("./Results_",prev,"/plot_fptpbyround2.tiff",sep=""),plot=figfptp,width=20,height=20, unit="cm")
  
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 
testcolrange <- colorRampPalette(c("white","#5c5c5c"))
legendtemp <- ggplot(filter(outdata,test != 2),aes(x=factor(cov,levels=c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')),y=median,fill=int)) +
  geom_bar(position=position_dodge(),stat="identity",colour="black") +
  geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
  scale_x_discrete("") +
  scale_y_continuous("",expand=c(0,0),labels=scales::comma) +
  scale_fill_manual("Round", values=testcolrange(durnum+1)[2:(durnum+1)]) +
  ggtitle("") +
  theme_classic(base_size=14) +
  theme(legend.position="bottom",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank()
  ) + #
  guides(fill=guide_legend(title.position="top",title.hjust=0.5))
fptpleg <- g_legend(legendtemp)

fig <- (figtp | figfp | figfptp) /
  fptpleg +
  plot_layout(heights = c(6, 0.2))
ggsave(filename=paste("./Results_",prev,"/plot_resbyround2.tiff",sep=""),plot=fig,width=60,height=30, unit="cm")