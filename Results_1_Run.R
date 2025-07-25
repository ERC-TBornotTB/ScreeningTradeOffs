#------------------------------------------------------------------------------#
# Results_1_Run.R                                                              #
# Analyse and present intervention results overall                             #
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
colnoint     <- '#5c5c5c'   #'#4A4E69'#'#283845' #'#899E8B'
colsympxpert <- '#BC3C29' #'#CB4C15' #'#1F4E79' '#4A4E69'
colcxrxpert  <- '#4A4E69'   #'#8E5572' '#C08497'
colxpert     <- '#E18727' #'#FFC325' #'#899E8B' '#899E8B'
colcxr       <- '#6F99AD' #'#537D8D' #'#4A4E69' #'#899E8B' #'#2C6E49' '#88665D' '#537D8D'

# Assign disease state colours
colsTB <- colorRampPalette(c("white",colcxr))(5)[[5]]  #'#8C2F39' 
colaTB <- colorRampPalette(c("white",colcxr))(5)[[3]]  #'#F29559'
colnTB <- colorRampPalette(c("white",colcxr))(5)[[2]]  #'#FFD166'
colnot <- colorRampPalette(c("white",colcxr))(5)[[5]]  #'#8AA29E'

# Set scenario parameters
pop      <- 100000  # Population
plotmos  <- 332     # Number of months of follow-up
datapts  <- 51      # Number of data points in model output
outnum   <- 29      # Number of outcomes 
algnum   <- 5       # Number of screening algorithms
covnum   <- 10      # Number of screening coverage levels
durnum   <- 5       # Number of screening years
prevlist <- c(250,500,1000)                             # List of prevalences
covlist  <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)  # List of coverages # 0.1,,1.0,
durlist  <- c(1,2,3,4,5)                                # List of durations

# Source functions and files
source("./Results_0_Functions.R")   

# Calculate outputs
for (prev in prevlist){

  # Create folder to store results
  dir.create(file.path(paste("./Results_",prev,sep="")), showWarnings=FALSE)

  # Set up output structure
  scrall      <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  conall      <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  scrposall   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  conposstb   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  conposatb   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  conposntb   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  conposinftb <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  conpostb    <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  conposnot   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  connegstb   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  connegatb   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  connegntb   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  conneginftb <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  connegtb    <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  connegnot   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  totnegstb   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  totnegatb   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  totnegntb   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  totneginftb <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  totnegtb    <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  totnegnot   <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  tbinc       <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  tbmort      <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  tbprev00    <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  tbprev10    <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  tbinfec     <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  tbepstb     <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  tbdeath     <- data.frame(none=0, sympxpert=0, cxrxpert=0, xpert=0, cxr=0)
  
  fullresults <- array(NA,dim=c(covnum,durnum,outnum,(algnum-1),3))
  
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
          
          # Number of individuals undergoing screening
          scrall[i,n]    <- sum(set[[n]][132,c("ScrNegS"  ,"ConNegS"  ,"ConPosS"  ,
                                               "ScrNegI"  ,"ConNegI"  ,"ConPosI"  ,
                                               "ScrNegnTB","ConNegnTB","ConPosnTB",
                                               "ScrNegaTB","ConNegaTB","ConPosaTB",
                                               "ScrNegsTB","ConNegsTB","ConPossTB",
                                               "ScrNegRn" ,"ConNegRn" ,"ConPosRn" ,
                                               "ScrNegRt" ,"ConNegRt" ,"ConPosRt" )])
  
          # Number of individuals undergoing confirmatory testing
          conall[i,n]    <- sum(set[[n]][132,c("ConNegS"  ,"ConPosS"  ,
                                               "ConNegI"  ,"ConPosI"  ,
                                               "ConNegnTB","ConPosnTB",
                                               "ConNegaTB","ConPosaTB",
                                               "ConNegsTB","ConPossTB",
                                               "ConNegRn" ,"ConPosRn" ,
                                               "ConNegRt" ,"ConPosRt" )])
          
          # Number of individuals with positive screening results
          scrposall[i,n] <- sum(set[[n]][132,c("ConPosS"  ,"ConNegS"  ,
                                               "ConPosI"  ,"ConNegI"  ,
                                               "ConPosnTB","ConNegnTB",
                                               "ConPosaTB","ConNegaTB",
                                               "ConPossTB","ConNegsTB",
                                               "ConPosRn" ,"ConNegRn" ,
                                               "ConPosRt" ,"ConNegRt" )])
  
          # Number of individuals with positive confirmatory results
          conposstb[i,n]   <- sum(set[[n]][132,c("ConPossTB")])
          conposatb[i,n]   <- sum(set[[n]][132,c("ConPosaTB")])
          conposntb[i,n]   <- sum(set[[n]][132,c("ConPosnTB")])
          conposinftb[i,n] <- sum(set[[n]][132,c("ConPossTB",
                                                 "ConPosaTB")])
          conpostb[i,n]    <- sum(set[[n]][132,c("ConPossTB",
                                                 "ConPosaTB",
                                                 "ConPosnTB")])
          conposnot[i,n]   <- sum(set[[n]][132,c("ConPosS",
                                                 "ConPosI",
                                                 "ConPosRn",
                                                 "ConPosRt")])
          
          # Number of individuals with negative confirmatory results
          connegstb[i,n]   <- sum(set[[n]][132,c("ConNegsTB")])
          connegatb[i,n]   <- sum(set[[n]][132,c("ConNegaTB")])
          connegntb[i,n]   <- sum(set[[n]][132,c("ConNegnTB")])
          conneginftb[i,n] <- sum(set[[n]][132,c("ConNegsTB",
                                                 "ConNegaTB")])
          connegtb[i,n]    <- sum(set[[n]][132,c("ConNegsTB",
                                                 "ConNegaTB",
                                                 "ConNegnTB")])
          connegnot[i,n]   <- sum(set[[n]][132,c("ConNegS",
                                                 "ConNegI",
                                                 "ConNegRn",
                                                 "ConNegRt")])
          
          # Number of individuals with negative screening or confirmatory results
          totnegstb[i,n]   <- sum(set[[n]][132,c("ScrNegsTB","ConNegsTB")])
          totnegatb[i,n]   <- sum(set[[n]][132,c("ScrNegaTB","ConNegaTB")])
          totnegntb[i,n]   <- sum(set[[n]][132,c("ScrNegnTB","ConNegnTB")])
          totneginftb[i,n] <- sum(set[[n]][132,c("ScrNegsTB","ConNegsTB",
                                                 "ScrNegaTB","ConNegaTB")])
          totnegtb[i,n]    <- sum(set[[n]][132,c("ScrNegsTB","ConNegsTB",
                                                 "ScrNegaTB","ConNegaTB",
                                                 "ScrNegnTB","ConNegnTB")])
          totnegnot[i,n]   <- sum(set[[n]][132,c("ScrNegS"  ,"ConNegS"  ,
                                                 "ScrNegI"  ,"ConNegI"  ,
                                                 "ScrNegRn" ,"ConNegRn" ,
                                                 "ScrNegRt" ,"ConNegRt" )])
  
          # TB burden indicators
          tbinc[i,n]    <- sum(set[[n]][121:132,c("Ic")])
          tbmort[i,n]   <- sum(set[[n]][121:132,c("MsTB")])
  
          tbinfec[i,n]  <- sum(set[[n]][1:132,c("Is","Ir","It")])
          tbepstb[i,n]  <- sum(set[[n]][1:132,c("Ic")])
          tbdeath[i,n]  <- sum(set[[n]][1:132,c("MsTB")])
          
          tbprev00[i,n] <- sum(set[[n]][  1,c("aTB","sTB","aTBx","sTBx")])
          tbprev10[i,n] <- sum(set[[n]][132,c("aTB","sTB","aTBx","sTBx")])
  
          # Generate incidence and mortality datasets
          fulldatainc  <- array(NA,dim=c(length(run_sets),algnum,plotmos,1))
          fulldatamort <- array(NA,dim=c(length(run_sets),algnum,plotmos,1))
          
        }
      }
      
      for (i in 1:length(run_sets)){    # Number of parameter sets
        for (j in 1:algnum){            # Number of interventions
          for (k in 1:plotmos){         # Number of years
            fulldatainc[i,j,k,1]  <- fulldata[i,j,k,c(30)]*12
            fulldatamort[i,j,k,1] <- fulldata[i,j,k,c(29)]*12
          }
        }
      }
          
      save(fulldatainc ,file=paste("./Results_",prev,"/res/fulldatainc_",cov,"_",int,".Rdata",sep=""))
      save(fulldatamort,file=paste("./Results_",prev,"/res/fulldatamort_",cov,"_",int,".Rdata",sep=""))

      ## SCREENING PERFORMANCE ##
      
      # Number of individuals with positive results: sTB
      fullresults[covn,int, 1,,] <- rbind(quantile(conposstb$sympxpert,probs=c(0.025,0.5,0.975)),
                                          quantile(conposstb$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                          quantile(conposstb$xpert    ,probs=c(0.025,0.5,0.975)),
                                          quantile(conposstb$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of individuals with positive results: aTB
      fullresults[covn,int, 2,,] <- rbind(quantile(conposatb$sympxpert,probs=c(0.025,0.5,0.975)),
                                          quantile(conposatb$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                          quantile(conposatb$xpert    ,probs=c(0.025,0.5,0.975)),
                                          quantile(conposatb$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of individuals with positive results: nTB
      fullresults[covn,int, 3,,] <- rbind(quantile(conposntb$sympxpert,probs=c(0.025,0.5,0.975)),
                                          quantile(conposntb$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                          quantile(conposntb$xpert    ,probs=c(0.025,0.5,0.975)),
                                          quantile(conposntb$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of individuals with positive results: infectious TB
      fullresults[covn,int, 4,,] <- rbind(quantile(conposinftb$sympxpert,probs=c(0.025,0.5,0.975)),
                                          quantile(conposinftb$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                          quantile(conposinftb$xpert    ,probs=c(0.025,0.5,0.975)),
                                          quantile(conposinftb$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of individuals with positive results: TB
      fullresults[covn,int, 5,,] <- rbind(quantile(conpostb$sympxpert,probs=c(0.025,0.5,0.975)),
                                          quantile(conpostb$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                          quantile(conpostb$xpert    ,probs=c(0.025,0.5,0.975)),
                                          quantile(conpostb$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of individuals with positive results: notTB
      fullresults[covn,int, 6,,] <- rbind(quantile(conposnot$sympxpert,probs=c(0.025,0.5,0.975)),
                                          quantile(conposnot$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                          quantile(conposnot$xpert    ,probs=c(0.025,0.5,0.975)),
                                          quantile(conposnot$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of individuals with negative results: sTB
      fullresults[covn,int, 7,,] <- rbind(quantile(totnegstb$sympxpert,probs=c(0.025,0.5,0.975)),
                                          quantile(totnegstb$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                          quantile(totnegstb$xpert    ,probs=c(0.025,0.5,0.975)),
                                          quantile(totnegstb$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of individuals with negative results: aTB
      fullresults[covn,int, 8,,] <- rbind(quantile(totnegatb$sympxpert,probs=c(0.025,0.5,0.975)),
                                          quantile(totnegatb$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                          quantile(totnegatb$xpert    ,probs=c(0.025,0.5,0.975)),
                                          quantile(totnegatb$cxr      ,probs=c(0.025,0.5,0.975)))
  
      # Number of individuals with negative results: nTB
      fullresults[covn,int, 9,,] <- rbind(quantile(totnegntb$sympxpert,probs=c(0.025,0.5,0.975)),
                                          quantile(totnegntb$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                          quantile(totnegntb$xpert    ,probs=c(0.025,0.5,0.975)),
                                          quantile(totnegntb$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of individuals with negative results: infectious TB
      fullresults[covn,int,10,,] <- rbind(quantile(totneginftb$sympxpert,probs=c(0.025,0.5,0.975)),
                                          quantile(totneginftb$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                          quantile(totneginftb$xpert    ,probs=c(0.025,0.5,0.975)),
                                          quantile(totneginftb$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of individuals with negative results: TB
      fullresults[covn,int,11,,] <- rbind(quantile(totnegtb$sympxpert,probs=c(0.025,0.5,0.975)),
                                          quantile(totnegtb$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                          quantile(totnegtb$xpert    ,probs=c(0.025,0.5,0.975)),
                                          quantile(totnegtb$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of individuals with negative results: notTB
      fullresults[covn,int,12,,] <- rbind(quantile(totnegnot$sympxpert,probs=c(0.025,0.5,0.975)),
                                          quantile(totnegnot$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                          quantile(totnegnot$xpert    ,probs=c(0.025,0.5,0.975)),
                                          quantile(totnegnot$cxr      ,probs=c(0.025,0.5,0.975)))
      
      ## SCREENING EFFICIENCY ##
      
      # Screening efficiency: Number needed to screen
  
      # Post-screening prevalence: sTB
      pspstb <- cbind(conposstb,connegstb,scrposall)
      names(pspstb) <- c(paste0("conpos_",names(conposstb)),paste0("conneg_",names(connegstb)),paste0("scrposall_",names(scrposall)))
      pspstb$scrposall_xpert <- pop
      fullresults[covn,int,13,,] <- rbind(quantile((pspstb$conpos_sympxpert+pspstb$conneg_sympxpert)/pspstb$scrposall_sympxpert,probs=c(0.025,0.5,0.975)),
                                          quantile((pspstb$conpos_cxrxpert +pspstb$conneg_cxrxpert )/pspstb$scrposall_cxrxpert ,probs=c(0.025,0.5,0.975)),
                                          quantile((pspstb$conpos_xpert    +pspstb$conneg_xpert    )/pspstb$scrposall_xpert    ,probs=c(0.025,0.5,0.975)),
                                          quantile((pspstb$conpos_cxr      +pspstb$conneg_cxr      )/pspstb$scrposall_cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Post-screening prevalence: aTB
      pspatb <- cbind(conposatb,connegatb,scrposall)
      names(pspatb) <- c(paste0("conpos_",names(conposatb)),paste0("conneg_",names(connegatb)),paste0("scrposall_",names(scrposall)))
      pspatb$scrposall_xpert <- pop
      fullresults[covn,int,14,,] <- rbind(quantile((pspatb$conpos_sympxpert+pspatb$conneg_sympxpert)/pspatb$scrposall_sympxpert,probs=c(0.025,0.5,0.975)),
                                          quantile((pspatb$conpos_cxrxpert +pspatb$conneg_cxrxpert )/pspatb$scrposall_cxrxpert ,probs=c(0.025,0.5,0.975)),
                                          quantile((pspatb$conpos_xpert    +pspatb$conneg_xpert    )/pspatb$scrposall_xpert    ,probs=c(0.025,0.5,0.975)),
                                          quantile((pspatb$conpos_cxr      +pspatb$conneg_cxr      )/pspatb$scrposall_cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Post-screening prevalence: nTB
      pspntb <- cbind(conposntb,connegntb,scrposall)
      names(pspntb) <- c(paste0("conpos_",names(conposntb)),paste0("conneg_",names(connegntb)),paste0("scrposall_",names(scrposall)))
      pspntb$scrposall_xpert <- pop
      fullresults[covn,int,15,,] <- rbind(quantile((pspntb$conpos_sympxpert+pspntb$conneg_sympxpert)/pspntb$scrposall_sympxpert,probs=c(0.025,0.5,0.975)),
                                          quantile((pspntb$conpos_cxrxpert +pspntb$conneg_cxrxpert )/pspntb$scrposall_cxrxpert ,probs=c(0.025,0.5,0.975)),
                                          quantile((pspntb$conpos_xpert    +pspntb$conneg_xpert    )/pspntb$scrposall_xpert    ,probs=c(0.025,0.5,0.975)),
                                          quantile((pspntb$conpos_cxr      +pspntb$conneg_cxr      )/pspntb$scrposall_cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Post-screening prevalence: infectious TB
      pspinftb <- cbind(conposinftb,conneginftb,scrposall)
      names(pspinftb) <- c(paste0("conpos_",names(conposinftb)),paste0("conneg_",names(conneginftb)),paste0("scrposall_",names(scrposall)))
      pspinftb$scrposall_xpert <- pop
      fullresults[covn,int,16,,] <- rbind(quantile((pspinftb$conpos_sympxpert+pspinftb$conneg_sympxpert)/pspinftb$scrposall_sympxpert,probs=c(0.025,0.5,0.975)),
                                          quantile((pspinftb$conpos_cxrxpert +pspinftb$conneg_cxrxpert )/pspinftb$scrposall_cxrxpert ,probs=c(0.025,0.5,0.975)),
                                          quantile((pspinftb$conpos_xpert    +pspinftb$conneg_xpert    )/pspinftb$scrposall_xpert    ,probs=c(0.025,0.5,0.975)),
                                          quantile((pspinftb$conpos_cxr      +pspinftb$conneg_cxr      )/pspinftb$scrposall_cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Post-screening prevalence: TB
      psptb <- cbind(conpostb,connegtb,scrposall)
      names(psptb) <- c(paste0("conpos_",names(conpostb)),paste0("conneg_",names(connegtb)),paste0("scrposall_",names(scrposall)))
      psptb$scrposall_xpert <- pop
      fullresults[covn,int,17,,] <- rbind(quantile((psptb$conpos_sympxpert+psptb$conneg_sympxpert)/psptb$scrposall_sympxpert,probs=c(0.025,0.5,0.975)),
                                          quantile((psptb$conpos_cxrxpert +psptb$conneg_cxrxpert )/psptb$scrposall_cxrxpert ,probs=c(0.025,0.5,0.975)),
                                          quantile((psptb$conpos_xpert    +psptb$conneg_xpert    )/psptb$scrposall_xpert    ,probs=c(0.025,0.5,0.975)),
                                          quantile((psptb$conpos_cxr      +psptb$conneg_cxr      )/psptb$scrposall_cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of screening tests
      scrall$xpert <- 0
      scrall$min   <- 0
      fullresults[covn,int,18,,] <- rbind(quantile(scrall$sympxpert,probs=c(0.025,0.5,0.975)),
                                          quantile(scrall$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                          quantile(scrall$xpert    ,probs=c(0.025,0.5,0.975)),
                                          quantile(scrall$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Number of confirmatory tests
      fullresults[covn,int,19,,] <- rbind(quantile(conall$sympxpert,probs=c(0.025,0.5,0.975)),
                                          quantile(conall$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                          quantile(conall$xpert    ,probs=c(0.025,0.5,0.975)),
                                          quantile(conall$cxr      ,probs=c(0.025,0.5,0.975)))
      
      ## IMPACT ON DISEASE BURDEN ##
      
      # Disease burden: sTB incidence in 2035
      fullresults[covn,int,20,,] <- rbind(quantile(with(tbinc,(none-sympxpert)/none*100),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbinc,(none-cxrxpert )/none*100),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbinc,(none-xpert    )/none*100),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbinc,(none-cxr      )/none*100),probs=c(0.025,0.5,0.975)))
      
      # Disease burden: sTB mortality in 2035
      fullresults[covn,int,21,,] <- rbind(quantile(with(tbmort,(none-sympxpert)/none*100),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbmort,(none-cxrxpert )/none*100),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbmort,(none-xpert    )/none*100),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbmort,(none-cxr      )/none*100),probs=c(0.025,0.5,0.975)))
      
      # Disease burden: infTB prevalence prior to screening
      fullresults[covn,int,22,,] <- rbind(quantile(with(tbprev00,(none-sympxpert)/none*100),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbprev00,(none-cxrxpert )/none*100),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbprev00,(none-xpert    )/none*100),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbprev00,(none-cxr      )/none*100),probs=c(0.025,0.5,0.975)))
      
      # Disease burden: infTB prevalence in 2035
      fullresults[covn,int,23,,] <- rbind(quantile(with(tbprev10,(none-sympxpert)/none*100),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbprev10,(none-cxrxpert )/none*100),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbprev10,(none-xpert    )/none*100),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbprev10,(none-cxr      )/none*100),probs=c(0.025,0.5,0.975)))
      
      # Disease burden: sTB episodes (n)
      fullresults[covn,int,24,,] <- rbind(quantile(tbepstb$sympxpert,probs=c(0.025,0.5,0.975)),
                                          quantile(tbepstb$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                          quantile(tbepstb$xpert    ,probs=c(0.025,0.5,0.975)),
                                          quantile(tbepstb$cxr      ,probs=c(0.025,0.5,0.975)))
      
      # Disease burden: sTB episodes averted (n)
      fullresults[covn,int,25,,] <- rbind(quantile(with(tbepstb,(none-sympxpert)),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbepstb,(none-cxrxpert )),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbepstb,(none-xpert    )),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbepstb,(none-cxr      )),probs=c(0.025,0.5,0.975)))
      
      # Disease burden: sTB episodes averted (%)
      fullresults[covn,int,26,,] <- rbind(quantile(with(tbepstb,(none-sympxpert)/none*100),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbepstb,(none-cxrxpert )/none*100),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbepstb,(none-xpert    )/none*100),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbepstb,(none-cxr      )/none*100),probs=c(0.025,0.5,0.975)))
  
      # Disease burden: sTB deaths (n)
      fullresults[covn,int,27,,] <- rbind(quantile(tbdeath$sympxpert,probs=c(0.025,0.5,0.975)),
                                          quantile(tbdeath$cxrxpert ,probs=c(0.025,0.5,0.975)),
                                          quantile(tbdeath$xpert    ,probs=c(0.025,0.5,0.975)),
                                          quantile(tbdeath$cxr      ,probs=c(0.025,0.5,0.975)))
        
      # Disease burden: sTB deaths averted (n)
      fullresults[covn,int,28,,] <- rbind(quantile(with(tbdeath,(none-sympxpert)),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbdeath,(none-cxrxpert )),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbdeath,(none-xpert    )),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbdeath,(none-cxr      )),probs=c(0.025,0.5,0.975)))
      
      # Disease burden: sTB deaths averted (%)
      fullresults[covn,int,29,,] <- rbind(quantile(with(tbdeath,(none-sympxpert)/none*100),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbdeath,(none-cxrxpert )/none*100),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbdeath,(none-xpert    )/none*100),probs=c(0.025,0.5,0.975)),
                                          quantile(with(tbdeath,(none-cxr      )/none*100),probs=c(0.025,0.5,0.975)))
      
      save(fullresults,file=paste("./Results_",prev,"/res/fullresults.Rdata",sep=""))
    }
  }  
  
  for (covn in c(1:covnum)){
    for (int in c(1:durnum)){
      for (out in c(1:outnum)){
        temp <- array(0,dim=c(4,3))
        temp <- fullresults[covn,int,out,,]
        write.csv(temp,file=paste("./Results_",prev,"/res/res_",out,"_",covn,"_",int,".csv",sep=""))
      }
    }
  }
} 

# Plot all outcomes
for (prev in prevlist){ 
  load(file=paste("./Results_",prev,"/res/fullresults.Rdata",sep=""))
  
  for (out in c(1:outnum)){
    if (out== 1)  measname <- "Number of individuals with positive results: sTB"
    if (out== 2)  measname <- "Number of individuals with positive results: aTB"
    if (out== 3)  measname <- "Number of individuals with positive results: nTB"
    if (out== 4)  measname <- "Number of individuals with positive results: infectious TB"
    if (out== 5)  measname <- "Number of individuals with positive results: TB"
    if (out== 6)  measname <- "Number of individuals with positive results: notTB"
    if (out== 7)  measname <- "Number of individuals with negative results: sTB"
    if (out== 8)  measname <- "Number of individuals with negative results: aTB"
    if (out== 9)  measname <- "Number of individuals with negative results: nTB"
    if (out==10)  measname <- "Number of individuals with negative results: infectious TB"
    if (out==11)  measname <- "Number of individuals with negative results: TB"
    if (out==12)  measname <- "Number of individuals with negative results: notTB"
    if (out==13)  measname <- "Post-screening prevalence: sTB"
    if (out==14)  measname <- "Post-screening prevalence: aTB"
    if (out==15)  measname <- "Post-screening prevalence: nTB"
    if (out==16)  measname <- "Post-screening prevalence: infectious TB"
    if (out==17)  measname <- "Post-screening prevalence: TB"
    if (out==18)  measname <- "Number of screening tests"
    if (out==19)  measname <- "Number of confirmatory tests"
    if (out==20)  measname <- "sTB incidence in 2035"
    if (out==21)  measname <- "sTB mortality in 2035"
    if (out==22)  measname <- "infTB prevalence prior to screening"
    if (out==23)  measname <- "infTB prevalence in 2035"
    if (out==24)  measname <- "sTB episodes (n)"
    if (out==25)  measname <- "sTB episodes averted (n)"
    if (out==26)  measname <- "sTB episodes averted (%)"
    if (out==27)  measname <- "sTB deaths (n)"
    if (out==28)  measname <- "sTB deaths averted (n)"
    if (out==29)  measname <- "sTB deaths averted (%)"

    plotlist <- list()
    ylimit <- max(fullresults[,,out,,3])
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
        for (int in durlist){
          if (covn== 1)  cov <- 0.1
          if (covn== 2)  cov <- 0.2
          if (covn== 3)  cov <- 0.3
          if (covn== 4)  cov <- 0.4
          if (covn== 5)  cov <- 0.5
          if (covn== 6)  cov <- 0.6
          if (covn== 7)  cov <- 0.7
          if (covn== 8)  cov <- 0.8
          if (covn== 9)  cov <- 0.9
          if (covn==10)  cov <- 1.0
          testn   <- test
          temp    <- t(fullresults[covn,int,out,testn,])
          outdata <- rbind(outdata,c(temp,covn,cov,int))
        }
      }
      outdata <- as.data.frame(outdata)
      colnames(outdata) <- c("lower","median","upper","covn","cov","int")
      outdata$cov <- as.character(outdata$cov)
      outdata$int <- as.character(outdata$int)
      outdataplot <- ggplot(filter(outdata,test != 2),aes(x=cov,y=median,fill=int)) +
                     geom_bar(position=position_dodge(),stat="identity",colour="black") +
                     geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(.9)) +
                     scale_x_discrete("Coverage") +
                     scale_y_continuous(measname,expand=c(0,0),labels=scales::comma,lim=c(0,ylimit)) +
                     scale_fill_manual("Rounds", values=testcolrange(durnum+1)[2:(durnum+1)]) +
                     ggtitle(testname) +
                     theme_classic() +
                     theme(legend.position="bottom",legend.key.size=unit(0.25,'cm'),plot.title=element_text(hjust=0.5),axis.title.y=element_blank()) + #
                     guides(fill=guide_legend(title.position="top",title.hjust=0.5))
      plotlist[[test]] <- outdataplot
    }
    fig <- grid.arrange(plotlist[[1]],plotlist[[3]],plotlist[[4]],ncol=3,top=measname)
    ggsave(filename=paste("./Results_",prev,"/plotv2_",out,".tiff",sep=""),plot=fig,width=60,height=15, unit="cm")
  }
}

f <- chromote::default_chromote_object() #get the f object
f$close()

# Number of individuals with true and false positive and negative results
postaball <- c()
for (prev in prevlist){
  for (covn in c(1:length(covlist))){
  for (int in durlist){

  load(file=paste("./Results_",prev,"/res/fullresults.Rdata",sep=""))
  
  # Positive results
  posres <- as.data.frame(rbind(fullresults[covn,int,1,c(1,3:4),],fullresults[covn,int,2,c(1,3:4),],fullresults[covn,int,3,c(1,3:4),],
                                fullresults[covn,int,4,c(1,3:4),],fullresults[covn,int,5,c(1,3:4),],fullresults[covn,int,6,c(1,3:4),]))
  colnames(posres) <- c("low","med","high")
  posres$state     <- c(rep("sTB",(algnum-2)),rep("aTB",(algnum-2)),rep("nTB",(algnum-2)),rep("infectious TB",(algnum-2)),rep("all TB",(algnum-2)),rep("not TB",(algnum-2)))
  posres$approach  <- c(rep(c("sympxpert","xpert","cxr")))

  plottruepos <- posres
  plottruepos$approach <- factor(plottruepos$approach,levels=c("sympxpert","xpert","cxr"),labels=c("Cough+Xpert","Xpert","CXR"))
  plottruepos$state <- factor(plottruepos$state,levels=c("sTB","aTB","nTB","infectious TB","all TB","not TB"))
  
  tabtruepostemp <- plottruepos[c(1:15),]
  tabtruepostemp$state <- factor(tabtruepostemp$state,levels=c("sTB","aTB","nTB","infectious TB","all TB","not TB"))
  for (i in levels(plottruepos$approach)){
    tempmed  <- 0
    templow  <- 0
    temphigh <- 0
    for (j in 1:nrow(tabtruepostemp)){
      if (tabtruepostemp$approach[j] == i){
        tempmed  <- tempmed  + as.numeric(tabtruepostemp$med[j])
        templow  <- templow  + as.numeric(tabtruepostemp$low[j])
        temphigh <- temphigh + as.numeric(tabtruepostemp$high[j])
      }
    }
    #tabtruepostemp <- rbind(tabtruepostemp,c(templow,tempmed,temphigh,"Total",i,"",""))
  }
  for (i in 1:nrow(tabtruepostemp)){
    tabtruepostemp$res[i] <- paste(scales::comma(round(as.numeric(tabtruepostemp$med[i]),  digits=0)), " (", 
                                   scales::comma(round(as.numeric(tabtruepostemp$low[i]),  digits=0)), "-",
                                   scales::comma(round(as.numeric(tabtruepostemp$high[i]), digits=0)), ")", sep="")
  }
  truepostab <- tabtruepostemp %>%
    filter(state != "not TB") %>%
    arrange(approach,factor(state,c("sTB","aTB","nTB","infectious TB","all TB"))) %>%
    relocate(5,4,6) %>%
    select(approach,state,res) %>%
    pivot_wider(names_from=approach,values_from=res)
  colnames(truepostab) <- c("State","Cough+Xpert\nmedian (95% UI)","Xpert\nmedian (95% UI)","CXR\nmedian (95% UI)")

  truepostab %>% 
    gt() %>% 
    tab_header(title = md("Number of individuals with true positive results by screening algorithm")) %>%
    cols_align(align = "center") %>%
    opt_table_font(font = 'Open sans') %>%
    tab_style(
      style=cell_text(whitespace="pre"),
      locations=cells_column_labels()) %>%
    gtsave(filename=paste("./Results_",prev,"/table_truepos_",covn,"_",int,".png",sep=""), vwidth = 1600, vheight = 600)

  plotfalsepos <- posres[posres$state == "not TB", ]
  plotfalsepos$approach <- factor(plotfalsepos$approach, levels=c("sympxpert","cxrxpert","xpert","cxr",""), labels=c("Cough+Xpert","CXR+Xpert","Xpert","CXR",""))
  
  tabfalsepostemp <- plotfalsepos
  for (i in 1:nrow(tabfalsepostemp)){
    tabfalsepostemp$res[i] <- paste(scales::comma(round(as.numeric(tabfalsepostemp$med[i]),  digits=0)), " (", 
                                    scales::comma(round(as.numeric(tabfalsepostemp$low[i]),  digits=0)), "-",
                                    scales::comma(round(as.numeric(tabfalsepostemp$high[i]), digits=0)), ")", sep="")
  }
  falsepostab <- tabfalsepostemp %>%
    relocate(5,4,6) %>%
    select(approach,state,res) %>%
    pivot_wider(names_from=approach,values_from=res)
  colnames(falsepostab) <- c("State","Cough+Xpert\nmedian (95% UI)","Xpert\nmedian (95% UI)","CXR\nmedian (95% UI)")
  
  falsepostab %>% 
   gt() %>% 
   tab_header(title = md("Number of individuals with false positive results by screening algorithm")) %>%
   cols_align(align = "center") %>%
   opt_table_font(font = 'Open sans')  %>%
   tab_style(
     style=cell_text(whitespace="pre"),
     locations=cells_column_labels()) %>%
   gtsave(filename = paste("./Results_",prev,"/table_falsepos_",covn,"_",int,".png",sep=""), vwidth = 1600, vheight = 600)
   
  postab <- rbind(truepostab,falsepostab)
   
  postab %>% 
    gt() %>% 
    tab_header(title = md("Number of individuals with positive results")) %>%
    cols_align(align = "center") %>%
    opt_table_font(font = 'Open sans')  %>%
    cols_width(everything() ~ px(200)) %>%
    cols_width(State ~ px(150)) %>%
    tab_style(
      style=cell_text(whitespace="pre"),
      locations=cells_column_labels()) %>%
    gtsave(filename = paste("./Results_",prev,"/table_pos_",covn,"_",int,".png",sep=""), vwidth = 1600, vheight = 600)
  
  posrestab <- plottruepos %>%
    filter(state %in% c("all TB","not TB")) %>%
    mutate(prev = prev, cov = cov, covn = covn, int = int)
  postaball <- rbind(postaball,posrestab)
  }
  }
}  

tab <- postaball %>%
  mutate(cov   = paste(covn*10,"%",sep=""),
         combo = paste(round(med ,digits=0), " \n(", 
                       round(low ,digits=0), "-",
                       round(high,digits=0), ")", sep="")) %>%
  select(state,approach,prev,cov,int,combo)
for (prevlev in prevlist){
  truetab <- tab %>%
    filter(state == "all TB") %>%
    filter(prev  == prevlev) %>%
    pivot_wider(names_from=cov,values_from=combo) %>%
    arrange(factor(approach,levels=c("Cough+Xpert","Xpert","CXR"))) %>%
    select(!c(state,prev))
  colnames(truetab) <- c("Algorithm","Rounds","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")
  write_xlsx(truetab,paste("./Results_",prevlev,"/tabletruepos.xlsx",sep=""))
}
for (prevlev in prevlist){
  falsetab <- tab %>%
    filter(state == "not TB") %>%
    filter(prev  == prevlev) %>%
    pivot_wider(names_from=cov,values_from=combo) %>%
    arrange(factor(approach,levels=c("Cough+Xpert","Xpert","CXR"))) %>%
    select(!c(state,prev))
  colnames(falsetab) <- c("Algorithm","Rounds","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")
  write_xlsx(falsetab,paste("./Results_",prevlev,"/tablefalsepos.xlsx",sep=""))
}

# Number of individuals with true and false positive and negative results
for (prev in prevlist){
  covn <- max(covlist)*10
  int  <- min(durlist)
  
  load(file=paste("./Results_",prev,"/res/fullresults.Rdata",sep=""))
  
  # Positive results
  posres <- as.data.frame(rbind(fullresults[covn,int,1,,],fullresults[covn,int,2,,],fullresults[covn,int,3,,],fullresults[covn,int,6,,]))
  colnames(posres) <- c("low","med","high")
  posres$state     <- c(rep("sTB",(algnum-1)),rep("aTB",(algnum-1)),rep("nTB",(algnum-1)),rep("not TB",(algnum-1)))
  posres$approach  <- c(rep(c("sympxpert","cxrxpert","xpert","cxr")))
  
  posres$lowadj <- 0
  posres$lowadj[1:(algnum-1)]                    <- posres$low[1:(algnum-1)]
  posres$lowadj[algnum:((algnum-1)*2)]           <- posres$low[algnum:((algnum-1)*2)]           + posres$med[1:(algnum-1)]
  posres$lowadj[((algnum-1)*2+1):((algnum-1)*3)] <- posres$low[((algnum-1)*2+1):((algnum-1)*3)] + posres$med[algnum:((algnum-1)*2)] + posres$med[1:(algnum-1)]
  posres$lowadj[((algnum-1)*3+1):((algnum-1)*4)] <- posres$low[((algnum-1)*3+1):((algnum-1)*4)]
  
  posres$highadj <- 0
  posres$highadj[1:(algnum-1)]                    <- posres$high[1:(algnum-1)]
  posres$highadj[algnum:((algnum-1)*2)]           <- posres$high[algnum:((algnum-1)*2)]           + posres$med[1:(algnum-1)]
  posres$highadj[((algnum-1)*2+1):((algnum-1)*3)] <- posres$high[((algnum-1)*2+1):((algnum-1)*3)] + posres$med[algnum:((algnum-1)*2)] + posres$med[1:(algnum-1)]
  posres$highadj[((algnum-1)*3+1):((algnum-1)*4)] <- posres$high[((algnum-1)*3+1):((algnum-1)*4)]
  
  plottruepos <- posres
  plottruepos[c(13:16),c(1:3,6:7)] <- 0
  plottruepos$approach <- factor(plottruepos$approach,levels=c("sympxpert","cxrxpert","xpert","cxr"),labels=c("Cough+Xpert","CXR+Xpert","Xpert","CXR"))
  plottruepos$state <- factor(plottruepos$state,levels=c("sTB","aTB","nTB","not TB"))
  
  ymax <- max(plottruepos$highadj)*1.05
  truepos <- ggplot(filter(plottruepos,approach != "CXR+Xpert"),aes(x=approach,y=med,fill=approach:state)) +
    geom_bar(position=position_stack(reverse = TRUE),stat="identity") +
    geom_errorbar(aes(ymin=lowadj,ymax=highadj),width=0.5,size=0.5) +
    scale_x_discrete("Algorithm") + 
    scale_y_continuous("Number of individuals",expand=c(0,0)) + 
    ggtitle("True positives") +
    scale_fill_manual(values=c(colorRampPalette(c("white",colsympxpert))(5)[[5]],colorRampPalette(c("white",colsympxpert))(5)[[3]],colorRampPalette(c("white",colsympxpert))(5)[[2]],colorRampPalette(c("white",colsympxpert))(5)[[4]],
                               colorRampPalette(c("white",colxpert))(5)[[5]]    ,colorRampPalette(c("white",colxpert))(5)[[3]]    ,colorRampPalette(c("white",colxpert))(5)[[2]]    ,colorRampPalette(c("white",colxpert))(5)[[4]]    ,
                               colorRampPalette(c("white",colcxr))(5)[[5]]      ,colorRampPalette(c("white",colcxr))(5)[[3]]      ,colorRampPalette(c("white",colcxr))(5)[[2]]      ,colorRampPalette(c("white",colcxr))(5)[[4]]      ,
                               colorRampPalette(c("white",colcxrxpert))(5)[[5]] ,colorRampPalette(c("white",colcxrxpert))(5)[[3]] ,colorRampPalette(c("white",colcxrxpert))(5)[[2]] ,colorRampPalette(c("white",colcxrxpert))(5)[[4]] )) + 
    theme_classic(base_size=16) +
    labs(fill="Algorithm &\ndisease state") +
    guides(fill=guide_legend(nrow=3,byrow=TRUE)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(filename=paste("./Results_",prev,"/plot_truepos_",covn,"_",int,".tiff",sep=""),plot=truepos,width=30,height=15, unit="cm")
  
  tabtruepostemp <- plottruepos[c(1:12),]
  tabtruepostemp$state <- factor(tabtruepostemp$state,levels=c("sTB","aTB","nTB","not TB","Total"))
  for (i in levels(plottruepos$approach)){
    tempmed  <- 0
    templow  <- 0
    temphigh <- 0
    for (j in 1:nrow(tabtruepostemp)){
      if (tabtruepostemp$approach[j] == i){
        tempmed  <- tempmed  + as.numeric(tabtruepostemp$med[j])
        templow  <- templow  + as.numeric(tabtruepostemp$low[j])
        temphigh <- temphigh + as.numeric(tabtruepostemp$high[j])
      }
    }
    tabtruepostemp <- rbind(tabtruepostemp,c(templow,tempmed,temphigh,"Total",i,"",""))
  }
  for (i in 1:nrow(tabtruepostemp)){
    tabtruepostemp$res[i] <- paste(round(as.numeric(tabtruepostemp$med[i]),  digits=0), " (", 
                                   round(as.numeric(tabtruepostemp$low[i]),  digits=0), "-",
                                   round(as.numeric(tabtruepostemp$high[i]), digits=0), ")", sep="")
  }
  truepostab <- tabtruepostemp %>%
    filter(state != "not TB") %>%
    arrange(approach,factor(state,c("sTB","aTB","nTB"))) %>%
    relocate(5,4,8) %>%
    select(approach,state,res)
  colnames(truepostab) <- c("Algorithm", "State", "Median (95% UI)")
  truepostab$Algorithm <- factor(truepostab$Algorithm,levels=c("Cough+Xpert","CXR+Xpert","Xpert","CXR",""),labels=c("Cough+Xpert","CXR+Xpert","Xpert","CXR",""))
  truepostab[c(2:4,6:8,10:12,14:16),1] <- ""
  truepostab %>% 
    gt() %>% 
    tab_header(title = md("Number of individuals with true positive results by screening algorithm")) %>%
    cols_align(align = "center") %>%
    opt_table_font(font = 'Open sans') %>%
    gtsave(filename = paste("./Results_",prev,"/table_truepos_",covn,"_",int,".png",sep=""), vwidth = 1600, vheight = 600)
  
  
  plotfalsepos <- posres[posres$state == "not TB", ]
  plotfalsepos$approach <- factor(plotfalsepos$approach, levels=c("sympxpert","cxrxpert","xpert","cxr",""), labels=c("Cough+Xpert","CXR+Xpert","Xpert","CXR",""))
  
  ymax <- max(plotfalsepos$highadj)*1.05
  falsepos <- ggplot(filter(plotfalsepos,approach != "CXR+Xpert"),aes(x=approach,y=med,fill=approach)) +
    geom_bar(position=position_stack(reverse = TRUE),stat="identity") +
    geom_errorbar(aes(ymin=lowadj,ymax=highadj),width=0.5,size=0.5) +
    scale_x_discrete("Algorithm") + 
    scale_y_continuous("Number of individuals",expand=c(0,0)) + 
    ggtitle("False positives") +
    scale_fill_manual(values=c(colorRampPalette(c("white",colsympxpert))(5)[[4]],
                               colorRampPalette(c("white",colxpert))(5)[[4]], 
                               colorRampPalette(c("white",colcxr))(5)[[4]],
                               colorRampPalette(c("white",colcxrxpert))(5)[[4]])) + 
    theme_classic(base_size=16) +
    labs(fill="Disease state") +
    guides(fill=guide_legend(nrow=3,byrow=TRUE)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(filename=paste("./Results_",prev,"/plot_falsepos_",covn,"_",int,".tiff",sep=""),plot=falsepos,width=30,height=15, unit="cm")
  
  tabfalsepostemp <- plotfalsepos
  for (i in 1:nrow(tabfalsepostemp)){
    tabfalsepostemp$res[i] <- paste(round(as.numeric(tabfalsepostemp$med[i]),  digits=0), " (", 
                                    round(as.numeric(tabfalsepostemp$low[i]),  digits=0), "-",
                                    round(as.numeric(tabfalsepostemp$high[i]), digits=0), ")", sep="")
  }
  falsepostab <- tabfalsepostemp %>%
    relocate(5,4,8) %>%
    select(approach,state,res)
  colnames(falsepostab) <- c("Algorithm", "State", "Median (95% UI)")
  falsepostab %>% 
    gt() %>% 
    tab_header(title = md("Number of individuals with false positive results by screening algorithm")) %>%
    cols_align(align = "center") %>%
    opt_table_font(font = 'Open sans') %>%
    gtsave(filename = paste("./Results_",prev,"/table_falsepos_",covn,"_",int,".png",sep=""), vwidth = 1600, vheight = 600)
  
  # Negative results
  covn <- max(covlist)*10
  int  <- min(durlist)
  
  load(file=paste("./Results_",prev,"/res/fullresults.Rdata",sep=""))
  
  negres <- as.data.frame(rbind(fullresults[covn,int,7,,],fullresults[covn,int,8,,],fullresults[covn,int,9,,],fullresults[covn,int,12,,]))
  colnames(negres) <- c("low","med","high")
  negres$state     <- c(rep("sTB",(algnum-1)),rep("aTB",(algnum-1)),rep("nTB",(algnum-1)),rep("not TB",(algnum-1)))
  negres$approach  <- c(rep(c("sympxpert","cxrxpert","xpert","cxr")))
  
  negres$lowadj <- 0
  negres$lowadj[1:(algnum-1)]                    <- negres$low[1:(algnum-1)]
  negres$lowadj[algnum:((algnum-1)*2)]           <- negres$low[algnum:((algnum-1)*2)]           + negres$med[1:(algnum-1)]
  negres$lowadj[((algnum-1)*2+1):((algnum-1)*3)] <- negres$low[((algnum-1)*2+1):((algnum-1)*3)] + negres$med[algnum:((algnum-1)*2)] + negres$med[1:(algnum-1)]
  negres$lowadj[((algnum-1)*3+1):((algnum-1)*4)] <- negres$low[((algnum-1)*3+1):((algnum-1)*4)]
  
  negres$highadj <- 0
  negres$highadj[1:(algnum-1)]                    <- negres$high[1:(algnum-1)]
  negres$highadj[algnum:((algnum-1)*2)]           <- negres$high[algnum:((algnum-1)*2)]           + negres$med[1:(algnum-1)]
  negres$highadj[((algnum-1)*2+1):((algnum-1)*3)] <- negres$high[((algnum-1)*2+1):((algnum-1)*3)] + negres$med[algnum:((algnum-1)*2)] + negres$med[1:(algnum-1)]
  negres$highadj[((algnum-1)*3+1):((algnum-1)*4)] <- negres$high[((algnum-1)*3+1):((algnum-1)*4)]
  
  plotfalseneg <- negres
  plotfalseneg[c(13:16),c(1:3,6:7)] <- 0
  plotfalseneg$approach <- factor(plotfalseneg$approach,levels=c("sympxpert","cxrxpert","xpert","cxr"),labels=c("Cough+Xpert","CXR+Xpert","Xpert","CXR"))
  plotfalseneg$state <- factor(plotfalseneg$state,levels=c("sTB","aTB","nTB","not TB"))
  
  ymax <- max(plotfalseneg$highadj)*1.05
  falseneg <- ggplot(filter(plotfalseneg,approach != "CXR+Xpert"),aes(x=approach,y=med,fill=approach:state)) +
    geom_bar(position=position_stack(reverse = TRUE),stat="identity") +
    geom_errorbar(aes(ymin=lowadj,ymax=highadj),width=0.5,size=0.5) +
    scale_x_discrete("Algorithm") + 
    scale_y_continuous("Number of individuals",expand=c(0,0)) + 
    ggtitle("False negatives") +
    scale_fill_manual(values=c(colorRampPalette(c("white",colsympxpert))(5)[[5]],colorRampPalette(c("white",colsympxpert))(5)[[3]],colorRampPalette(c("white",colsympxpert))(5)[[2]],colorRampPalette(c("white",colsympxpert))(5)[[4]],
                               colorRampPalette(c("white",colxpert))(5)[[5]]    ,colorRampPalette(c("white",colxpert))(5)[[3]]    ,colorRampPalette(c("white",colxpert))(5)[[2]]    ,colorRampPalette(c("white",colxpert))(5)[[4]]    ,
                               colorRampPalette(c("white",colcxr))(5)[[5]]      ,colorRampPalette(c("white",colcxr))(5)[[3]]      ,colorRampPalette(c("white",colcxr))(5)[[2]]      ,colorRampPalette(c("white",colcxr))(5)[[4]]      ,
                               colorRampPalette(c("white",colcxrxpert))(5)[[5]] ,colorRampPalette(c("white",colcxrxpert))(5)[[3]] ,colorRampPalette(c("white",colcxrxpert))(5)[[2]] ,colorRampPalette(c("white",colcxrxpert))(5)[[4]] )) + 
    theme_classic(base_size=16) +
    labs(fill="Algorithm &\ndisease state") +
    guides(fill=guide_legend(nrow=3,byrow=TRUE)) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename=paste("./Results_",prev,"/plot_falseneg_",covn,"_",int,".tiff",sep=""),plot=falseneg,width=30,height=15, unit="cm")
  
  tabfalsenegtemp <- plotfalseneg[c(1:12),]
  tabfalsenegtemp$state <- factor(tabfalsenegtemp$state,levels=c("sTB","aTB","nTB","not TB","Total"))  
  for (i in levels(plotfalseneg$approach)){
    tempmed  <- 0
    templow  <- 0
    temphigh <- 0
    for (j in 1:nrow(plotfalseneg)){
      if (plotfalseneg$approach[j] == i){
        tempmed  <- tempmed  + plotfalseneg$med[j]
        templow  <- templow  + plotfalseneg$low[j]
        temphigh <- temphigh + plotfalseneg$high[j]
      }
    }
    tabfalsenegtemp <- rbind(tabfalsenegtemp,c(templow,tempmed,temphigh,"Total",i,"",""))
  }
  for (i in 1:nrow(tabfalsenegtemp)){
    tabfalsenegtemp$res[i] <- paste(round(as.numeric(tabfalsenegtemp$med[i]),  digits=0), " (", 
                                    round(as.numeric(tabfalsenegtemp$low[i]),  digits=0), "-",
                                    round(as.numeric(tabfalsenegtemp$high[i]), digits=0), ")", sep="")
  }
  falsenegtab <- tabfalsenegtemp %>%
    filter(state != "not TB") %>%
    arrange(approach,factor(state,c("sTB","aTB","nTB"))) %>%
    relocate(5,4,8) %>%
    select(approach,state,res)
  colnames(falsenegtab) <- c("Algorithm", "State", "Median (95% UI)")
  falsenegtab$Algorithm <- factor(falsenegtab$Algorithm,levels=c("Cough+Xpert","CXR+Xpert","Xpert","CXR",""),labels=c("Cough+Xpert","CXR+Xpert","Xpert","CXR",""))
  falsenegtab[c(2:4,6:8,10:12,14:16),1] <- ""
  falsenegtab %>% 
    gt() %>% 
    tab_header(title = md("Number of individuals with false negative results by screening algorithm")) %>%
    cols_align(align = "center") %>%
    opt_table_font(font = 'Open sans') %>%
    gtsave(filename = paste("./Results_",prev,"/table_falseneg_",covn,"_",int,".png",sep=""), vwidth = 1600, vheight = 600)
  
  plottrueneg <- negres[negres$state == "not TB", ]
  plottrueneg$approach <- factor(plottrueneg$approach,levels=c("sympxpert","cxrxpert","xpert","cxr",""),labels=c("Cough+Xpert","CXR+Xpert","Xpert","CXR",""))
  
  ymax <- max(plottrueneg$highadj)*1.05
  trueneg <- ggplot(filter(plottrueneg,approach != "CXR+Xpert"),aes(x=approach,y=med,fill=approach)) +
    geom_bar(position=position_stack(reverse = TRUE),stat="identity") +
    geom_errorbar(aes(ymin=lowadj,ymax=highadj),width=0.5,size=0.5) +
    scale_x_discrete("Algorithm") + 
    scale_y_continuous("Number of individuals",expand=c(0,0)) + 
    ggtitle("True negatives") +
    scale_fill_manual(values=c(colorRampPalette(c("white",colsympxpert))(5)[[4]],
                               colorRampPalette(c("white",colxpert))(5)[[4]], 
                               colorRampPalette(c("white",colcxr))(5)[[4]],
                               colorRampPalette(c("white",colcxrxpert))(5)[[4]])) + 
    theme_classic(base_size=16) +
    labs(fill="Disease state") +
    guides(fill=guide_legend(nrow=3,byrow=TRUE)) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename=paste("./Results_",prev,"/plot_trueneg_",covn,"_",int,".tiff",sep=""),plot=trueneg,width=30,height=15, unit="cm")
  
  tabtruenegtemp <- plottrueneg
  for (i in 1:nrow(tabtruenegtemp)){
    tabtruenegtemp$res[i] <- paste(round(as.numeric(tabtruenegtemp$med[i]),  digits=0), " (", 
                                   round(as.numeric(tabtruenegtemp$low[i]),  digits=0), "-",
                                   round(as.numeric(tabtruenegtemp$high[i]), digits=0), ")", sep="")
  }
  truenegtab <- tabtruenegtemp %>%
    relocate(5,4,8) %>%
    select(approach,state,res)
  colnames(truenegtab) <- c("Algorithm", "State", "Median (95% UI)")
  truenegtab %>% 
    gt() %>% 
    tab_header(title = md("Number of individuals with true negative results by screening algorithm")) %>%
    cols_align(align = "center") %>%
    opt_table_font(font = 'Open sans') %>%
    gtsave(filename = paste("./Results_",prev,"/table_trueneg_",covn,"_",int,".png",sep=""), vwidth = 1600, vheight = 600)
  
  testres <- ggarrange(truepos, falsepos, falseneg, trueneg, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
  ggsave(filename=paste("./Results_",prev,"/plot_testres_",covn,"_",int,".tiff",sep=""),plot=testres,width=30,height=20, unit="cm", bg="white")
  
}  

# Number of screening and confirmatory tests
for (cov in covlist){
  if (cov == 0.1)  covn <-  1
  if (cov == 0.2)  covn <-  2
  if (cov == 0.3)  covn <-  3
  if (cov == 0.4)  covn <-  4
  if (cov == 0.5)  covn <-  5
  if (cov == 0.6)  covn <-  6
  if (cov == 0.7)  covn <-  7
  if (cov == 0.8)  covn <-  8
  if (cov == 0.9)  covn <-  9
  if (cov == 1.0)  covn <- 10
  
  for (int in durlist){
    
    figtestlist <- list()
    tabtesttemp <- c()
    for (prev in prevlist){
      if (prev == 250)  prevnum  <- 1
      if (prev == 500)  prevnum  <- 2
      if (prev == 1000) prevnum  <- 3
      
      load(file=paste("./Results_",prev,"/res/fullresults.Rdata",sep=""))
      
      plottest <- as.data.frame(rbind(fullresults[covn,int,18,,],fullresults[covn,int,19,,]))
      plottest[4,c(1:3)] <- c(0,0,0)
      colnames(plottest) <- c("low","med","high")
      plottest$test      <- c(rep("Screening",(algnum-1)),rep("Confirmatory",(algnum-1)))
      plottest$approach  <- c(rep(c("sympxpert","cxrxpert","xpert","cxr")))
      plottest$approach  <- factor(plottest$approach,levels=c("sympxpert","cxrxpert","xpert","cxr"),labels=c("Cough+Xpert","CXR+Xpert","Xpert","CXR"))
      
      plottest$lowadj <- 0
      plottest$lowadj[1:(algnum-1)]          <- plottest$low[algnum:((algnum-1)*2)] + plottest$med[1:(algnum-1)]
      plottest$lowadj[algnum:((algnum-1)*2)] <- plottest$low[algnum:((algnum-1)*2)]
      
      plottest$highadj <- 0
      plottest$highadj[1:(algnum-1)]          <- plottest$high[algnum:((algnum-1)*2)] + plottest$med[1:(algnum-1)]
      plottest$highadj[algnum:((algnum-1)*2)] <- plottest$high[algnum:((algnum-1)*2)]
    
      ymax <- max(plottest$highadj)*1.0005
      if (prev == 1000){
        figtest <- ggplot(filter(plottest,approach != "CXR+Xpert"),aes(x=approach,y=med,fill=interaction(approach,test),label=round(med,0))) +
                   geom_bar(position=position_stack(reverse = TRUE),stat="identity") +
                   geom_errorbar(aes(ymin=lowadj,ymax=highadj),width=0.5,size=0.5) +
                   scale_x_discrete(paste("Prevalence ",prev," per 100,000",sep="")) + 
                   scale_y_continuous(str_wrap("Number of tests",width=28),expand=c(0,0),breaks=seq(0,125000,25000),limits=c(0,ymax)) + 
                   scale_fill_manual(values=c(colsympxpert,colxpert,colcxr,colorRampPalette(c("white",colsympxpert))(5)[[2]],colorRampPalette(c("white",colxpert))(5)[[2]],colorRampPalette(c("white",colcxr))(5)[[2]]),
                                     labels=c("Cough+Xpert: Confirmatory (GeneXpert Ultra)","Xpert: Confirmatory (GeneXpert Ultra)","CXR: Confirmatory (CXR)",
                                              "Cough+Xpert: Initial screen (Cough)","Xpert: Screening (none)","CXR: Screening (none)")) +                theme_classic(base_size=16) +
                   theme_classic(base_size=16) +
                   labs(fill="Algorithm &\nTest type") +
                   guides(fill=guide_legend(nrow=2,byrow=TRUE))
      }else{
        figtest <- ggplot(filter(plottest,approach != "CXR+Xpert"),aes(x=approach,y=med,fill=interaction(approach,test),label=round(med,0))) +
                   geom_bar(position=position_stack(reverse = TRUE),stat="identity") +
                   geom_errorbar(aes(ymin=lowadj,ymax=highadj),width=0.5,size=0.5) +
                   scale_x_discrete(paste("Prevalence ",prev," per 100,000",sep="")) + 
                   scale_y_continuous("",expand=c(0,0),breaks=seq(0,125000,25000),limits=c(0,ymax)) + 
                   scale_fill_manual(values=c(colsympxpert,colxpert,colcxr,colorRampPalette(c("white",colsympxpert))(5)[[2]],colorRampPalette(c("white",colxpert))(5)[[2]],colorRampPalette(c("white",colcxr))(5)[[2]]),
                                     labels=c("Cough+Xpert: Confirmatory (GeneXpert Ultra)","Xpert: Confirmatory (GeneXpert Ultra)","CXR: Confirmatory (CXR)",
                                              "Cough+Xpert: Initial screen (Cough)","Xpert: Screening (none)","CXR: Screening (none)")) +                theme_classic(base_size=16) +
                   theme(axis.text.y = element_blank()) +
                   labs(fill="Algorithm &\nTest type") +
                   guides(fill=guide_legend(nrow=2,byrow=TRUE))
      }  
      figtestlist[[prevnum]] <- figtest
      plottest$prev <- prev
      tabtesttemp <- rbind(tabtesttemp, plottest)
    }
    figtests <- ggarrange(figtestlist[[3]], figtestlist[[2]], figtestlist[[1]], ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
    ggsave(filename=paste("./Results/plottests","_",covn,"_",int,".tiff",sep=""),plot=figtests,width=40,height=20, unit="cm", bg="white")

    for (i in 1:nrow(tabtesttemp)){
      tabtesttemp$res[i] <- paste(scales::comma(round(as.numeric(tabtesttemp$med[i]) , digits=0)), " (", 
                                  scales::comma(round(as.numeric(tabtesttemp$low[i]) , digits=0)), "-",
                                  scales::comma(round(as.numeric(tabtesttemp$high[i]), digits=0)), ")", sep="")
    }
    tabtest <- tabtesttemp %>%
      filter(approach != "CXR+Xpert") %>%
      relocate(8,5,4,9) %>%
      select(test,prev,approach,res) %>% 
      pivot_wider(names_from=test,values_from=res) %>%
      arrange(prev)
    colnames(tabtest) <- c("Prevalence", "Algorithm", "Screening tests", "Confirmatory tests")
    tabtest$Prevalence <- as.character(tabtest$Prevalence)
    tabtest[c(2:3,5:6,8:9),1] <- ""
    
    tabtest %>% 
      gt() %>% 
      tab_header(title = md("Number of screening and confirmatory tests")) %>%
      cols_align(align = "center") %>%
      opt_table_font(font = 'Open sans') %>%
      gtsave(filename = paste("./Results/table_tests_",covn,"_",int,".png",sep=""), vwidth = 1600, vheight = 600)
  }
}

# Plot trends in incidence and mortality
for (prev in prevlist){
  
  plotlistinc  <- measplot("inc",covlist,durlist,algnum,plotmos)
  figinctemp   <- (plotlistinc[[ 1]] | plotlistinc[[ 6]] | plotlistinc[[11]] | plotlistinc[[16]] | plotlistinc[[21]] | plotlistinc[[26]] | plotlistinc[[31]] | plotlistinc[[36]] | plotlistinc[[41]] | plotlistinc[[46]]) / 
                  (plotlistinc[[ 2]] | plotlistinc[[ 7]] | plotlistinc[[12]] | plotlistinc[[17]] | plotlistinc[[22]] | plotlistinc[[27]] | plotlistinc[[32]] | plotlistinc[[37]] | plotlistinc[[42]] | plotlistinc[[47]]) / 
                  (plotlistinc[[ 3]] | plotlistinc[[ 8]] | plotlistinc[[13]] | plotlistinc[[18]] | plotlistinc[[23]] | plotlistinc[[28]] | plotlistinc[[33]] | plotlistinc[[38]] | plotlistinc[[43]] | plotlistinc[[48]]) / 
                  (plotlistinc[[ 4]] | plotlistinc[[ 9]] | plotlistinc[[14]] | plotlistinc[[19]] | plotlistinc[[24]] | plotlistinc[[29]] | plotlistinc[[34]] | plotlistinc[[39]] | plotlistinc[[44]] | plotlistinc[[49]]) / 
                  (plotlistinc[[ 5]] | plotlistinc[[10]] | plotlistinc[[15]] | plotlistinc[[20]] | plotlistinc[[25]] | plotlistinc[[30]] | plotlistinc[[35]] | plotlistinc[[40]] | plotlistinc[[45]] | plotlistinc[[50]]) +
                  plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  figinc       <- figinctemp + plot_annotation(title=paste("Baseline prevalence of ",prev," per 100,000",sep=""),
                                               theme=theme(plot.title=element_text(hjust=0.5,size=30)))
  ggsave(filename=paste("./Results_",prev,"/plotinc.tiff",sep=""),plot=figinc,width=90,height=45, unit="cm")
  
  plotlistmort <- measplot("mort",covlist,durlist,algnum,plotmos)
  figmorttemp  <- (plotlistmort[[ 1]] | plotlistmort[[ 6]] | plotlistmort[[11]] | plotlistmort[[16]] | plotlistmort[[21]] | plotlistmort[[26]] | plotlistmort[[31]] | plotlistmort[[36]] | plotlistmort[[41]] | plotlistmort[[46]]) / 
                  (plotlistmort[[ 2]] | plotlistmort[[ 7]] | plotlistmort[[12]] | plotlistmort[[17]] | plotlistmort[[22]] | plotlistmort[[27]] | plotlistmort[[32]] | plotlistmort[[37]] | plotlistmort[[42]] | plotlistmort[[47]]) / 
                  (plotlistmort[[ 3]] | plotlistmort[[ 8]] | plotlistmort[[13]] | plotlistmort[[18]] | plotlistmort[[23]] | plotlistmort[[28]] | plotlistmort[[33]] | plotlistmort[[38]] | plotlistmort[[43]] | plotlistmort[[48]]) / 
                  (plotlistmort[[ 4]] | plotlistmort[[ 9]] | plotlistmort[[14]] | plotlistmort[[19]] | plotlistmort[[24]] | plotlistmort[[29]] | plotlistmort[[34]] | plotlistmort[[39]] | plotlistmort[[44]] | plotlistmort[[49]]) / 
                  (plotlistmort[[ 5]] | plotlistmort[[10]] | plotlistmort[[15]] | plotlistmort[[20]] | plotlistmort[[25]] | plotlistmort[[30]] | plotlistmort[[35]] | plotlistmort[[40]] | plotlistmort[[45]] | plotlistmort[[50]]) +
                  plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  figmort      <- figmorttemp + plot_annotation(title=paste("Baseline prevalence of ",prev," per 100,000",sep=""),
                                                theme=theme(plot.title=element_text(hjust=0.5,size=30)))
  ggsave(filename=paste("./Results_",prev,"/plotmort.tiff",sep=""),plot=figmort,width=90,height=45, unit="cm")
  
}

# Table of % decline in incidence
for (prev in prevlist){
  pctinc  <- as.data.frame(meascalc("inc",covlist,durlist,algnum,plotmos))

  pctinc$"Cough+Xpert" <- paste(round(pctinc$sympxpert_med *100,digits=1), "%\n(", 
                            round(pctinc$sympxpert_low *100,digits=1), "-",
                            round(pctinc$sympxpert_high*100,digits=1), "%)", sep="")
  pctinc$"CXR+Xpert"     <- paste(round(pctinc$cxrxpert_med *100,digits=1), "%\n(", 
                            round(pctinc$cxrxpert_low *100,digits=1), "-",
                            round(pctinc$cxrxpert_high*100,digits=1), "%)", sep="")
  pctinc$"Xpert"    <- paste(round(pctinc$xpert_med *100,digits=1), "%\n(", 
                            round(pctinc$xpert_low *100,digits=1), "-",
                            round(pctinc$xpert_high*100,digits=1), "%)", sep="")
  pctinc$"CXR"      <- paste(round(pctinc$cxr_med *100,digits=1), "%\n(", 
                            round(pctinc$cxr_low *100,digits=1), "-",
                            round(pctinc$cxr_high*100,digits=1), "%)", sep="")

  pctinc_temp <- pctinc %>%
    select(cov,int,"Cough+Xpert","Xpert","CXR") %>%
    pivot_longer(!c(int,cov), names_to = "approach", values_to = "res") %>%
    pivot_wider(names_from=cov,values_from=res) %>%
    arrange(factor(approach,levels=c("Cough+Xpert","Xpert","CXR"))) %>%
    relocate(2)
  colnames(pctinc_temp) <- c("Algorithm","Rounds","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")
  write_xlsx(pctinc_temp,paste("./Results_",prev,"/tableincpct.xlsx",sep=""))
  
  pctinctab <- pctinc_temp %>%
    gt(groupname_col="Algorithm") %>% 
    tab_header(title = md("Projected reduction in 2035 sTB incidence (%)")) %>%
    tab_spanner(label="Population coverage",columns=c(3:12)) %>%
    cols_align(align = "center") %>%
    opt_table_font(font = 'Open sans') %>%
    cols_width(everything() ~ px(100)) %>%
    cols_width(Rounds ~ px(80)) %>%
    tab_style(
      style=cell_text(whitespace="pre"),
      locations=cells_body(columns=c(2:12),rows=everything())) %>%
      gtsave(filename = paste("./Results_",prev,"/tableincpct.html",sep=""))
}

# Table of % decline in mortality
for (prev in prevlist){
  pctmort  <- as.data.frame(meascalc("mort",covlist,durlist,algnum,plotmos))
  
  pctmort$"Cough+Xpert" <- paste(round(pctmort$sympxpert_med *100,digits=1), "%\n(", 
                             round(pctmort$sympxpert_low *100,digits=1), "-",
                             round(pctmort$sympxpert_high*100,digits=1), "%)", sep="")
  pctmort$"CXR+Xpert"     <- paste(round(pctmort$cxrxpert_med *100,digits=1), "%\n(", 
                             round(pctmort$cxrxpert_low *100,digits=1), "-",
                             round(pctmort$cxrxpert_high*100,digits=1), "%)", sep="")
  pctmort$"Xpert"    <- paste(round(pctmort$xpert_med *100,digits=1), "%\n(", 
                             round(pctmort$xpert_low *100,digits=1), "-",
                             round(pctmort$xpert_high*100,digits=1), "%)", sep="")
  pctmort$"CXR"      <- paste(round(pctmort$cxr_med *100,digits=1), "%\n(", 
                             round(pctmort$cxr_low *100,digits=1), "-",
                             round(pctmort$cxr_high*100,digits=1), "%)", sep="")
  
  pctmort_temp <- pctmort %>%
    select(cov,int,"Cough+Xpert","Xpert","CXR") %>%
    pivot_longer(!c(int,cov), names_to = "approach", values_to = "res") %>%
    pivot_wider(names_from=cov,values_from=res) %>%
    arrange(factor(approach,levels=c("Cough+Xpert","Xpert","CXR"))) %>%
    relocate(2)
  colnames(pctmort_temp) <- c("Algorithm","Rounds","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")
  write_xlsx(pctmort_temp,paste("./Results_",prev,"/tablemortpct.xlsx",sep=""))
  
 # pctmort_temp %>%
 #   gt(groupname_col="Algorithm") %>% 
 #   tab_header(title = md("Projected reduction in 2035 TB-associated mortality (%)")) %>%
 #   tab_spanner(label="Population coverage",columns=c(3:12)) %>%
 #   cols_align(align = "center") %>%
 #   opt_table_font(font = 'Open sans') %>%
 #   cols_width(everything() ~ px(100)) %>%
 #   cols_width(Rounds ~ px(80)) %>%
 #   tab_style(
 #     style=cell_text(whitespace="pre"),
 #     locations=cells_body(columns=c(2:12),rows=everything())) %>%
 #   gtsave(filename = paste("./Results_",prev,"/tablemortpct.html",sep=""))
 # webshot2::webshot(url=paste("./Results_",prev,"/tablemortpct.html",sep=""),
 #                   file=paste("./Results_",prev,"/tablemortpct.png",sep=""),
 #                   zoom = 1)
}

# sTB episodes and deaths averted (n)
for (prev in prevlist){ 
  load(file=paste("./Results_",prev,"/res/fullresults.Rdata",sep=""))
  plotlist <- list()
  for (outnum in c(1:2)){
    outdatalist <- c()
    out <- 25 + outnum ^ 2
    
    if (out==26)  measname <- "sTB episodes averted (%)"
    if (out==29)  measname <- "TB-associated deaths averted (%)"
    
    ylimit <- max(fullresults[,,out,,3])
    for (test in c(1:4)){
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
        testname <- "Xpert only"
      }
      if (test==4){
        testcol  <- colcxr
        testname <- "CXR only"
      }
      for (covn in c(1:covnum)){
        for (int in durlist){
          if (covn== 1)  cov <- 0.1
          if (covn== 2)  cov <- 0.2
          if (covn== 3)  cov <- 0.3
          if (covn== 4)  cov <- 0.4
          if (covn== 5)  cov <- 0.5
          if (covn== 6)  cov <- 0.6
          if (covn== 7)  cov <- 0.7
          if (covn== 8)  cov <- 0.8
          if (covn== 9)  cov <- 0.9
          if (covn==10)  cov <- 1.0
          testn   <- test
          temp    <- t(fullresults[covn,int,out,testn,])
          outdata <- rbind(outdata,c(temp,covn,cov,int,testname))
        }
      }
      outdata <- as.data.frame(outdata)
      colnames(outdata) <- c("lower","median","upper","covn","cov","int","approach")
      outdata$cov <- as.character(outdata$cov)
      if (outnum == 1){
        if (test == 1){
          outdataplot <- ggplot(outdata,aes(x=int,y=-as.numeric(median),fill=cov)) +
            geom_bar(position=position_dodge(),stat="identity",colour="black") +
            geom_errorbar(aes(ymin=-as.numeric(upper),ymax=-as.numeric(lower)),width=.2,position=position_dodge(.9)) +
            scale_x_discrete("Rounds and coverage",position = "top") +
            scale_y_continuous(measname,expand=c(0,0),labels=scales::comma,lim=c(-100,0)) +
            scale_fill_manual("Coverage", values=colorRampPalette(c("white",testcol))(10),
                              labels=c("10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")) +
            ggtitle(testname) +
            theme_classic(base_size=16) +
            theme(legend.position="bottom",plot.title=element_text(hjust=0.5))+
            guides(fill=guide_legend(nrow=1,byrow=TRUE))
          ggsave(filename=paste("./Results_",prev,"/plot_incmort_legend.tiff",sep=""),plot=outdataplot,width=30,height=20, unit="cm")
          
        }else{
          outdataplot <- ggplot(outdata,aes(x=int,y=-as.numeric(median),fill=cov)) +
            geom_bar(position=position_dodge(),stat="identity",colour="black") +
            geom_errorbar(aes(ymin=-as.numeric(upper),ymax=-as.numeric(lower)),width=.2,position=position_dodge(.9)) +
            scale_x_discrete("Rounds and coverage",position = "top") +
            scale_y_continuous("",expand=c(0,0),labels=scales::comma,lim=c(-100,0)) +
            scale_fill_manual("Coverage", values=colorRampPalette(c("white",testcol))(10)) +
            ggtitle(testname) +
            theme_classic(base_size=16) +
            theme(legend.position="",plot.title=element_text(hjust=0.5),axis.text.y=element_blank())
        }
      }
      if (outnum == 2){
        if (test == 1){
          outdataplot <- ggplot(outdata,aes(x=int,y=-as.numeric(median),fill=cov)) +
            geom_bar(position=position_dodge(),stat="identity",colour="black") +
            geom_errorbar(aes(ymin=-as.numeric(upper),ymax=-as.numeric(lower)),width=.2,position=position_dodge(.9)) +
            scale_x_discrete("",position = "top") +
            scale_y_continuous(measname,expand=c(0,0),labels=scales::comma,lim=c(-100,0)) +
            scale_fill_manual("Coverage", values=colorRampPalette(c("white",testcol))(10)) +
            ggtitle("") +
            theme_classic(base_size=16) +
            theme(legend.position="",plot.title=element_text(hjust=0.5))
        }else{
          outdataplot <- ggplot(outdata,aes(x=int,y=-as.numeric(median),fill=cov)) +
            geom_bar(position=position_dodge(),stat="identity",colour="black") +
            geom_errorbar(aes(ymin=-as.numeric(upper),ymax=-as.numeric(lower)),width=.2,position=position_dodge(.9)) +
            scale_x_discrete("",position = "top") +
            scale_y_continuous("",expand=c(0,0),labels=scales::comma,lim=c(-100,0)) +
            scale_fill_manual("Coverage", values=colorRampPalette(c("white",testcol))(10)) +
            ggtitle("") +
            theme_classic(base_size=16) +
            theme(legend.position="",plot.title=element_text(hjust=0.5),axis.text.y=element_blank())
        }
      }
      plotlist[[(algnum-1)*(outnum-1)+test]] <- outdataplot
      outdatalist <- rbind(outdatalist,outdata)
    }
    if (outnum == 1){
      fig <- ggarrange(plotlist[[1]],plotlist[[2]],plotlist[[4]],ncol=3) 
      ggsave(filename=paste("./Results_",prev,"/plot_stbeps_n.tiff",sep=""),plot=fig,width=60,height=20, unit="cm")
    }
    if (outnum == 2){
      fig <- ggarrange(plotlist[[5]],plotlist[[6]],plotlist[[8]],ncol=3) 
      ggsave(filename=paste("./Results_",prev,"/plot_tbmort_n.tiff",sep=""),plot=fig,width=60,height=20, unit="cm")
      figall <- ggarrange(plotlist[[1]],plotlist[[2]],plotlist[[4]], 
                          plotlist[[5]],plotlist[[6]],plotlist[[8]], 
                          ncol=3,nrow=2)
      ggsave(filename=paste("./Results_",prev,"/plot_stbeps+tbmort_n.tiff",sep=""),plot=figall,width=60,height=30, unit="cm")
    }
    outdatalist$res <- paste(round(as.numeric(outdatalist$median),digits=1), "%\n(", 
                             round(as.numeric(outdatalist$lower ),digits=1), "-",
                             round(as.numeric(outdatalist$upper ),digits=1), "%)", sep="")
    #colnames(outdatalist) <- c("lower","median","upper","covn","cov","int","approach","res")
    outdata_temp <- outdatalist %>%
      select(cov,int,approach,res) %>%
      pivot_wider(names_from=cov,values_from=res) %>%
      #arrange(factor(approach,levels=c("Cough+Xpert","CXR+Xpert","Xpert","CXR"))) %>%
      relocate(2) %>%
      filter(approach != "CXR+Xpert")
    colnames(outdata_temp) <- c("Algorithm","Rounds","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")
    if (outnum==1){
      write_xlsx(outdata_temp,paste("./Results_",prev,"/table_stbeps_n.xlsx",sep=""))
      outdata_temp %>%
        filter(Algorithm != "CXR+Xpert") %>% 
        gt(groupname_col="Algorithm") %>% 
        tab_header(title = md("Infections averted, 2025-2035 (%)")) %>%
        tab_spanner(label="Population coverage",columns=c(3:12)) %>%
        cols_align(align = "center") %>%
        opt_table_font(font = 'Open sans') %>%
        cols_width(everything() ~ px(100)) %>%
        cols_width(Rounds ~ px(80)) %>%
        tab_style(
          style=cell_text(whitespace="pre"),
          locations=cells_body(columns=c(2:12),rows=everything())) %>%
        gtsave(filename = paste("./Results_",prev,"/table_stbeps_n.html",sep=""))
      webshot2::webshot(url=paste("./Results_",prev,"/table_stbeps_n.html",sep=""),
                        file=paste("./Results_",prev,"/table_stbeps_n.png",sep=""),
                        zoom = 0.9)
    }
    if (outnum==2){
      write_xlsx(outdata_temp,paste("./Results_",prev,"/table_deaths_n.xlsx",sep=""))
      outdata_temp %>%
        filter(Algorithm != "CXR+Xpert") %>% 
        gt(groupname_col="Algorithm") %>% 
        tab_header(title = md("TB-associated deaths averted, 2025-2035 (%)")) %>%
        tab_spanner(label="Population coverage",columns=c(3:12)) %>%
        cols_align(align = "center") %>%
        opt_table_font(font = 'Open sans') %>%
        cols_width(everything() ~ px(100)) %>%
        cols_width(Rounds ~ px(80)) %>%
        tab_style(
          style=cell_text(whitespace="pre"),
          locations=cells_body(columns=c(2:12),rows=everything())) %>%
        gtsave(filename = paste("./Results_",prev,"/table_deaths_n.html",sep=""))
      webshot2::webshot(url=paste("./Results_",prev,"/table_deaths_n.html",sep=""),
                        file=paste("./Results_",prev,"/table_deaths_n.png",sep=""),
                        zoom=1)
    }
  }
}

# sTB episodes and deaths averted
for (prev in prevlist){ 
  load(file=paste("./Results_",prev,"/res/fullresults.Rdata",sep=""))
  plotlist <- list()
  for (outnum in c(1:2)){
    outdatalist <- c()
    out <- 25 + outnum ^ 2
    
    if (out==26)  measname <- "sTB episodes averted (%)"
    if (out==29)  measname <- "TB-associated deaths averted (%)"
    
    ylimit <- max(fullresults[,,out,,3])
    for (test in c(1:4)){
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
        testname <- "Xpert only"
      }
      if (test==4){
        testcol  <- colcxr
        testname <- "CXR only"
      }
      for (covn in c(1:covnum)){
        for (int in durlist){
          if (covn== 1)  cov <- 0.1
          if (covn== 2)  cov <- 0.2
          if (covn== 3)  cov <- 0.3
          if (covn== 4)  cov <- 0.4
          if (covn== 5)  cov <- 0.5
          if (covn== 6)  cov <- 0.6
          if (covn== 7)  cov <- 0.7
          if (covn== 8)  cov <- 0.8
          if (covn== 9)  cov <- 0.9
          if (covn==10)  cov <- 1.0
          testn   <- test
          temp    <- t(fullresults[covn,int,out,testn,])
          outdata <- rbind(outdata,c(temp,covn,cov,int,testname))
        }
      }
      outdata <- as.data.frame(outdata)
      colnames(outdata) <- c("lower","median","upper","covn","cov","int","approach")
      outdata$cov <- as.character(outdata$cov)
      if (outnum == 1){
        if (test == 1){
          outdataplot <- ggplot(outdata,aes(x=int,y=-as.numeric(median),fill=cov)) +
            geom_bar(position=position_dodge(),stat="identity",colour="black") +
            geom_errorbar(aes(ymin=-as.numeric(upper),ymax=-as.numeric(lower)),width=.2,position=position_dodge(.9)) +
            scale_x_discrete("Rounds and coverage",position = "top") +
            scale_y_continuous(measname,expand=c(0,0),labels=scales::comma,lim=c(-100,0)) +
            scale_fill_manual("Coverage", values=colorRampPalette(c("white",testcol))(10),
                              labels=c("10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")) +
            ggtitle(testname) +
            theme_classic(base_size=16) +
            theme(legend.position="bottom",plot.title=element_text(hjust=0.5))+
            guides(fill=guide_legend(nrow=1,byrow=TRUE))
          ggsave(filename=paste("./Results_",prev,"/plot_incmort_legend.tiff",sep=""),plot=outdataplot,width=30,height=20, unit="cm")
          
        }else{
          outdataplot <- ggplot(outdata,aes(x=int,y=-as.numeric(median),fill=cov)) +
            geom_bar(position=position_dodge(),stat="identity",colour="black") +
            geom_errorbar(aes(ymin=-as.numeric(upper),ymax=-as.numeric(lower)),width=.2,position=position_dodge(.9)) +
            scale_x_discrete("Rounds and coverage",position = "top") +
            scale_y_continuous("",expand=c(0,0),labels=scales::comma,lim=c(-100,0)) +
            scale_fill_manual("Coverage", values=colorRampPalette(c("white",testcol))(10)) +
            ggtitle(testname) +
            theme_classic(base_size=16) +
            theme(legend.position="",plot.title=element_text(hjust=0.5),axis.text.y=element_blank())
        }
      }
      if (outnum == 2){
        if (test == 1){
          outdataplot <- ggplot(outdata,aes(x=int,y=-as.numeric(median),fill=cov)) +
            geom_bar(position=position_dodge(),stat="identity",colour="black") +
            geom_errorbar(aes(ymin=-as.numeric(upper),ymax=-as.numeric(lower)),width=.2,position=position_dodge(.9)) +
            scale_x_discrete("",position = "top") +
            scale_y_continuous(measname,expand=c(0,0),labels=scales::comma,lim=c(-100,0)) +
            scale_fill_manual("Coverage", values=colorRampPalette(c("white",testcol))(10)) +
            ggtitle("") +
            theme_classic(base_size=16) +
            theme(legend.position="",plot.title=element_text(hjust=0.5))
        }else{
          outdataplot <- ggplot(outdata,aes(x=int,y=-as.numeric(median),fill=cov)) +
            geom_bar(position=position_dodge(),stat="identity",colour="black") +
            geom_errorbar(aes(ymin=-as.numeric(upper),ymax=-as.numeric(lower)),width=.2,position=position_dodge(.9)) +
            scale_x_discrete("",position = "top") +
            scale_y_continuous("",expand=c(0,0),labels=scales::comma,lim=c(-100,0)) +
            scale_fill_manual("Coverage", values=colorRampPalette(c("white",testcol))(10)) +
            ggtitle("") +
            theme_classic(base_size=16) +
            theme(legend.position="",plot.title=element_text(hjust=0.5),axis.text.y=element_blank())
        }
      }
      plotlist[[(algnum-1)*(outnum-1)+test]] <- outdataplot
      outdatalist <- rbind(outdatalist,outdata)
    }
    if (outnum == 1){
      fig <- ggarrange(plotlist[[1]],plotlist[[2]],plotlist[[4]],ncol=3) 
      ggsave(filename=paste("./Results_",prev,"/plot_stbeps.tiff",sep=""),plot=fig,width=60,height=20, unit="cm")
    }
    if (outnum == 2){
      fig <- ggarrange(plotlist[[5]],plotlist[[6]],plotlist[[8]],ncol=3) 
      ggsave(filename=paste("./Results_",prev,"/plot_tbmort.tiff",sep=""),plot=fig,width=60,height=20, unit="cm")
      figall <- ggarrange(plotlist[[1]],plotlist[[2]],plotlist[[4]], 
                          plotlist[[5]],plotlist[[6]],plotlist[[8]], 
                          ncol=3,nrow=2)
      ggsave(filename=paste("./Results_",prev,"/plot_stbeps+tbmort.tiff",sep=""),plot=figall,width=60,height=30, unit="cm")
    }
    outdatalist$res <- paste(round(as.numeric(outdatalist$median),digits=1), "%\n(", 
                             round(as.numeric(outdatalist$lower ),digits=1), "-",
                             round(as.numeric(outdatalist$upper ),digits=1), "%)", sep="")
    #colnames(outdatalist) <- c("lower","median","upper","covn","cov","int","approach","res")
    outdata_temp <- outdatalist %>%
      select(cov,int,approach,res) %>%
      pivot_wider(names_from=cov,values_from=res) %>%
      #arrange(factor(approach,levels=c("Cough+Xpert","CXR+Xpert","Xpert","CXR"))) %>%
      relocate(2) %>%
      filter(approach != "CXR+Xpert")
    colnames(outdata_temp) <- c("Algorithm","Rounds","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")
    if (outnum==1){
      write_xlsx(outdata_temp,paste("./Results_",prev,"/table_stbeps.xlsx",sep=""))
      outdata_temp %>%
        filter(Algorithm != "CXR+Xpert") %>% 
        gt(groupname_col="Algorithm") %>% 
        tab_header(title = md("Infections averted, 2025-2035 (%)")) %>%
        tab_spanner(label="Population coverage",columns=c(3:12)) %>%
        cols_align(align = "center") %>%
        opt_table_font(font = 'Open sans') %>%
        cols_width(everything() ~ px(100)) %>%
        cols_width(Rounds ~ px(80)) %>%
        tab_style(
          style=cell_text(whitespace="pre"),
          locations=cells_body(columns=c(2:12),rows=everything())) %>%
        gtsave(filename = paste("./Results_",prev,"/table_stbeps.html",sep=""))
      webshot2::webshot(url=paste("./Results_",prev,"/table_stbeps.html",sep=""),
                        file=paste("./Results_",prev,"/table_stbeps.png",sep=""),
                        zoom = 0.9)
    }
    if (outnum==2){
      write_xlsx(outdata_temp,paste("./Results_",prev,"/table_deaths.xlsx",sep=""))
      outdata_temp %>%
        filter(Algorithm != "CXR+Xpert") %>% 
        gt(groupname_col="Algorithm") %>% 
        tab_header(title = md("TB-associated deaths averted, 2025-2035 (%)")) %>%
        tab_spanner(label="Population coverage",columns=c(3:12)) %>%
        cols_align(align = "center") %>%
        opt_table_font(font = 'Open sans') %>%
        cols_width(everything() ~ px(100)) %>%
        cols_width(Rounds ~ px(80)) %>%
        tab_style(
          style=cell_text(whitespace="pre"),
          locations=cells_body(columns=c(2:12),rows=everything())) %>%
        gtsave(filename = paste("./Results_",prev,"/table_deaths.html",sep=""))
      webshot2::webshot(url=paste("./Results_",prev,"/table_deaths.html",sep=""),
                        file=paste("./Results_",prev,"/table_deaths.png",sep=""),
                        zoom=1)
    }
  }
}

# Table of maximum impact on incidence and mortality
for (prev in prevlist){
  
  imptemp <- c()
  for (cov in covlist){
    for (int in durlist){
      load(file=paste("./Results_",prev,"/res/fulldatainc_",cov,"_",int,".Rdata",sep=""))
      load(file=paste("./Results_",prev,"/res/fulldatamort_",cov,"_",int,".Rdata",sep=""))
      
      incsum  <- array(NA,dim=c(algnum,plotmos,3))
      mortsum <- array(NA,dim=c(algnum,plotmos,3))
      for (i in c(1:algnum)){ 
        for (j in c(1:plotmos)){
          incsum[i,j,1:3] <- quantile((fulldatainc[,1,j,1]-fulldatainc[,i,j,1])/fulldatainc[,1,j,1],probs=c(0.5,0.025,0.975))
          mortsum[i,j,1:3] <- quantile((fulldatamort[,1,j,1]-fulldatamort[,i,j,1])/fulldatamort[,1,j,1],probs=c(0.5,0.025,0.975))
        }
        imptemp <- rbind(imptemp,c(cov,int,i,max(incsum[i,,1]),max(incsum[i,,2]),max(incsum[i,,3]),max(mortsum[i,,1]),max(mortsum[i,,2]),max(mortsum[i,,3])))
      }
    }
  }
  imp <- as.data.frame(imptemp)
  colnames(imp) <- c("cov","int","alg","inc_med","inc_low","inc_high","mort_med","mort_low","mort_high")
  for (j in 1:nrow(imp)){
    if (imp$alg[j] == 1) imp$alg[j] <- "None"
    if (imp$alg[j] == 2) imp$alg[j] <- "Cough+Xpert"
    if (imp$alg[j] == 3) imp$alg[j] <- "CXR+Xpert"
    if (imp$alg[j] == 4) imp$alg[j] <- "Xpert"
    if (imp$alg[j] == 5) imp$alg[j] <- "CXR"
    imp$inc[j]  <- paste(round(imp$inc_med[j]*100 ,digits=1), "%\n(", 
                         round(imp$inc_low[j]*100 ,digits=1), "-",
                         round(imp$inc_high[j]*100 ,digits=1), "%)", sep="")
    imp$mort[j] <- paste(round(imp$mort_med[j] *100,digits=1), "%\n(", 
                         round(imp$mort_low[j] *100,digits=1), "-",
                         round(imp$mort_high[j]*100,digits=1), "%)", sep="")
  }
  
  impactinc <- imp %>%
    filter(!alg %in% c("None","CXR+Xpert")) %>%
    select(1:3,10) %>%
    pivot_wider(names_from=cov,values_from=inc) %>%
    relocate(2) %>%
    arrange(factor(alg,levels=c("Cough+Xpert","Xpert","CXR")))
  colnames(impactinc) <- c("Algorithm","Rounds","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")
  write_xlsx(impactinc,paste("./Results_",prev,"/tablemaxinc.xlsx",sep=""))
  
  impactmort <- imp %>%
    filter(!alg %in% c("None","CXR+Xpert")) %>%
    select(1:3,11) %>%
    pivot_wider(names_from=cov,values_from=mort) %>%
    relocate(2) %>%
    arrange(factor(alg,levels=c("Cough+Xpert","Xpert","CXR")))
  colnames(impactmort) <- c("Algorithm","Rounds","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")
  write_xlsx(impactmort,paste("./Results_",prev,"/tablemaxmort.xlsx",sep=""))
  
}

# Table of rebound for incidence and mortality
for (prev in prevlist){
  
  rebtemp <- c()
  for (cov in covlist){
    for (int in durlist){
      load(file=paste("./Results_",prev,"/res/fulldatainc_",cov,"_",int,".Rdata",sep=""))
      load(file=paste("./Results_",prev,"/res/fulldatamort_",cov,"_",int,".Rdata",sep=""))
      
      incsum  <- array(NA,dim=c(algnum,1000,3))
      mortsum <- array(NA,dim=c(algnum,1000,3))
      
      for (i in c(1:algnum)){ 
        for (m in 1:1000){
          incsum[i,m,1]  <- max((fulldatainc[m,1, ,1]-fulldatainc[m,i, ,1])/fulldatainc[m,1, ,1])
          mortsum[i,m,1] <- max((fulldatamort[m,1, ,1]-fulldatamort[m,i, ,1])/fulldatamort[m,1, ,1])

          incsum[i,m,2]  <- ((sum(fulldatainc[m,1,c(121:132),1])/12) -(sum(fulldatainc[m,i,c(121:132),1])/12)) /(sum(fulldatainc[m,1,c(121:132),1])/12)
          mortsum[i,m,2] <- ((sum(fulldatamort[m,1,c(121:132),1])/12)-(sum(fulldatamort[m,i,c(121:132),1])/12))/(sum(fulldatamort[m,1,c(121:132),1])/12)
          
          incsum[i,m,3]  <- (incsum[i,m,1]  - incsum[i,m,2] )/incsum[i,m,1]  
          mortsum[i,m,3] <- (mortsum[i,m,1] - mortsum[i,m,2])/mortsum[i,m,1] 
        }
        if (i != 1)   rebtemp  <- rbind(rebtemp,c(prev,cov,int,i,quantile(incsum[i,,3],probs=c(0.5,0.025,0.975)),quantile(mortsum[i,,3],probs=c(0.5,0.025,0.975))))
      }
    }
  }
  reb <- as.data.frame(rebtemp)
  colnames(reb) <- c("prev","cov","int","alg","inc_med","inc_low","inc_high","mort_med","mort_low","mort_high")
  for (j in 1:nrow(reb)){
    if (reb$alg[j] == 1) reb$alg[j] <- "None"
    if (reb$alg[j] == 2) reb$alg[j] <- "Cough+Xpert"
    if (reb$alg[j] == 3) reb$alg[j] <- "CXR+Xpert"
    if (reb$alg[j] == 4) reb$alg[j] <- "Xpert"
    if (reb$alg[j] == 5) reb$alg[j] <- "CXR"
    reb$inc[j]  <- paste(round(reb$inc_med[j]*100 ,digits=1), "%\n(", 
                         round(reb$inc_low[j]*100 ,digits=1), "-",
                         round(reb$inc_high[j]*100 ,digits=1), "%)", sep="")
    reb$mort[j] <- paste(round(reb$mort_med[j] *100,digits=1), "%\n(", 
                         round(reb$mort_low[j] *100,digits=1), "-",
                         round(reb$mort_high[j]*100,digits=1), "%)", sep="")
  }
  
  reboundinc <- reb %>%
    filter(!alg %in% c("None","CXR+Xpert")) %>%
    select(2:4,11) %>%
    pivot_wider(names_from=cov,values_from=inc) %>%
    relocate(2) %>%
    arrange(factor(alg,levels=c("Cough+Xpert","Xpert","CXR")))
  colnames(reboundinc) <- c("Algorithm","Rounds","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")
  write_xlsx(reboundinc,paste("./Results_",prev,"/tablerebinc.xlsx",sep=""))
  
  reboundmort <- reb %>%
    filter(!alg %in% c("None","CXR+Xpert")) %>%
    select(2:4,12) %>%
    pivot_wider(names_from=cov,values_from=mort) %>%
    relocate(2) %>%
    arrange(factor(alg,levels=c("Cough+Xpert","Xpert","CXR")))
  colnames(reboundmort) <- c("Algorithm","Rounds","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")
  write_xlsx(reboundmort,paste("./Results_",prev,"/tablerebmort.xlsx",sep=""))
  
}

# Table of calibration epidemiology
estims <- read.csv("./Results/tbestims_Screening.csv")

estims$res <- paste(round(estims$Median,digits=1), "\n(", 
                    round(estims$Lower ,digits=1), "-",
                    round(estims$Upper ,digits=1), ")", sep="")
estims$res[c(4,14,24)] <- paste(round(estims$Median[c(4,14,24)],digits=2), "\n(", 
                                round(estims$Lower [c(4,14,24)],digits=2), "-",
                                round(estims$Upper [c(4,14,24)],digits=2), ")", sep="")


estims_temp <- estims %>%
  select(Estimate,Prev,res) %>%
  pivot_wider(names_from=Prev,values_from=res)
estims_temp <- estims_temp[-c(6,7,9,10),]
estims_temp$Estimate <- c("infTB prevalence\n(per 100,000)","TB-associated mortality\n(per 100,000)",
                          "TB notifications\n(per 100,000)","Proportion aTB\namong infTB",
                          "nTB:infTB ratio","sTB incidence\n(per 100,000)")
estims_temp$Estimate <- factor(estims_temp$Estimate,levels=c("infTB prevalence\n(per 100,000)","sTB incidence\n(per 100,000)",
                         "TB-associated mortality\n(per 100,000)","TB notifications\n(per 100,000)",
                         "Proportion aTB\namong infTB","nTB:infTB ratio"))
colnames(estims_temp) <- c("Estimate","250","500","1,000")
estims_temp %>%
  arrange(Estimate) %>%
  gt() %>%
  tab_header(title = md("Epidemiological estimates from calibrated models")) %>%
  tab_spanner(label="Target prevalence per 100,000",columns=c(2:4)) %>%
  cols_align(align = "center") %>%
  opt_table_font(font = 'Open sans') %>%
  cols_width(everything() ~ px(200)) %>%
  cols_width(Estimate ~ px(300)) %>%
  #tab_style(
  #  style=cell_text(whitespace="pre"),
  #  locations=cells_body(columns=c(2:12),rows=everything())) %>%
  gtsave(filename = paste("./Results/epiestims.html",sep=""))
webshot2::webshot(url=paste("./Results/epiestims.html",sep=""),
                  file=paste("./Results/epiestims.png",sep=""),
                  zoom = 1)
