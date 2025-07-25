#------------------------------------------------------------------------------#
# Intervention_1_RunModel.R                                                    #
# Main script to run interventions for calibrated model                        #
# Last updated 2025-01-17 by KCH                                               #
#------------------------------------------------------------------------------#

# Clear environment
rm(list = ls())

# Load required libraries
library(deSolve)
library(ggplot2)
library(tictoc)

# Set working directory (to this file's location)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# Set scenario parameters
pop      <- 100000                                      # Population
datapts  <- 51
prevlist <- c(1000,500,250)                             # List of prevalences
covlist  <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)  # List of coverages
durlist  <- c(1,2,3,4,5)                                # List of durations

# Set up test performance data
testdata <- read.csv(file="./testperformance.csv")
no_screen     <- testdata[testdata$id == 'none',]
sputum_screen <- testdata[testdata$id == 'sputum',]
cxr_screen    <- testdata[testdata$id == 'cxr',]
symp_screen   <- testdata[testdata$id == 'symptom',]
xpert_screen  <- testdata[testdata$id == 'xpert',]
xpert_confirm <- testdata[testdata$id == 'xpert_confirm',]

for (prev in prevlist){
  print(prev)

  # Source functions and files
  source("./Intervention_0_Functions.R")
  load(file=paste("./Output_",prev,"/param_options_1000.Rdata",sep=""))
  
  # Following line is used to rename posterior parameters from model run prior to renaming of disease states in Oct 2024
  colnames(param_options_1000) <- c('beta','t','p','r','infclr','ntbrec','infntb','infatb','ntbatb','atbntb','atbstb','stbatb','stbtrt','mustb','P1')
  
  # Run interventions across coverages and rounds
  for (cov in covlist){ 
    print(cov)
    for (int in durlist){ 
      print(int)
      tic()
      
      none         <- list()
      sympxpert    <- list()
      cxrxpert     <- list()
      xpert        <- list()
      cxr          <- list()
      run_sets_rec <- list()
      run_sets     <- list()
      
      for(i in 1:nrow(param_options_1000)){
        print(i)
        
        cov_param  <- c(popcov = cov)
        inttimes   <- int
        parameters <- c(param_options_1000[i,c(1:14)], mu = 0.0137/12, trtrec = 1/6)
        names(parameters)[1:14] <- colnames(param_options_1000)[1:14]
          
        none[[i]]      <- run_base_set(parameters,screen_parameters=NA,confirm_parameters=NA,cov_param,pop,prev,intyears=inttimes)
        baseres        <- as.data.frame(none[[i]])[1:4670,1:52]
        sympxpert[[i]] <- run_int_set(baseres,parameters,uncert_screen(symp_screen)  ,uncert_confirm(xpert_confirm),cov_param,pop,prev,inttimes)
        cxrxpert[[i]]  <- run_int_set(baseres,parameters,uncert_screen(cxr_screen)   ,uncert_confirm(xpert_confirm),cov_param,pop,prev,inttimes)
        xpert[[i]]     <- run_int_set(baseres,parameters,uncert_screen(sputum_screen),uncert_confirm(xpert_screen ),cov_param,pop,prev,inttimes)
        cxr[[i]]       <- run_int_set(baseres,parameters,uncert_screen(no_screen)    ,uncert_confirm(cxr_screen)   ,cov_param,pop,prev,inttimes)
        run_sets[[i]]  <- list(as.list(as.data.frame(     none[[i]])[4670:5001,1:datapts+1]),
                               as.list(as.data.frame(sympxpert[[i]])[4670:5001,1:datapts+1]),
                               as.list(as.data.frame( cxrxpert[[i]])[4670:5001,1:datapts+1]), 
                               as.list(as.data.frame(    xpert[[i]])[4670:5001,1:datapts+1]),
                               as.list(as.data.frame(      cxr[[i]])[4670:5001,1:datapts+1]))
      }
      
      save(run_sets,file=paste("./Output_",prev,"/run_sets_",int,"_",cov,".Rdata",sep=""))
      toc()
    }
  }
}