#------------------------------------------------------------------------------#
# Calibration_2_Expand.R                                                      #
# Expand calibration to increase parameter sets                                #
# Last updated 2025-01-17 by KCH                                               #
#------------------------------------------------------------------------------#

# Clear environment
rm(list = ls())

# Load required libraries
library(hmer)
library(deSolve)
library(hexbin)
library(bayesplot)
library(beepr)

# Set working directory (to this file's location)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# Source model and hmer functions
source("./Calibration_0_Functions.R")

# Set model run parameters
prev    <- 500     # infTB prevalence
pop     <- 100000  # Population size
timerun <- 5000    # Run length
k <- 18

# Set calibration targets
targets <- list(
  P1 = c(prev*0.9, prev*1.1)
)

### Load previous output
#load(file=paste("./Output_",prev,"/ranges_",k,".Rdata",sep=''))
load(file=paste("./Output_",prev,"/non_imp_list_",k,".Rdata",sep=''))
#load(file=paste("./Output_",prev,"/non_imp_pts_",k,".Rdata",sep=''))
#load(file=paste("./Output_",prev,"/ems_",k,".Rdata",sep=''))

# If 10% above then this next line should give ~1000 calibrated parameter sets; if it's <1000 then just run more
param_options <- generate_new_design(non_imp_list, 1100, targets, verbose = TRUE)  
save(param_options,file=paste("./Output_",prev,"/param_options.Rdata",sep=''))
#load(file=paste("./Output_",prev,"/param_options.Rdata",sep=''))
# Generate final calibration points
results_f <- calibration_function(param_options,prev,pop,timerun)
# Pair up parameter sets with calibration points
save(results_f,file=paste("./Output_",prev,"/results_f.Rdata",sep=''))
#load(file=paste("./Output_",prev,"/results_f.Rdata",sep=''))

param_options_tot <-cbind(param_options, results_f)

# Select calibrated runs (can match to criteria above - could be 100 and 500?)
param_options_tot <- param_options_tot[which(param_options_tot$P1 > prev*0.9),]
param_options_tot <- param_options_tot[which(param_options_tot$P1 < prev*1.1),]

save(param_options_tot, file=paste("./Output_",prev,"/param_options_tot.Rdata",sep=''))
#load(file=paste("./Output_",prev,"/param_options_tot.Rdata",sep=''))

# Subset to 1000 random parameter sets
param_options_1000 <- param_options_tot[sample(nrow(param_options_tot), size=1000), ]
save(param_options_1000, file=paste("./Output_",prev,"/param_options_1000.Rdata",sep=''))

# Multiply to annual
param_options_1000[, 5] = param_options_1000[, 5]*12
param_options_1000[, 6] = param_options_1000[, 6]*12
param_options_1000[, 7] = param_options_1000[, 7]*12
param_options_1000[, 8] = param_options_1000[, 8]*12
param_options_1000[, 9] = param_options_1000[, 9]*12
param_options_1000[,10] = param_options_1000[,10]*12
param_options_1000[,11] = param_options_1000[,11]*12
param_options_1000[,12] = param_options_1000[,12]*12
param_options_1000[,13] = param_options_1000[,13]*12
param_options_1000[,14] = param_options_1000[,14]*12

# Write params (annual) to a csv which then gets used by intervention files
postparams <- as.data.frame(t(apply(param_options_1000,2,quantile,probs=c(0.025,0.50,0.975))))
postparams$param <- row.names(postparams)  
write.csv(postparams, file=paste("./Output_",prev,"/postparams.csv",sep=''), row.names = FALSE)

