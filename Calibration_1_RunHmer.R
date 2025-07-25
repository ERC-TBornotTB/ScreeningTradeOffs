#------------------------------------------------------------------------------#
# Calibration_1_RunHmer.R                                                      #
# Run calibration using hmer                                                   #
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

# Create folder to store results
dir.create(file.path(paste("./Output_",prev,sep="")), showWarnings=FALSE)

# Create empty lists to store data from the different waves
ems         <- list() # ems[[k]] will contain wave-k emulators
wave_data   <- list() # wave_data[[k]] will contain the data used to train and validate wave-k emulators
ranges      <- list() # ranges[[k]] will contain the parameter ranges used to train wave-k emulators
non_imp_pts <- list() # non_imp_pts[[k]] will contain the non-implausible points generated at the end of wave k

# Import and format parameters
params_round <- as.data.frame(read.csv("./params_round.csv")[,c(2,4,7)])
postparams   <- params_round[,-1]
rownames(postparams) <- params_round[,1]
postparams   <- as.data.frame(t(postparams))

ranges[[1]] <- list(beta = c(0,12),                                                    # Transmission parameter (1000)
                    t    = c(0.62,1.00),                                               # Relative transmission from aTB
                                                                
                    p    = c(0.14,0.30),                                               # Protection for recovered
                    r    = c(2.14,4.27),                                               # Risk for treated
                       
                    infclr  = c(postparams$infclr[1]/12,postparams$infclr[2]/12),      # Infection to susceptible
                    ntbrec  = c(postparams$minrec[1]/12,postparams$minrec[2]/12),      # nTB to recovery
                         
                    infntb  = c(postparams$infmin[1]/12,postparams$infmin[2]/12),      # Infection to nTB
                    infatb  = c(postparams$infsub[1]/12,postparams$infsub[2]/12),      # Infection to aTB
                       
                    ntbatb  = c(postparams$min_sub[1]/12,postparams$min_sub[2]/12),    # nTB to aTB
                    atbntb  = c(postparams$sub_min[1]/12,postparams$sub_min[2]/12),    # aTB to nTB
                    
                    atbstb  = c(postparams$sub_clin[1]/12,postparams$sub_clin[2]/12),  # aTB to sTB
                    stbatb  = c(postparams$clin_sub[1]/12,postparams$clin_sub[2]/12),  # sTB to aTB
                     
                    stbtrt  = c(0.57/12, 0.77/12),                                     # sTB to treatment

                    mus     = c(postparams$mort[1]/12,postparams$mort[2]/12))          # TB mortality
                                                                                       # Background mortality and treatment duration are defined 
                                                                                       #  in line 55 of Calibration_0_Functions.R

# Set calibration targets
targets <- list(
  P1 = c(prev*0.9,prev*1.1)
)

# Define Latin hypercube design
k <- 1

############################################## Define a latin hypercube design ##############################################

# This can be done through the function `randomLHS`, which assumes that each parameter is distributed on [0,1]
initial_LHS <- lhs::randomLHS(20 * length(ranges[[1]]), length(ranges[[1]]))
# Adjust each parameter range to be the corrected one (instead of [0,1]) and add columns names to identify the parameters
initial_points <- setNames(data.frame(t(apply(initial_LHS, 1, 
                                              function(x) x*unlist(lapply(ranges[[1]], function(x) x[2]-x[1])) + 
                                                unlist(lapply(ranges[[1]], function(x) x[1]))))), names(ranges[[1]]))

initial_results <- calibration_function(initial_points,prev,pop,timerun)

wave_data[[1]]<-cbind(initial_points, initial_results)

t_sample <- sample(1:nrow(wave_data[[1]]), round(length(wave_data[[1]][,1])/2))
training <- wave_data[[1]][t_sample,]
validation <- wave_data[[1]][-t_sample,]

ems[[1]] <- emulator_from_data(training, names(targets), ranges[[1]])

emplot <- emulator_plot(ems[[k]], cb=TRUE)
png(paste("./Output_",prev,"/EmulatorPlot_",k,".png",sep=''),width=1200, height=800)
print(emplot)
dev.off()

#emplotsd <- emulator_plot(ems[[k]], 'sd', cb=TRUE)
#png(paste("./Output_",prev,"/EmulatorPlotSD_",k,".png",sep=''),width=1200, height=800)
#print(emplotsd)
#dev.off()

# Plot the three validation tests for each emulator in `ems[[k]]`
#validation_diagnostics(ems[[k]], validation = validation, targets = targets, plt=TRUE)

#  classification_diag(ems[[1]],targets,validation)
#  validation_pairs(ems[[1]], validation, targets, cb=TRUE)

non_imp_list <- c(ems[[1]])
non_imp_pts[[1]] <- generate_new_design(non_imp_list, 20 * length(ranges[[1]]), targets, verbose=TRUE)

fitpct <- c(1,"NA","NA")
save(ranges,file=paste("./Output_",prev,"/ranges_",k,".Rdata",sep=''))
save(non_imp_pts,file=paste("./Output_",prev,"/non_imp_pts_",k,".Rdata",sep=''))
save(ems,file=paste("./Output_",prev,"/ems_",k,".Rdata",sep=''))

## Subsequent waves
#load(file=paste("./Output_",prev,"/ranges.Rdata",sep=''))
#load(file=paste("./Output_",prev,"/non_imp_list.Rdata",sep=''))
#load(file=paste("./Output_",prev,"/non_imp_pts.Rdata",sep=''))

# The code below can be used for all waves after wave 1. Start with k equal to 2 to perform 
# the second wave, then change k to 3 to perform the third wave, and so on.
k <- k+1
for (k in 2:40){
# Define the new parameters' ranges
# Wave-k emulators will be trained only on the non-implausible region found in wave (k-1). 
# To do this, we define new ranges for the parameters, identifying the smallest hyper-rectangle 
# containing all points in `non_imp_pts[[k-1]]`. When defining such a hyper-rectangle, we first 
# find the minimum and maximum value of each parameter and then add (resp. subtract) 5% of the 
# obtained range to the maximum (resp. minimum). This is to provide a safety margin, and help 
# ensure that we do not discard any non-implausible point. The new ranges are then stored in
#`ranges[[k]]`. 
min_val <- list()
max_val <- list()
ranges[[k]] = ranges[[k-1]]
for (i in 1:length(ranges[[k-1]])) {
    par <- names(ranges[[1]])[[i]]
    min_val[[par]] <- max(min(non_imp_pts[[k-1]][,par])-0.05*diff(range(non_imp_pts[[k-1]][,par])), 
                          ranges[[1]][[par]][1])
    max_val[[par]] <- min(max(non_imp_pts[[k-1]][,par])+0.05*diff(range(non_imp_pts[[k-1]][,par])),
                          ranges[[1]][[par]][2])
    ranges[[k]][[par]] <- c(min_val[[par]], max_val[[par]])
}

results_k <- calibration_function(non_imp_pts[[k-1]],prev,pop,timerun)

wave_data[[k]]<-cbind(non_imp_pts[[k-1]], results_k)
t_sample <- sample(1:nrow(wave_data[[k]]), round(length(wave_data[[k]][,1])/2))
training <- wave_data[[k]][t_sample,]
validation <- wave_data[[k]][-t_sample,]

ems[[k]] <- emulator_from_data(training, names(targets), ranges[[k]])

emplot <- emulator_plot(ems[[k]], cb=TRUE)
png(paste("./Output_",prev,"/EmulatorPlot_",k,".png",sep=''),width=1200, height=800)
print(emplot)
dev.off()

emplotsd <- emulator_plot(ems[[k]], 'sd', cb=TRUE)
png(paste("./Output_",prev,"/EmulatorPlotSD_",k,".png",sep=''),width=1200, height=800)
print(emplotsd)
dev.off()

#vd <- validation_diagnostics(ems[[k]], validation = validation, targets = targets, plt=TRUE)

ems[[k]]$P1 = ems[[k]]$P1$mult_sigma(0.9)

non_imp_list <- c(ems[[k]], non_imp_list)

non_imp_pts[[k]] <- generate_new_design(non_imp_list, 20 * length(ranges[[k]]), targets, verbose=TRUE)

save(ranges,file=paste("./Output_",prev,"/ranges_",k,".Rdata",sep=''))
save(non_imp_list,file=paste("./Output_",prev,"/non_imp_list_",k,".Rdata",sep=''))
save(non_imp_pts,file=paste("./Output_",prev,"/non_imp_pts_",k,".Rdata",sep=''))
save(ems,file=paste("./Output_",prev,"/ems_",k,".Rdata",sep=''))

if (k > 3){
  simplot <- simulator_plot(list(wave_data[[k-2]], wave_data[[k-1]], wave_data[[k]]), targets)
  png(paste("./Output_",prev,"/SimulatorPlot_",k,".png",sep=''),width=1200, height=800)
  print(simplot)
  dev.off()
}

# Cannot run next step without 3+ waves

# Remultiply parameters by 12 to be annual rather than monthly
apply(non_imp_pts[[k-1]], 2, function(x) quantile(x, c(0.025, 0.5, 0.975))) * c(rep(1,3),rep(1,3),rep(1,3),rep(1,3),rep(12,3),rep(12,3),rep(12,3),rep(12,3),rep(12,3),rep(12,3),rep(12,3),rep(12,3),rep(12,3),rep(12,3))

# How many rows from output are within certain criteria (results_k defined above)
# mcmc_pairs prints plots for each parameter pairing 
#mcmcpairs <- mcmc_pairs(test_data_full, diag_fun = "dens", off_diag_fun = "hex", fill = "purple")
#png(paste("./Output_",prev,"/MCMCPairs_",k,".png",sep=''),width=1200, height=800)
#print(mcmcpairs)
#dev.off()

test <- results_k
# Narrows down based on various criteria - could just use most recent prevalence within 10% of prevalence target.
# Have likely removed measures in 202/206/7/etc. (can check whether prevalence is actually flat or not)
# Can this be run with different prevalence targets rather than running everything above for each target?
test <- test[which(test$P1 > prev*0.9 & test$P1 < prev*1.1),]

# Tells how many runs fit - check whether this is ~10% from hmer
# k
# nrow(test)/(20*length(ranges[[1]]))
fitpct <- cbind(fitpct,c(k,length(test),length(test)/(20*length(ranges[[1]]))))
write.csv(fitpct,file=paste("./Output_",prev,"/fitpct.csv",sep=''))

k <- k+1
}
