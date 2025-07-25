#------------------------------------------------------------------------------#
# Calibration_0_Functions.R                                                    #
# Define functions to run calibration                                          #
# Last updated 2025-01-17 by KCH                                               #
#------------------------------------------------------------------------------#

# Define model structure and trackers
transmission_model <- function(times, states, parameters) {
  with(as.list(c(states, parameters)), {
    N = S + I + nTB + aTB + sTB + Tp + Rm + Rt  
    
    dS   = - (beta * ((t * aTB + sTB) / N) + mu) * S + infclr * I + mu * N + mus * sTB
    dI   = - (infclr + infmin + infsub + mu) * I + beta * ((t * aTB + sTB) / N) * (S + p * Rm + r * Rt)
    dnTB = - (minrec + minsub + mu) * nTB + infmin * I + submin * aTB
    daTB = - (submin + subclin + mu) * aTB + infsub * I + minsub * nTB + clinsub * sTB
    dsTB = - (clinsub + clintrt + mus + mu) * sTB + subclin * aTB
    dTp  = - (trtrec + mu) * Tp + clintrt * sTB
    dRm  = - (beta * ((t * aTB + sTB) / N) * p + mu) * Rm + minrec * nTB
    dRt  = - (beta * ((t * aTB + sTB) / N) * r + mu) * Rt + trtrec * Tp
    
    dIs  = - Is + (beta * ((t * aTB + sTB) / N)) * S 
    dIr  = - Ir + (beta * ((t * aTB + sTB) / N) * p + mu) * Rm
    dIt  = - It + (beta * ((t * aTB + sTB) / N) * r + mu) * Rt
    dRi  = infclr * I
    dMs  = mus * sTB
    
    list(c(dS , dI , dnTB, daTB, dsTB, dTp, dRm, dRt, 
           dIs, dIr, dIt , dRi , dMs))
  })
}

# Define calibration function
calibration_function <- function(parameters,prev,pop,timerun){
  states <- c(S   = pop*.75 - prev*5,
              I   = pop*.25,
              nTB = prev*4,
              aTB = prev/2,
              sTB = prev/2,
              Tp  = 0,
              Rm  = 0,
              Rt  = 0,
              Is  = 0,
              Ir  = 0,
              It  = 0,
              Ri  = 0,
              Ms  = 0)
    
  output <- as.data.frame(matrix(nrow = nrow(parameters), ncol = 1))
  names(output) <- c("P1")
  
  for(i in 1:nrow(parameters)){
      print(i)
      params <- parameters[i,]
      names(params) <- names(parameters[i,])
      params <- c(params, mu = 0.0137/12, trtrec = 1/6)
      run <- as.data.frame(ode(func=transmission_model,y=states,parms=params,times=seq(0,timerun,1)))
      output[i,1] <- run$aTB[4200] + run$sTB[4200]
  }
  
  return(output)

}

run_params <- function(model_parameters,prev,pop,timerun){
  states <- c(S   = pop*.75 - prev*5,
              I   = pop*.25,
              nTB = prev*4,
              aTB = prev/2,
              sTB = prev/2,
              Tp  = 0,
              Rm  = 0,
              Rt  = 0,
              Is  = 0,
              Ir  = 0,
              It  = 0,
              Ri  = 0,
              Ms  = 0)
  output <- ode(func  = transmission_model, 
                y     = states, 
                parms = model_parameters,
                times = seq(0,timerun,1))
}