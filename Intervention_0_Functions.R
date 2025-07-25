#------------------------------------------------------------------------------#
# Intervention_0_Functions.R                                                   #
# Define functions to run interventions for calibrated model                   #
# Last updated 2025-01-17 by KCH                                               #
#------------------------------------------------------------------------------#

# Define function to sample probability of a positive screening test for each state
uncert_screen <- function(test){
  val_S   <- runif(1,test$S_min  ,test$S_max)
  val_I   <- runif(1,test$I_min  ,test$I_max)
  val_nTB <- runif(1,test$nTB_min,test$nTB_max)
  val_aTB <- runif(1,test$aTB_min,test$aTB_max)
  val_sTB <- runif(1,test$sTB_min,test$sTB_max)
  val_Rn  <- runif(1,test$Rn_min ,test$Rn_max)
  val_Rt  <- runif(1,test$Rt_min ,test$Rt_max)
  
  values <- c(pos_screen_S   = val_S,
              pos_screen_I   = val_I,
              pos_screen_nTB = val_nTB,
              pos_screen_aTB = val_aTB,
              pos_screen_sTB = val_sTB,
              pos_screen_Rn  = val_Rn,
              pos_screen_Rt  = val_Rt)
  
  return(values)
}

# Define function to sample probability of a positive screening test for each state
uncert_confirm <- function(test){
  val_S   <- runif(1,test$S_min, test$S_max)
  val_I   <- runif(1,test$I_min, test$I_max)
  val_nTB <- runif(1,test$nTB_min,test$nTB_max)
  val_aTB <- runif(1,test$aTB_min,test$aTB_max)
  val_sTB <- runif(1,test$sTB_min,test$sTB_max)
  val_Rn  <- runif(1,test$Rn_min,test$Rn_max)
  val_Rt  <- runif(1,test$Rt_min,test$Rt_max)
  
  values <- c(pos_confirm_S   = val_S,
              pos_confirm_I   = val_I,
              pos_confirm_nTB = val_nTB,
              pos_confirm_aTB = val_aTB,
              pos_confirm_sTB = val_sTB,
              pos_confirm_Rn  = val_Rn,
              pos_confirm_Rt  = val_Rt)
  
  return(values)
}

# Model equations
# S = susceptible, I = infected, nTB = non-infectious disease, aTB = asymptomatic disease, sTB = symptomatic disease, 
# Tp = on treatment through passive case detection, 
# Ti = on treatment through intervention from infection, 
# Td = on treatment through intervention from disease, 
# Tr = on treatment through intervention from recovered, 
# Tt = on treatment through intervention from treated, 
# Rn = recovered from non-infectious, Rt = recovered from treatment
# Suffix "x" indicates screening compartment
intervention_model <- function(times, states, parameters) {
  with(as.list(c(states, parameters)), {
    # Sum all compartments for total population N
    N = S  + I  + nTB  + aTB  + sTB  + Tp  + Ti  + Td  + Tr  + Tt  + Rn  + Rt + 
        Sx + Ix + nTBx + aTBx + sTBx + Tpx + Tix + Tdx + Trx + Ttx + Rnx + Rtx
    
    # Equations for main model compartments
    dS   = - (beta * (((t * aTB + sTB) + (t * aTBx + sTBx)) / N) + mu) * S + infclr * I + trtrec * Ti + mu * N + mustb * (sTB + sTBx)
    dI   = - (infclr + infntb + infatb + mu) * I + beta * (((t * aTB + sTB) + (t * aTBx + sTBx)) / N) * (S + p * Rn + r * Rt)
    dnTB = - (ntbrec + ntbatb + mu) * nTB + infntb * I + atbntb * aTB
    daTB = - (atbntb + atbstb + mu) * aTB + infatb * I + ntbatb * nTB + stbatb * sTB
    dsTB = - (stbatb + stbtrt + mustb + mu) * sTB + atbstb * aTB 
    dTp  = - (trtrec + mu) * Tp + stbtrt * sTB 
    dTi  = - (trtrec + mu) * Ti
    dTd  = - (trtrec + mu) * Td
    dTr  = - (trtrec + mu) * Tr
    dTt  = - (trtrec + mu) * Tt
    dRn  = - (beta * (((t * aTB + sTB) + (t * aTBx + sTBx)) / N) * p + mu) * Rn + ntbrec * nTB + trtrec * Tr
    dRt  = - (beta * (((t * aTB + sTB) + (t * aTBx + sTBx)) / N) * r + mu) * Rt + trtrec * (Tp + Td + Tt)
    
    # Equations for screening compartments
    dSx   = - (beta * (((t * aTB + sTB) + (t * aTBx + sTBx)) / N) + mu) * Sx + infclr * Ix + trtrec * Tix
    dIx   = - (infclr + infntb + infatb + mu) * Ix + beta * (((t * aTB + sTB) + (t * aTBx + sTBx)) / N) * (Sx + p * Rnx + r * Rtx)
    dnTBx = - (ntbrec + ntbatb + mu) * nTBx + infntb * Ix + atbntb * aTBx
    daTBx = - (atbntb + atbstb + mu) * aTBx + infatb * Ix + ntbatb * nTBx + stbatb * sTBx
    dsTBx = - (stbatb + stbtrt + mustb + mu) * sTBx + atbstb * aTBx 
    dTpx  = - (trtrec + mu) * Tpx + stbtrt * sTBx 
    dTix  = - (trtrec + mu) * Tix
    dTdx  = - (trtrec + mu) * Tdx
    dTrx  = - (trtrec + mu) * Trx
    dTtx  = - (trtrec + mu) * Ttx
    dRnx  = - (beta * (((t * aTB + sTB) + (t * aTBx + sTBx)) / N) * p + mu) * Rnx + ntbrec * nTBx + trtrec * Trx
    dRtx  = - (beta * (((t * aTB + sTB) + (t * aTBx + sTBx)) / N) * r + mu) * Rtx + trtrec * (Tpx + Tdx + Ttx)

    # Additional compartments to track specific outcomes
    dIs   = - Is + (beta * (((t * aTB + sTB) + (t * aTBx + sTBx)) / N)) * (S + Sx)             # Infections from susceptible
    dIr   = - Ir + (beta * (((t * aTB + sTB) + (t * aTBx + sTBx)) / N) * p + mu) * (Rn + Rnx)  # Infections from recovered
    dIt   = - It + (beta * (((t * aTB + sTB) + (t * aTBx + sTBx)) / N) * r + mu) * (Rt + Rtx)  # Infections from treated
    dRi   = - Ri + infclr * (I + Ix)                                                           # Recovered from infection (cleared)
    dMsTB = - MsTB + mustb * (sTB + sTBx)                                                      # Deaths attributable to TB
    dIc   = - Ic + atbstb * (aTB + aTBx)                                                       # Incident sTB

    dScrNegS   = 0  # Number who test negative on screening test from each compartment
    dScrNegI   = 0
    dScrNegnTB = 0
    dScrNegaTB = 0
    dScrNegsTB = 0
    dScrNegRn  = 0
    dScrNegRt  = 0
    
    dConPosS   = 0  # Number who test positive on confirmatory test (after screening positive) from each compartment
    dConPosI   = 0 
    dConPosnTB = 0
    dConPosaTB = 0
    dConPossTB = 0
    dConPosRn  = 0
    dConPosRt  = 0

    dConNegS   = 0  # Number who test negative on confirmatory test (after screening positive) from each compartment
    dConNegI   = 0
    dConNegnTB = 0
    dConNegaTB = 0
    dConNegsTB = 0
    dConNegRn  = 0
    dConNegRt  = 0
        
    list(c(dS , dI , dnTB , daTB , dsTB , dTp , dTi , dTd , dTr , dTt , dRn , dRt ,
           dSx, dIx, dnTBx, daTBx, dsTBx, dTpx, dTix, dTdx, dTrx, dTtx, dRnx, dRtx, 
           dIs, dIr, dIt  , dRi  , dMsTB, dIc,
           dScrNegS, dScrNegI, dScrNegnTB, dScrNegaTB, dScrNegsTB, dScrNegRn, dScrNegRt,
           dConPosS, dConPosI, dConPosnTB, dConPosaTB, dConPossTB, dConPosRn, dConPosRt,
           dConNegS, dConNegI, dConNegnTB, dConNegaTB, dConNegsTB, dConNegRn, dConNegRt))
  })
}

# Define function to run intervention for each time step in the intervention period and return everyone 
# in screening compartments to respective main model compartments at end of intervention period
intervention <- function(times, states, parameters) {
  with(as.list(c(states, parameters)),{
    intervention_step <- times %% 12
    prop <- 12 - intervention_step
    
    if(intervention_step == 0 && Sx != 0){
    # Used to return everyone in screening compartments to respective main model compartments at end of intervention period
      
      S_hold   = S
      I_hold   = I 
      nTB_hold = nTB
      aTB_hold = aTB
      sTB_hold = sTB
      Tp_hold  = Tp
      Ti_hold  = Ti
      Td_hold  = Td
      Tr_hold  = Tr
      Tt_hold  = Tt
      Rn_hold  = Rn
      Rt_hold  = Rt

      S   = S   + Sx
      I   = I   + Ix
      nTB = nTB + nTBx
      aTB = aTB + aTBx
      sTB = sTB + sTBx
      Tp  = Tp  + Tpx
      Ti  = Ti  + Tix
      Td  = Td  + Tdx
      Tr  = Tr  + Trx
      Tt  = Tt  + Ttx
      Rn  = Rn  + Rnx
      Rt  = Rt  + Rtx
      
      # Set screening compartments to 0
      Sx   = 0
      Ix   = 0
      nTBx = 0
      aTBx = 0
      sTBx = 0
      Tpx  = 0
      Tix  = 0
      Tdx  = 0
      Trx  = 0
      Ttx  = 0
      Rnx  = 0
      Rtx  = 0
      
      # No changes to count compartments
      Is   = Is
      Ir   = Ir
      It   = It
      Ri   = Ri
      MsTB = MsTB
      Ic   = Ic
      
      # No change to number who test negative on screening test for each state
      ScrNegS   = ScrNegS
      ScrNegI   = ScrNegI 
      ScrNegnTB = ScrNegnTB
      ScrNegaTB = ScrNegaTB
      ScrNegsTB = ScrNegsTB
      ScrNegRn  = ScrNegRn
      ScrNegRt  = ScrNegRt
      
      # No change to number who test positive on confirmatory test for each state
      ConPosS   = ConPosS
      ConPosI   = ConPosI 
      ConPosnTB = ConPosnTB
      ConPosaTB = ConPosaTB
      ConPossTB = ConPossTB
      ConPosRn  = ConPosRn
      ConPosRt  = ConPosRt
            
      # No change to number who test negative on confirmatory test for each state
      ConNegS   = ConNegS
      ConNegI   = ConNegI 
      ConNegnTB = ConNegnTB
      ConNegaTB = ConNegaTB
      ConNegsTB = ConNegsTB
      ConNegRn  = ConNegRn
      ConNegRt  = ConNegRt
      
    } else {
    # Used to run intervention for each time step in the intervention period
      
      # Define population to be screened from each compartment at each time step
      # population to screen = compartment size * coverage / proportion screened each month
      S_screen   = (S   * popcov) / prop
      I_screen   = (I   * popcov) / prop
      nTB_screen = (nTB * popcov) / prop
      aTB_screen = (aTB * popcov) / prop
      sTB_screen = (sTB * popcov) / prop
      Tp_screen  = (Tp  * popcov) / prop
      Ti_screen  = (Ti  * popcov) / prop
      Td_screen  = (Td  * popcov) / prop
      Tr_screen  = (Tr  * popcov) / prop
      Tt_screen  = (Tt  * popcov) / prop
      Rn_screen  = (Rn  * popcov) / prop
      Rt_screen  = (Rt  * popcov) / prop
        
      # (1/12th of population)*(1-coverage) moves to the screening compartments without being screened
      S_move   = (S   * (1 - popcov)) / prop
      I_move   = (I   * (1 - popcov)) / prop
      nTB_move = (nTB * (1 - popcov)) / prop
      aTB_move = (aTB * (1 - popcov)) / prop
      sTB_move = (sTB * (1 - popcov)) / prop
      Tp_move  = (Tp  * (1 - popcov)) / prop
      Ti_move  = (Ti  * (1 - popcov)) / prop
      Td_move  = (Td  * (1 - popcov)) / prop
      Tr_move  = (Tr  * (1 - popcov)) / prop
      Tt_move  = (Tt  * (1 - popcov)) / prop
      Rn_move  = (Rn  * (1 - popcov)) / prop
      Rt_move  = (Rt  * (1 - popcov)) / prop
      
      # Remove those who have been screened or moved from respective compartments
      S   = S   - S_move   - S_screen
      I   = I   - I_move   - I_screen
      nTB = nTB - nTB_move - nTB_screen
      aTB = aTB - aTB_move - aTB_screen
      sTB = sTB - sTB_move - sTB_screen
      Tp  = Tp  - Tp_move  - Tp_screen
      Ti  = Ti  - Ti_move  - Ti_screen
      Td  = Td  - Td_move  - Td_screen
      Tr  = Tr  - Tr_move  - Tr_screen
      Tt  = Tt  - Tt_move  - Tt_screen
      Rn  = Rn  - Rn_move  - Rn_screen
      Rt  = Rt  - Rt_move  - Rt_screen
      
      # Add those who have been screened or moved to respective compartments
      # Also move those who screen positive to respective treatment compartments
      Sx   = Sx   + S_move   + S_screen   * ((1 - pos_screen_S)   + pos_screen_S   * (1 - pos_confirm_S)) 
      Ix   = Ix   + I_move   + I_screen   * ((1 - pos_screen_I)   + pos_screen_I   * (1 - pos_confirm_I)) 
      nTBx = nTBx + nTB_move + nTB_screen * ((1 - pos_screen_nTB) + pos_screen_nTB * (1 - pos_confirm_nTB))
      aTBx = aTBx + aTB_move + aTB_screen * ((1 - pos_screen_aTB) + pos_screen_aTB * (1 - pos_confirm_aTB))
      sTBx = sTBx + sTB_move + sTB_screen * ((1 - pos_screen_sTB) + pos_screen_sTB * (1 - pos_confirm_sTB))
      Tpx  = Tpx + Tp_move + Tp_screen  # Individuals on treatment are not screened
      Tdx  = Tdx + Td_move + Td_screen +
             nTB_screen * (pos_screen_nTB * pos_confirm_nTB) +
             aTB_screen * (pos_screen_aTB * pos_confirm_aTB) +
             sTB_screen * (pos_screen_sTB * pos_confirm_sTB)
      Tix  = Tix + Ti_move + Ti_screen +                                                                       
             S_screen * (pos_screen_S * pos_confirm_S) + 
             I_screen * (pos_screen_I * pos_confirm_I)
      Trx  = Trx + Tr_move + Tr_screen +                                                                       
             Rn_screen * (pos_screen_Rn * pos_confirm_Rn)
      Ttx  = Ttx + Tt_move + Tt_screen +                                                                       
             Rt_screen * (pos_screen_Rt * pos_confirm_Rt)
      Rnx  = Rnx + Rn_move + Rn_screen * ((1 - pos_screen_Rn) + pos_screen_Rn * (1 - pos_confirm_Rn))
      Rtx  = Rtx + Rt_move + Rt_screen * ((1 - pos_screen_Rt) + pos_screen_Rt * (1 - pos_confirm_Rt))
      
      # No changes to count compartments
      Is   = Is
      Ir   = Ir
      It   = It
      Ri   = Ri
      MsTB = MsTB
      Ic   = Ic
      
      # Update number who test negative on screening test for each state
      ScrNegS   = ScrNegS   + S_screen   * (1 - pos_screen_S  )
      ScrNegI   = ScrNegI   + I_screen   * (1 - pos_screen_I  )
      ScrNegnTB = ScrNegnTB + nTB_screen * (1 - pos_screen_nTB)
      ScrNegaTB = ScrNegaTB + aTB_screen * (1 - pos_screen_aTB)
      ScrNegsTB = ScrNegsTB + sTB_screen * (1 - pos_screen_sTB)
      ScrNegRn  = ScrNegRn  + Rn_screen  * (1 - pos_screen_Rn )
      ScrNegRt  = ScrNegRt  + Rt_screen  * (1 - pos_screen_Rt )
      
      # Update number who test positive on confirmatory test for each state
      ConPosS   = ConPosS   + S_screen   * (pos_screen_S   * pos_confirm_S  )
      ConPosI   = ConPosI   + I_screen   * (pos_screen_I   * pos_confirm_I  )
      ConPosnTB = ConPosnTB + nTB_screen * (pos_screen_nTB * pos_confirm_nTB)
      ConPosaTB = ConPosaTB + aTB_screen * (pos_screen_aTB * pos_confirm_aTB)
      ConPossTB = ConPossTB + sTB_screen * (pos_screen_sTB * pos_confirm_sTB)
      ConPosRn  = ConPosRn  + Rn_screen  * (pos_screen_Rn  * pos_confirm_Rn )
      ConPosRt  = ConPosRt  + Rt_screen  * (pos_screen_Rt  * pos_confirm_Rt )
      
      # Update number who test negative on confirmatory test for each state
      ConNegS   = ConNegS   + S_screen   * (pos_screen_S   * (1 - pos_confirm_S  ))
      ConNegI   = ConNegI   + I_screen   * (pos_screen_I   * (1 - pos_confirm_I  ))
      ConNegnTB = ConNegnTB + nTB_screen * (pos_screen_nTB * (1 - pos_confirm_nTB))
      ConNegaTB = ConNegaTB + aTB_screen * (pos_screen_aTB * (1 - pos_confirm_aTB))
      ConNegsTB = ConNegsTB + sTB_screen * (pos_screen_sTB * (1 - pos_confirm_sTB))
      ConNegRn  = ConNegRn  + Rn_screen  * (pos_screen_Rn  * (1 - pos_confirm_Rn ))
      ConNegRt  = ConNegRt  + Rt_screen  * (pos_screen_Rt  * (1 - pos_confirm_Rt ))
    }
        
    return(c(S , I , nTB , aTB , sTB , Tp , Ti , Td , Tr , Tt , Rn , Rt ,
             Sx, Ix, nTBx, aTBx, sTBx, Tpx, Tix, Tdx, Trx, Ttx, Rnx, Rtx,
             Is, Ir, It  , Ri  , MsTB, Ic,
             ScrNegS, ScrNegI, ScrNegnTB, ScrNegaTB, ScrNegsTB, ScrNegRn, ScrNegRt,
             ConPosS, ConPosI, ConPosnTB, ConPosaTB, ConPossTB, ConPosRn, ConPosRt,
             ConNegS, ConNegI, ConNegnTB, ConNegaTB, ConNegsTB, ConNegRn, ConNegRt))
  })
}

# Define function to run no-intervention scenario
run_base_set <- function(model_parameters,screen_parameters=NA,confirm_parameters=NA,coverage_parameters,pop,prev,intyears){
  
  # Define initial state values 
  states <- c(S    = pop*.75 - prev*5,
              I    = pop*.25,
              nTB  = prev*4,
              aTB  = prev/2,
              sTB  = prev/2,
              Tp   = 0,
              Ti   = 0,
              Td   = 0,
              Tr   = 0,
              Tt   = 0,
              Rn   = 0,
              Rt   = 0,
              Sx   = 0,
              Ix   = 0,
              nTBx = 0,
              aTBx = 0,
              sTBx = 0,
              Tpx  = 0,
              Tix  = 0,
              Tdx  = 0,
              Trx  = 0,
              Ttx  = 0,
              Rnx  = 0,
              Rtx  = 0,
              Is   = 0,
              Ir   = 0,
              It   = 0,
              Ri   = 0,
              MsTB = 0,
              Ic   = 0,
              ScrNegS   = 0,
              ScrNegI   = 0,
              ScrNegnTB = 0,
              ScrNegaTB = 0,
              ScrNegsTB = 0,
              ScrNegRn  = 0,
              ScrNegRt  = 0,
              ConPosS   = 0, 
              ConPosI   = 0, 
              ConPosnTB = 0, 
              ConPosaTB = 0,
              ConPossTB = 0, 
              ConPosRn  = 0, 
              ConPosRt  = 0, 
              ConNegS   = 0, 
              ConNegI   = 0, 
              ConNegnTB = 0, 
              ConNegaTB = 0, 
              ConNegsTB = 0, 
              ConNegRn  = 0, 
              ConNegRt  = 0
  )
  
  output <- as.data.frame(ode(func  = intervention_model,  # Define function as intervention_model
                              y     = states,              # Define initial state values
                              parms = model_parameters,    # Define parameter values
                              times = seq(0,5000,1)))      # Define length of model run
  return(c(output,screen_parameters,confirm_parameters))
}

# Define function to run intervention scenarios
# To improve speed, the first 4670 months are copied from the no-intervention output
# Intervention models run from month 4670 with interventions beginning in month 4680
run_int_set <- function(baseset,model_parameters,screen_parameters=NA,confirm_parameters=NA,coverage_parameters,pop,prev,intyears){
  
  # Define initial state values as those from month 4670 of the no-intervention model
  states <- c(S    = baseset$S[4670],
              I    = baseset$I[4670],
              nTB  = baseset$nTB[4670],
              aTB  = baseset$aTB[4670],
              sTB  = baseset$sTB[4670],
              Tp   = baseset$Tp[4670],
              Ti   = baseset$Ti[4670],
              Td   = baseset$Td[4670],
              Tr   = baseset$Tr[4670],
              Tt   = baseset$Tt[4670],
              Rn   = baseset$Rn[4670],
              Rt   = baseset$Rt[4670],
              Sx   = baseset$Sx[4670],
              Ix   = baseset$Ix[4670],
              nTBx = baseset$nTBx[4670],
              aTBx = baseset$aTBx[4670],
              sTBx = baseset$sTBx[4670],
              Tpx  = baseset$Tpx[4670],
              Tix  = baseset$Tix[4670],
              Tdx  = baseset$Tdx[4670],
              Trx  = baseset$Trx[4670],
              Ttx  = baseset$Ttx[4670],
              Rnx  = baseset$Rnx[4670],
              Rtx  = baseset$Rtx[4670],
              Is   = baseset$Is[4670],
              Ir   = baseset$Ir[4670],
              It   = baseset$It[4670],
              Ri   = baseset$Ri[4670],
              MsTB = baseset$MsTB[4670],
              Ic   = baseset$Ic[4670],
              ScrNegS   = baseset$ScrNegS[4670],
              ScrNegI   = baseset$ScrNegI[4670],
              ScrNegnTB = baseset$ScrNegnTB[4670],
              ScrNegaTB = baseset$ScrNegaTB[4670],
              ScrNegsTB = baseset$ScrNegsTB[4670],
              ScrNegRn  = baseset$ScrNegRn[4670],
              ScrNegRt  = baseset$ScrNegRt[4670],
              ConPosS   = baseset$ConPosS[4670],
              ConPosI   = baseset$ConPosI[4670],
              ConPosnTB = baseset$ConPosnTB[4670],
              ConPosaTB = baseset$ConPosaTB[4670],
              ConPossTB = baseset$ConPossTB[4670],
              ConPosRn  = baseset$ConPosRn[4670],
              ConPosRt  = baseset$ConPosRt[4670],
              ConNegS   = baseset$ConNegS[4670],
              ConNegI   = baseset$ConNegI[4670],
              ConNegnTB = baseset$ConNegnTB[4670],
              ConNegaTB = baseset$ConNegaTB[4670],
              ConNegsTB = baseset$ConNegsTB[4670],
              ConNegRn  = baseset$ConNegRn[4670],
              ConNegRt  = baseset$ConNegRt[4670]
  )
  
  if (intyears == 1){           # Run intervention with screening in Year 1
    intoutput <- as.data.frame(ode(func   = intervention_model,                                                            # Define function as intervention_model
                                   y      = states,                                                                        # Define initial state values
                                   parms  = c(model_parameters,screen_parameters,confirm_parameters,coverage_parameters),  # Define parameter values
                                   times  = seq(4670,5000,1),                                                              # Define length of model run
                                   events = list(func = intervention, time = c(seq(4681, 4692, 1)))))                      # "Events" option" is used to run intervention
  } else if (intyears == 2){    # Run intervention with screening in Years 1, 2
    intoutput <- as.data.frame(ode(func   = intervention_model, 
                                   y      = states, 
                                   parms  = c(model_parameters,screen_parameters,confirm_parameters,coverage_parameters),
                                   times  = seq(4670,5000,1),
                                   events = list(func = intervention, time = c(seq(4681, 4692, 1),
                                                                               seq(4693, 4704, 1)))))
  } else if (intyears == 3){    # Run intervention with screening in Years 1, 2, 3
    intoutput <- as.data.frame(ode(func   = intervention_model, 
                                   y      = states, 
                                   parms  = c(model_parameters,screen_parameters,confirm_parameters,coverage_parameters),
                                   times  = seq(4670,5000,1),
                                   events = list(func = intervention, time = c(seq(4681, 4692, 1),
                                                                               seq(4693, 4704, 1),
                                                                               seq(4705, 4716, 1)))))
  } else if (intyears == 4){    # Run intervention with screening in Years 1, 2, 3, 4
    intoutput <- as.data.frame(ode(func   = intervention_model, 
                                   y      = states, 
                                   parms  = c(model_parameters,screen_parameters,confirm_parameters,coverage_parameters),
                                   times  = seq(4670,5000,1),
                                   events = list(func = intervention, time = c(seq(4681, 4692, 1),
                                                                               seq(4693, 4704, 1),
                                                                               seq(4705, 4716, 1),
                                                                               seq(4717, 4728, 1)))))
  } else if (intyears == 5){    # Run intervention with screening in Years 1, 2, 3, 4, 5
    intoutput <- as.data.frame(ode(func   = intervention_model, 
                                   y      = states, 
                                   parms  = c(model_parameters,screen_parameters,confirm_parameters,coverage_parameters),
                                   times  = seq(4670,5000,1),
                                   events = list(func = intervention, time = c(seq(4681, 4692, 1),
                                                                               seq(4693, 4704, 1),
                                                                               seq(4705, 4716, 1),
                                                                               seq(4717, 4728, 1),
                                                                               seq(4729, 4740, 1)))))
  }
  
  output <- rbind(baseset[1:4670,],intoutput)
  return(c(output,screen_parameters,confirm_parameters))
}
