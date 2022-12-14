from cobra.flux_analysis import flux_variability_analysis
import numpy as np
import cobra
from scipy.integrate import solve_ivp
import progressbar

def isfeasible(model):
    sol = model.optimize()
    if sol.status == 'optimal':
        feasibility = 1
    else:
        feasibility = 0
    return feasibility


def IntSigMet(CBmodel,Obj,OdeSig,odeargs,init,dt,T,cal_frac,beta,Se,Sm,SigToReac,FVA=0,MAX=0):

  ##INPUT
  #   CBmodel: COBRA model of the metabolic network 
  #
  #   Obj: list of strings of the reactions whose flux has to be maximised
  #
  #   OdeSig: A function defining the ODE based model of signalling network
  #
  #   odeargs:  Other arguments required for OdeSig (Will be external inputs to the signalling network)
  #
  #   init: A list denoting the initial conditions of the ODE based model of the signalling network
  #
  #   dt: Float value indicating time step of integrating the ODE model
  #
  #   T: Float value indicating total time period of the simulation
  #
  #   cal_frac: A funtion that takes concentration of signals as input and calculates fractional 
  #             concentration values of the signals that affect the metabolic network
  #
  #   beta: numpy array denoting the tuning parameters for integrating metabolic exchange fluxes with signalling network
  #
  #   Se: dictionary that indicates enzymes in metabolic network acting as a signal.
  #       ['signal_ids']: list denoting the indices of the signalling proteins
  #       ['rxn_ids']: list denoting the indices of the reactions catalysed by enzymes that take part in signalling network
  #
  #   Sm: dictionary that indicates metabolites in metabolic network acting as an external input to the signalling network.
  #       ['input_ids']: list denoting the indices of the external input to signalling network
  #       ['rxn_ids']: list denoting the indices of the exchange reactions of the corresponding metabolite
  #
  #   SigToReac: Numpy matrix with n_rows = size of output of cal_frac; n_columns = number of reactions in the metabolic network;
  #              Elements of the matrix (SigToReac_ij) are 1,0 and -1 indicating activation, noeffect and inhibition of the j^th reaction by
  #              the i^th fractional concentration of signalling molecule respectively.
                
  ##OPTIONAL INPUTS        
  #   FVA: Boolean value indicating whether to use flux variability analysis or flux balance analysis for reactions/enzymes having
  #        effect on signalling molecules
  #        true: Use flux variability analysis
  #        false: Use flux balance analysis (default value)
  #
  #   MAX: Boolean value indicating whether to use max or average of fractional concentration values when more than one signalling 
  #        molecule has effect on a reaction
  #        true: maximum of fractional concentration values
  #        false: average of fractional concentration values (default value)

  ##OUTPUT
  #   C_signal: Numpy matrix with n_rows = number of time steps; n_columns = total number of signalling molecules
  #             each element in the matrix indicates concentration of the signalling molecule
  #
  #   F_reactions: Numpy matrix with n_rows = number of time steps; n_columns = total number of reactions
  #             each element in the matrix indicates flux through the corresponding reaction
  #
  #   concUnint: Numpy matrix with n_rows = number of time steps; n_columns = total number of signalling molecules
  #             each element in the matrix indicates concentration of the signalling molecule (without metabolic network integration)
  #
  #   F: Solution to FBA obtained at each iteration
  #   
  #   fracInt: Factional activation values of the signalling molecule obtained at each iteration
  #
  #   lbnew: Modified lower bounds obtained at each time step of integration
  #
  #   ubnew: Modified upper bounds obtained at each time step of integration
  
  
  ##AUTHOR
  #       Pavan Kumar S, BioSystems Engineering and Control (BiSECt) lab, IIT Madras

  CBmodel.objective = []
  fva_df = flux_variability_analysis(CBmodel)
  minFlux = fva_df['minimum'].values
  maxFlux = fva_df['maximum'].values
  CBmodel.objective = Obj

  

  tspan = np.arange(0,T+dt,dt)

  # preallocation of variables
  F = np.zeros(len(tspan)) # value of objective function at every time step
  F_reactions = np.zeros((len(tspan),len(CBmodel.reactions))) # fluxes through the reactions at every time step
  ext_inp = [] # to store the external inputs at each step
  ext_inp.append(odeargs.copy())

  
  lbnew = np.zeros((len(tspan),len(CBmodel.reactions))) # modified lower bounds at each time step
  ubnew = np.zeros((len(tspan),len(CBmodel.reactions))) # modified upper bounds at each time step
 
  fracInt = []

  x = solve_ivp(OdeSig,[0,T],init,args=(odeargs,),method='BDF',t_eval=tspan)
  H = x.y.T
  concUnint = H
  concInt = np.zeros((len(tspan), len(init)))
  H_1 = init
  k = -1
  
  widgets = [' Progress: ',progressbar.Percentage(),progressbar.Bar('*'),] 
  bar = progressbar.ProgressBar(max_value=T+dt,widgets=widgets).start()
  for time in tspan:
    k = k+1
    concInt[k,:] = H_1;
    if time == 0:
        fractSig = cal_frac(H[1,:])
    else:
        fractSig = cal_frac(H_1)
    
    fracInt.append(fractSig)
    # Signalling to metabolic
    for i in range(SigToReac.shape[1]):
        if any(SigToReac[:,i]):
            if MAX:
                tempVar = fractSig[SigToReac[:,i]!=0]
                sigidx = tempVar==max(abs(tempVar))
                frac = fractSig[sigidx]*SigToReac[sigidx,i] 
            elif not MAX:
                sigidx = SigToReac[:,i]!=0
                frac = np.mean(fractSig[sigidx]*SigToReac[sigidx,i])
                    
            model_temp = CBmodel
            if frac>0:
                lbnew[k,i] = minFlux[i]+((maxFlux[i]-minFlux[i])*frac)
                model_temp.reactions[i].lower_bound=lbnew[k,i]
            else:
                ubnew[k,i] = maxFlux[i]+((maxFlux[i]-minFlux[i])*frac)
                model_temp.reactions[i].upper_bound=ubnew[k,i]
            
            if isfeasible(model_temp) == 1:
                if frac>0:
                    lbnew[k,i] = minFlux[i]+((maxFlux[i]-minFlux[i])*frac)
                    CBmodel.reactions[i].lower_bound=lbnew[k,i]
                else:
                    ubnew[k,i] = maxFlux[i]+((maxFlux[i]-minFlux[i])*frac)
                    CBmodel.reactions[i].upper_bound=ubnew[k,i]
    CBmodel.objective = Obj

    soln = CBmodel.optimize()
    F[k] = soln.objective_value
    F_reactions[k,:] = soln.fluxes.values
    
    
    # metabolite to signalling
    if ~FVA:
      if len(Sm) !=0:
          odeargs[Sm['input_ids']] = odeargs[Sm['input_ids']]+(soln.fluxes.values[Sm['rxn_ids']]*beta*dt)
          
      if len(Se) !=0:
          H_1[Se['signal_ids']] = (soln.fluxes.values[Se['rxn_ids']]-minFlux[Se['rxn_ids']])/(maxFlux[Se['rxn_ids']]-minFlux[Se['rxn_ids']]);
      
    elif FVA:
      fva_df1 = flux_variability_analysis(CBmodel)
      minF=fva_df1['minimum'].values
      maxF=fva_df1['maximum'].values
      if len(Sm) !=0:
          odeargs[Sm['input_ids']] = odeargs[Sm['input_ids']]+(soln.fluxes.values[Sm['rxn_ids']]*beta*dt)
      
      if len(Se) !=0:
          first_term = (minF[Se['rxn_ids']]-minFlux[Se['rxn_ids']])/(maxFlux[Se['rxn_ids']]-minFlux[Se['rxn_ids']])
          second_term = (maxF[Se['rxn_ids']]-minFlux[Se['rxn_ids']])/(maxFlux[Se['rxn_ids']]-minFlux[Se['rxn_ids']])
          H_1[Se['signal_ids']] = (first_term+second_term)/2;
    tintspan = [0,dt]
    ext_inp.append(odeargs.copy())
    x = solve_ivp(OdeSig,tintspan,H_1,args=(odeargs,),method='BDF')
    Con = x.y.T
    H_1 = Con[-1,:]
    bar.update(time+dt)
  return concInt, F_reactions, concUnint, F, fracInt, lbnew, ubnew, ext_inp
