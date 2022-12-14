def InsGluSignal(t,H,u):
  '''This code is adapted from 
  "Somvanshi, P. R., Tomar, M., & Kareenhalli, V. (2019). 
  computational Analysis of insulin-Glucagon Signalling network: implications of 
  Bistability to Metabolic Homeostasis and Disease states. Scientific reports, 9(1), 1-11"
  Inputs:
          H: State variables of the signalling network
          t: Time
          u: External inputs to the signalling network
  '''
  import math
  # external inputs
  Gpi,FF,AAm = u[0],u[1],u[2]

  # State variables
  IR      = H[0]
  IRp     = H[1]
  IRS     = H[2]
  IRSp    = H[3]
  PI3K    = H[4]
  PI3Kp   = H[5]
  AKT     = H[6]
  AKTp    = H[7]
  PKC     = H[8]
  PKCp    = H[9]
  GLUT4c  = H[10]
  GLUT4s  = H[11]
  mTOR    = H[12]
  mTORa   = H[13]
  S6K1    = H[14]
  S6K1a   = H[15]
  InsV    = H[16]
  InsL    = H[17]
  InsP    = H[18]
  Glnp    = H[19]
  Glcn    = H[20]
  FR      = H[21]
  RS      = H[22]
  LR      = H[23]
  LRS     = H[24]
  LRP     = H[25]
  DAG     = H[26]
  PKA     = H[27]
  cAMP    = H[28]
  PDE3a   = H[29]
  AMPK    = H[30]
  # Plasma macronutrients 
  Ca_Glu=Gpi
  Ca_FFA=FF
  Ca_Ala=AAm
  
  Gnp=Glnp 
  insulinp=((InsP/0.05)*1e-12)*2*((25e-12)**2/(Gnp**2+(25e-12)**2))
  
  #-----------------------Parameters for Insulin Signalling Model----------------------

  #Nyman, E. et al. A hierarchical whole-body modeling approach elucidates the link between in vitro insulin signaling and in vivo glucose homeostasis. Journal of Biological Chemistry, doi:10.1074/jbc.M110.188987 (2011).
  #Giri, L., Mutalik, V. K. & Venkatesh, K. V. A steady state analysis
  #indicates that negative feedback regulation of PTP1B by Akt elicits bistability in insulin-stimulated GLUT4 translocation. Theoretical Biology and Medical Modelling, doi:10.1186/1742-4682-1-2 (2004)


  kn=7
  nakt = 1.8
  npkc=1
  kakt=4.25*1
  k1aBasic = 137775                                                                   
  k1f = 673.445*(10**10)*1                                                      
  k1b = 543752*1                                                               
  k2f = 141.254                                                                 
  k21f = 9178.23*7    
  k2b = 3533330*1             
  k2c=  3533330
  k4f = 226300*1   #phosphorylation of PI3K by IRSp                                                              
  k4b = 1493960*1     
  k4f1 = 513844*1  #phosphorylation of AKT by PI3Kp                                                                 
  k4b1 = 320822*1 
  k6f = 1*1337290/2                                                            
  k6b = 1*1629.22*20                                                                
  k5f = 0.06308*1                                                         
  k5b = 0.212*1 
  k8f = 856353*1 # mTOR activation by AKTp
  k8f1 = 1*856353/2
  k8b = 1430390*1.5*1
  Vi=0.05
  # S6K negative feedback on IRS1
  S6K_Ntv_IRS = (S6K1a**4)/((6**4.0 + S6K1a**4.0))  

  # rate of phosphorylation of IR 
  v1f = 1*k1f*insulinp*1000*IR + k1aBasic*IR

    
  # rate of dephosphorylation of IRp
  v1b = k1b*IRp        

  PKCp_Ntv_IRS = (kn/(kn+(PKCp/1.5)**npkc)) 

  AKTp_Ptv_IRS = AKTp**nakt/(AKTp**nakt+kakt**nakt)  

  # rate of phosphorylation of IRSp which is activated by AKTp and inhibited by PKCp

  v2f = (k2f + k21f*AKTp_Ptv_IRS)*IRp*IRS*PKCp_Ntv_IRS # normal one
  # v2f = (k2f + 0)*IRp*IRS*PKCp_Ntv_IRS # insulin resistance

  # rate of dephosphorelation of IRSp

  v2b = k2b*IRSp*(1 + 0.2*S6K_Ntv_IRS)

  # rate of phosphorylation of PI3K by IRSp

  v4f = 1*k4f*PI3K*IRSp

  # rate of dephosphorylation of PI3Kp

  v4b = k4b*PI3Kp

  # rate of phosphorylation of AKT by PI3Kp

  v4f1 = k4f1*AKT*PI3Kp

  # rate of dephosphorylation of AKTp

  v4b1 = k4b1*AKTp

  # Positive feedback of FFA on PKC
  FFA=FF
  FFA_Ptv_PKC = 1 + 0.5*(FFA**3/(FFA**3+1.5**3))

  DAG_PKC_Ptv = 1+1*(DAG**3/(7**3+DAG**3)) # normal 
  # DAG_PKC_Ptv = 2  # insulin resistance

  # rate of phosphorylation of PKC by PI3Kp , it is activated by positive feedbacks from free fatty acids (FFA) and DAG

  v6f = k6f*(FFA_Ptv_PKC*DAG_PKC_Ptv)*PKC*PI3Kp

  # rate of dephosphorylation of PKCp

  v6b = k6b*PKCp

  # rate of translocation of 4 from cytosol to cell surface, activated by AKTp and PKCp

  v5f = k5f*(0.2*AKTp+0.8*PKCp)*GLUT4c 

  # rate of translocation of GLUT4 from cell surface to cytosol
  v5b = k5b*GLUT4s

  # mTOR_Raptor activation by AA in absence of insulin
  # Vinod, P. K. U. & Venkatesh, K. V. Quantification of the effect of amino acids on an integrated mTOR and insulin signaling pathway. Molecular BioSystems, doi:10.1039/b816965a (2009).

  AA=AAm
  ns=3.00
  faa=0.75
  AA_mTOR_eff = 2.5*((AA**ns)/(AA**ns + faa**ns))

  # rate of mTOR activation by AKTp and amino acids (AA)

  v8f = k8f*mTOR*AKTp+k8f1*(1+ AA_mTOR_eff)*mTOR

  # rate of mTORa deactivation
  v8b = k8b*mTORa

  # S6K1 phosphorylation-dephosphorylation  
  PP2Amax = 5 
  AMPK_Eff=1*AMPK 
  f20=1.50e-3*1
  f21=0.5e-3*1
  fpp=6
  fto=2
  np=1.2
  nt=4
  PDK1a = 3*(PI3K/(4+PI3K))

  # Regulation factor of PP2A on S6K1-phosphorylation
  PP2A =  PP2Amax*(fpp**np/(fpp**np + mTORa**np))
  AMPK_Ntv_S6K = ((0.75**2.0)/(0.75**2.0+AMPK_Eff**2.0))
  reg_1 = (fto**nt/(fto**nt + PP2A**nt))

  # rate of S6K1 phosphorylation; activated by PI3K & mTORa, inhibited by PP2a and AMPK                       

  v9f = f20*(PDK1a)*(mTORa)*reg_1*(S6K1)*AMPK_Ntv_S6K

  # rate of S6K1a dephosphorylation
  v9b = f21*(PP2A)*(S6K1a)

  IRprime = v1b-v1f                                                             
  IRpprime = -v1b+v1f                                                                
  IRSprime= v2b-v2f                                                                 
  IRSpprime= -v2b+v2f     
  PI3Kprime = v4b - v4f                                       
  PI3Kpprime = v4f - v4b                                        
  AKTprime = v4b1-v4f1                                                                 
  AKTpprime = -v4b1+v4f1 
  PKCprime = v6b - v6f                                        
  PKCpprime = v6f - v6b                                          
  GLUT4cprime = v5b-v5f                                                              
  GLUT4sprime = -v5b+v5f
  mTORprime = v8b - v8f                                        
  mTORaprime = v8f - v8b                                          
  S6K1prime = v9b - v9f                                        
  S6K1aprime = v9f - v9b                                          

  #---Insulin kinetics Parameters--------
  Vi=0.05
  m1=0.190
  m2=0.484
  m4=0.194
  m5=0.0304*Vi
  m6=0.6471
  gamma=0.5

  ## Insulin secretion
  #Somvanshi, P. R., Patel, A. K., Bhartiya, S. & Venkatesh, K. V. Influence of plasma macronutrient levels on hepatic metabolism: role of regulatory networks in homeostasis and disease states. RSC Advances 6, 14344-14371, doi:10.1039/C5RA18128C (2016).

  K_AA=0.6
  K_FF=1.6
  na=5.8
  nf=4.8
  n=4.65
  K_Glu=8.9
  V_Glu=48e-12
  V_Ala=25e-12
  V_FFA=25e-12

  # Insulin secreted by pancreas in response to plasma glucose, alanine and free fatty acids in plasma
  Ins_sec=V_Glu*(Ca_Glu**n/(Ca_Glu**n+K_Glu**n))+(V_Ala*((Ca_Ala**na)/(((Ca_Ala**na)+K_AA**na))))+(V_FFA*(Ca_FFA**nf)/((Ca_FFA**nf)+K_FF**nf))

  ISR=(Ins_sec/1e-12)


  ## InsV
  InsVprime=-(gamma*InsV)+ISR  

  InSec=gamma*InsV
  
  #HE-Hepatic Extraction
  HE=(-m5*InSec)+ m6

  m3=((HE*m1)/(1-HE))
  # InsL
  InsLprime=(-(m1+m3)*InsL)+(m2*InsP)+InSec


  # InsP
  InsPprime=(-(m2+m4)*InsP)+(m1*InsL) # *10^-10


  ################ Glucagon Secretion and kinetics ################
  ################--------Parameters for Glucagon Receptror Model------------################

  C1=1*10*60
  c2=100*60
  c3a=5.2*60e-3
  c4=4*60e-3
  c5=5.2*60e-3
  c8=650*60
  A0=3
  B1=100
  B2=1e6
  V=11
  Vp=3

  
  # Gucagon balance-
  # Liu, W., Hsin, C. C. & Tang, F. A molecular mathematical model of glucose mobilization and uptake. Mathematical Biosciences, doi:10.1016/j.mbs.2009.07.005 (2009).

  Gm=1*1.350e-10
  p1=0.005*180
  q1=10

  if (Ca_Ala-0.25)>=0:
    Am_Glcn=(Ca_Ala-0.25)
  else:
    Am_Glcn=0

  Km_AA=1
  nA=4.5
  Vm_AAg=2*125e-12

  AA_Glcn_Eff=Vm_AAg*((Am_Glcn**nA)/(Am_Glcn**nA+Km_AA**nA))

  if Ca_Glu <=5:
    Glcn_Sec= (Gm/(1*(Ca_Glu/5)+(q1*math.exp(p1*(Ca_Glu-5)))))+AA_Glcn_Eff 
  else:
    Glcn_Sec= (Gm/(1+(q1*math.exp(p1*(Ca_Glu-5)))))+AA_Glcn_Eff 


  # Plasma Glucagon
  a1=0.12*1.25
  a2=0.3*1
  Glnpprime=(-(a1+a2)*Glnp)+Glcn_Sec

  #---------------Glucagon receptor model---------------

  #Riccobene, T. A., Omann, G. M. & Linderman, J. J. Modeling activation and desensitization of G-protein coupled receptors provides insight into ligand efficacy. J Theor Biol, doi:10.1006/jtbi.1999.0988 (1999).

  Gln_recp=9e-13
  a4=6e7
  a1=0.12*1
  a3=0.1
  Rt=126500
  alf=10000
  Glcnprime=(a1*Vp*Gnp/V)-(a3*Glcn)-((a4)*Glcn*(FR/Rt)*Gln_recp)
  GCN=Glcn*1e6
  c3=c3a*(2*(GCN)/(32.0e-6+GCN))

  # The equations given below describe the dynamic model of glucagon binding to G-Protein Coupled Receptors
  # that includes G protein activation and receptor desensitization

  ## FR-free receptor
  FRprime= (C1*LR)-(c2*GCN*FR)-(c3*FR)+(c4*RS)

  ## RS-Sequesterd receptor
  RSprime=(C1*LRS)-(alf*c2*GCN*RS)+(c3*FR)-(c4*RS)

  ## LR-Ligand bound receptor
  LRprime=(c2*GCN*FR)-(C1*LR)+((c4/alf)*LRS)-(c3*LR)+ (c5*LRP)

  # G-protein and calcium dynamics
  # Kummer, U. et al. Switching from simple to complex oscillations in calcium signaling. Biophysical Journal, doi:10.1016/S0006-3495(00)76373-9 (2000).
  # Derived steady state relationships using Kummer et al. (2000) and Somvanshi et al.,(2016)

  # Due to the faster timescales of the calcium signalling with respect to other signalling events,
  # the Hill formulations have been used under pseudo steady state approximations in order to describe 
  # the input-output relationship in form of dose response curve for the above mentioned 4 signalling components – GPRT, PLC, IP3i, Cal

  GPRT=26*(LRS/(95000+LRS))   
  plc=233*(GPRT**2/(32**2+GPRT**2))
  IP3i=16.3*(plc**4/((42)**4)+plc**4)
  Cal=35*(IP3i**1.8)/(9**1.8+IP3i**1.8)

  ## LRS-ligand bound sequestered receptor
  # GPCR
  LRSprime=(c3*LR)-((c4/alf)*LRS)+(alf*c2*GCN*RS)-(C1*LRS)-(c8*(1+((A0*(GPRT))/(B1+(GPRT))))*(LRS/(B2+LRS)))

  ## LRP
  LRPprime=(c8*(1+((A0*(GPRT))/(B1+(GPRT))))*(LRS/(B2+LRS)))- (c5*LRP)

  ## DAG-diacylglycerol
  # DAG is activated by Cal, PKA and Fatty Acids
  vi=1*2*60
  vf=1.00*60
  bd=0.500*60
  ks=1e-4
  Kc2p= 1-1*(PKA**4/(ks**4+PKA**4))
  DAGprime=vi*((0.01*Cal*plc)/(Kc2p + 0.1*Cal))+vf*(FFA)- bd*DAG

  ##########--------CAMP-PKA signaling-------##########
  # Siso-Nadal, F., Fox, J. J., Laporte, S. A., Hébert, T. E. & Swain, P. S. Cross-Talk between Signaling Pathways Can Generate Robust Oscillations in Calcium and cAMP. PLOS ONE 4, e7189, doi:10.1371/journal.pone.0007189 (2009).

  kc1= 2.0*1e-6
  kcm1 = 25.0e-12
  kc2 = 1.18          
  kcm2 = 1.00   
  pn = 1
  pe = 2  
  cqi = 1*2.0 

  Kck = 1 + 0.5*(cqi*(Cal**3/(9**3+Cal**3)))
  PDE3_Ntv_cAMP = PDE3a**pn/(kcm2**pn+PDE3a**pn) 
    
  Jg12= 1*(Kck*kc1*(Gnp**pe/(kcm1**pe+Gnp**pe)))-1*kc2*PDE3_Ntv_cAMP*cAMP
  Va1=0.9
  Va2=8.0
  kcamp1=2*3.2e-6
  kcamp2=4*3.2e-6
  PKAt=0.6e-3
  cAMP_ptv_PKA = cAMP**3/(kcamp1**3+cAMP**3)
  cAMP_ntv_PKA = kcamp2/(kcamp2+cAMP)
  Jg13= 1*(Va1*(cAMP_ptv_PKA)*(PKAt-PKA))-(Va2*PKA*cAMP_ntv_PKA)  
    
  ## dPKA
  PKAprime = Jg13  

  ## dcAMP 
  cAMPprime = Jg12 - 2*PKAprime 

  ## PDE3a
  kakt=1*0.3
  kpde=1*1.50
  PDE3t=5
  PKA_L=PKA/8e-6

  AKT_ptv_PDE3 = kakt*(AKTp**2/(0.1**2+AKTp**2))
  PKA_ntv_PDE3 = (PKA_L**1.0)/(12**1.0+(PKA_L**1.0))
  PDE3prime=(AKT_ptv_PDE3)*(PDE3t-PDE3a) - (kpde*PKA_ntv_PDE3*PDE3a)

  #AMPK
  #(Somvanshi et al., 2016)

  AMP_ATP_AMPK=1
  AKT_Ntv_AMPK=5*(AKTp**2/(AKTp**2+0.05**2))
  PKA_Ntv_AMPK=(3**2/(3**2+PKA_L**2))

  AMPKt=1
  kam1=1
  Kam2=2.25
  AMPKprime=1*(kam1*AMP_ATP_AMPK*(PKA_Ntv_AMPK)*(AMPKt-AMPK))-Kam2*(AMPK)*AKT_Ntv_AMPK       

  return [IRprime, IRpprime, IRSprime, IRSpprime, PI3Kprime, PI3Kpprime, AKTprime, AKTpprime,
          PKCprime, PKCpprime, GLUT4cprime, GLUT4sprime, mTORprime, mTORaprime, S6K1prime, S6K1aprime,
          InsVprime, InsLprime, InsPprime, Glnpprime, Glcnprime, FRprime, RSprime, LRprime, LRSprime, 
          LRPprime, DAGprime, PKAprime, cAMPprime, PDE3prime, AMPKprime]