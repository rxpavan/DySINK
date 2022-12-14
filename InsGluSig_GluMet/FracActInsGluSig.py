import numpy as np
def FracActInsGluSig(ConSignal):
  # Given the concentration of the signalling molecules this function returns
  # the fractional activation values of the signals that effect
  # metabolic network
  
  Cmax_akt = 10
  Cmax_pka = 6.0e-4

  akt_frac = ConSignal[7]/Cmax_akt
  pka_frac = ConSignal[27]/Cmax_pka
  glut_frac = ConSignal[11]/(ConSignal[10]+ConSignal[11]);
  return np.array([akt_frac, pka_frac,glut_frac])