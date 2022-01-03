# Criticality function using nuFiss, Absorption, geometric buckling

def calc_keff_buck(nuFiss, absorb, diffusionCoe, geometryPara):
    import numpy as np
    k_inf = nuFiss/absorb
    L2 = diffusionCoe/absorb
    B2 = (np.pi/geometryPara)**2

    return k_inf/(1+L2*B2)