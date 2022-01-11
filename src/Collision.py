# Collision using Monte Carlo method

import numpy as np
import matplotlib.pyplot as plt

def slab_transmission(sigma_t,thickness,N,isotropic=False):
    """
    sigma_t total cross section
    thickness of slab
    N number of neutrons

    solve the fraction of neutrons that make it through
    """

    theta=np.random.random(N)
    if isotropic==False:
        mu=np.ones(N)
    else:
        mu=np.random.random(N)
    x = -np.log(1-theta)/sigma_t
    transmission = np.sum(x>thickness/mu)/N

    if (N<=10000):
        plt.scatter(x*mu,np.arange(N))
        plt.xlabel("distance to collision")
        plt.ylabel("Neutron number")
        plt.show()
    return transmission


if __name__=="__main__":
    print(slab_transmission(2.0,3.0,int(1e7),True))

