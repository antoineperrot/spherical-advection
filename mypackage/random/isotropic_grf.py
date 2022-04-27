import numpy as np
import pyshtools as sh  #py spherical harmonics tools package https://shtools.github.io/SHTOOLS/index.html
from scipy.special import factorial

def random_isotropic_ensemble(ensemble_size,m,n,Lmax=20, alpha=0.7):
    phi,theta = np.linspace(-np.pi, np.pi, m+1)[:-1], np.linspace(-np.pi/2,np.pi/2,n+2)[1:-1]
    
    #phi,theta = np.linspace(-np.pi,np.pi,m), np.linspace(-np.pi/2,np.pi/2,n)
    ls = [ [index] * (index+1) for index in range(0,Lmax+1)]
    ms=  [  list(range(0,l+1)) for l in range(0,Lmax+1)  ]

    ls = np.array(np.sum(ls))
    ms = np.array(np.sum(ms))

    Llm =    (  (factorial(ls-ms)/factorial(ls+ms) * (2*ls+1)/(4*np.pi))**.5  ) * \
        np.array([sh.legendre.legendre_lm(ls,ms, np.sin(_theta),normalization='unnorm') for _theta in theta]) 

    # influence the length-scales by putting more or less weigth on the spherical harmonics of higher order.
    Al = alpha**np.arange(0,Lmax+1,1)

    Ll0_index = [i*(i+1)//2 for i in range(Lmax+1)]

    ensemble = np.zeros((ensemble_size,m,n))
    for l, Al_coef in enumerate(Al):
        ensemble += Al_coef**.5 * Llm.T[Ll0_index[l]] * np.random.normal(0,1,ensemble_size)[:,np.newaxis,np.newaxis]

        for k in range(1,l+1):
            mm  = Ll0_index[l] + k
            ensemble += (2*Al_coef)**.5 * Llm.T[mm] * (np.random.normal(0,1,ensemble_size)[:,np.newaxis,np.newaxis] * np.cos(k*phi)[:,np.newaxis] + 
                                              np.random.normal(0,1,ensemble_size)[:,np.newaxis,np.newaxis] * np.sin(k*phi)[:,np.newaxis])
    
    return ensemble
