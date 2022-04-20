import numpy as np


def make_samples(ensemble_size, domain, lh=.06):
    gauss = lambda dx : np.exp(-0.5*dx**2/lh**2) # lh has been previously specified 
    correlation = gauss(domain-domain[domain.shape[0]//2])
    spectrum = np.abs(np.fft.fft(correlation))

    std_spectrum= np.sqrt(spectrum)
    samples = []
    for i in range(ensemble_size):
        zeta = np.random.normal(size=domain.shape)
        zeta = np.fft.fft(zeta)
        ef = np.fft.ifft(std_spectrum * zeta)
        ef = np.real(ef)
        samples.append(ef)
    samples = np.vstack(samples)
    return samples



def make_bivariate_ensemble(size, initial_setting, xx):
    # defining initial state :
    A0 = np.ones(len(xx)) * initial_setting['A']
    B0 = np.ones(len(xx)) * initial_setting['B']

    # generating initial error fields :
    perturbationsA = make_samples(size, xx, initial_setting['Length-scale A'])\
                        * initial_setting['std A']
    perturbationsB = make_samples(size, xx, initial_setting['Length-scale B'])\
                        * initial_setting['std B']

    As = A0 + perturbationsA
    Bs = B0 + perturbationsB
    ensemble0 = np.swapaxes(np.array([As,Bs]),0,1)
    return ensemble0




import pyshtools as sh  #py spherical harmonics tools package https://shtools.github.io/SHTOOLS/index.html
from scipy.special import factorial

def random_isotropic_field(m,n, Lmax=20,alpha=0.7):
    phi,theta = np.linspace(-np.pi,np.pi,m), np.linspace(-np.pi/2,np.pi/2,n)
    ls = [ [index] * (index+1) for index in range(0,Lmax+1)]
    ms=  [  list(range(0,l+1)) for l in range(0,Lmax+1)  ]

    ls = np.array(np.sum(ls))
    ms = np.array(np.sum(ms))

    Llm =    (  (factorial(ls-ms)/factorial(ls+ms) * (2*ls+1)/(4*np.pi))**.5  ) * \
        np.array([sh.legendre.legendre_lm(ls,ms, np.sin(_theta),normalization='unnorm') for _theta in theta]) 

    # influence the length-scales by putting more or less weigth on the spherical harmonics of higher order.
    Al = .7**np.arange(0,Lmax+1,1)

    Ll0_index = [n*(n+1)//2 for n in range(Lmax+1)]

    T = np.zeros((m,n))*1.0
    for l, Al_coef in enumerate(Al):
        T += Al_coef**.5 * Llm.T[Ll0_index[l]] * np.random.normal(0,1)

        for m in range(1,l):
            mm  = Ll0_index[l] + m
            T += (2*Al_coef)**.5 * Llm.T[mm] * (np.random.normal(0,1) * np.cos(m*phi)[:,np.newaxis] + 
                                              np.random.normal(0,1) * np.sin(m*phi)[:,np.newaxis])
    return T


import pyshtools as sh  #py spherical harmonics tools package https://shtools.github.io/SHTOOLS/index.html
from scipy.special import factorial

def random_isotropic_ensemble(ensemble_size,m,n,Lmax=20, alpha=0.7):
    phi,theta = np.linspace(-np.pi, np.pi - 2*np.pi / m, m), np.linspace(-np.pi/2,np.pi/2,n)
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

        for k in range(1,l):
            mm  = Ll0_index[l] + k
            ensemble += (2*Al_coef)**.5 * Llm.T[mm] * (np.random.normal(0,1,ensemble_size)[:,np.newaxis,np.newaxis] * np.cos(k*phi)[:,np.newaxis] + 
                                              np.random.normal(0,1,ensemble_size)[:,np.newaxis,np.newaxis] * np.sin(k*phi)[:,np.newaxis])

    ensemble = ensemble - ensemble.mean(axis=0)

    ensemble = ensemble / ensemble.std(axis=0)
    
    return ensemble


def correlation_matrix(ensemble):
    tmp = np.hstack([ensemble[i] for i in range(ensemble.shape[0])]).T
    return np.corrcoef(tmp)

def covariance_matrix(ensemble):
    tmp = np.hstack([ensemble[i] for i in range(ensemble.shape[0])]).T
    return np.cov(tmp)
    
def get_correlation_matrix(ensemble, species_name):
    Nx = ensemble.shape[2]
    ensemble = np.hstack([*ensemble])
    cm = np.corrcoef(ensemble.T)
    cor_matrix={}
    for i,field in enumerate(species_name):
        for j,field2 in enumerate(species_name):
            cor_matrix[field+'/'+field2] = cm[i*Nx:(i+1)*Nx,j*Nx:(j+1)*Nx]
    return cor_matrix


