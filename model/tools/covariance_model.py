import numpy as np

def autocorrelation_model(s):
    Nx = len(s)
    X,Y = np.meshgrid(np.arange(Nx),np.arange(Nx))
    abss = np.linspace(0,1,Nx)
    template_dist = np.abs(abss - 0.5)
    dist = np.array([np.roll(template_dist, r+1) for r in range(-Nx//2,Nx//2)])
    output = (s[X]*s[Y])**.25/(0.5*(s[X]+s[Y]))**.5 * \
        np.exp(- dist**2/ (s[X]+s[Y]))
    return output


def autocovariance_model(V,s):
    Nx = len(s)
    X,Y = np.meshgrid(np.arange(Nx),np.arange(Nx))
    abss = np.linspace(0,1,Nx)
    template_dist = np.abs(abss - 0.5)
    dist = np.array([np.roll(template_dist, r+1) for r in range(-Nx//2,Nx//2)])
    output = (V[X]*V[Y])**.5 *\
        (s[X]*s[Y])**.25/(0.5*(s[X]+s[Y]))**.5 * \
        np.exp(- dist**2/ (s[X]+s[Y]))
    return output

def crosscorrelation_model(Vab,Va,Vb,sA,sB):
    Nx = len(Vab)
    X,Y = np.meshgrid(np.arange(Nx),np.arange(Nx))
    abss = np.linspace(0,1,Nx)
    template_dist = np.abs(abss - 0.5)
    dist = np.array([np.roll(template_dist, r+1) for r in range(-Nx//2,Nx//2)])
    output = 0.5 * (Vab[X]/(Va[X]*Vb[X])**.5 + Vab[Y]/(Va[Y]*Vb[Y])**.5) * \
            np.exp(-dist**2/( 0.5*(sA[X]+sA[Y] +sB[X]+sB[Y]) ))
    return output

    

def crosscovariance_model(Vab,Va,Vb,sA,sB):
    Nx = len(Vab)
    X,Y = np.meshgrid(np.arange(Nx),np.arange(Nx))
    abss = np.linspace(0,1,Nx)
    template_dist = np.abs(abss - 0.5)
    dist = np.array([np.roll(template_dist, r+1) for r in range(-Nx//2,Nx//2)])
    output = 0.5 * (Vab[X]/(Va[X]*Vb[X])**.5 + Vab[Y]/(Va[Y]*Vb[Y])**.5) * \
            np.exp(-dist**2/( 0.5*(sA[X]+sA[Y] +sB[X]+sB[Y]) )) * \
            np.sqrt(Va[X]*Vb[Y])
    return output
