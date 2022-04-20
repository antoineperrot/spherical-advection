from ..model import Model, Scheme
from .multivariate_1d_enkf import *
import numpy as np

from ..tools.domain import Domain1D


from scipy.interpolate import interp1d

k3_values =[0,0,0,0,0,0.00675528,.19,.3910734,.5074326,.5755002,.611526,.622824,.622824,.611526,.5755002,.5074326,.3910734,.1972314,.00675528,0,0,0,0,0]

k3_temp = interp1d(range(25),k3_values+ [0],kind='quadratic')

sigma=0.02
def smooth(x):
    return np.maximum(x,0)**2 * np.exp(-x**2/sigma**2)+ np.maximum(x,0)* (1-np.exp(-x**2/sigma**2))

min_to_hour = 60
def k3(t):
    return smooth(k3_temp(t%24.0)) * min_to_hour


# approx for smoothing k3 func :
def k3(t):
    return np.exp( -np.abs((t%24.0)-12)**3/10**2) * 0.6239692319730302 * min_to_hour


def k1(t): return 0.00152 * k3(t)


k2 = 12.3* min_to_hour
k4 = 0.275* min_to_hour
k5 = 10.2* min_to_hour 
k6 = 0.12* min_to_hour



day_to_hour = 1./24

eROC = 0.0235 * day_to_hour
eNO = .243 * day_to_hour
eNO2 = .027 * day_to_hour
l = .02 * day_to_hour


class AdvectionGRS(Model, MultivariateEnKF1D):
    def __init__(self, u,cadastre, dt = 1e-1, scheme = 'RK4', name= 'Undefined model name'):
        Model.__init__(self, dt, scheme, name)
        self.u = u
        self.domain = Domain1D(len(u))
        self.du = self.domain.derive(u)
        sup_u = np.linalg.norm(u,ord=np.inf)
        self.dt = self.domain.dx/sup_u
        
        self.cadastre = cadastre
        self.km = np.arange(-1, self.domain.Nx-1)%self.domain.Nx

        self.kp = np.arange(1, self.domain.Nx+1 )%self.domain.Nx
        MultivariateEnKF1D.__init__(self, ['ROC', 'RP', 'NO', 'NO2', 'O3', 'SN(G)N'], self.domain)
        
    def trend(self, state,t):
        #transport :
        dstate = self.domain.derive(state,axis = 2)
#         dstate = (state[:,:,self.kp] - state[:,:,self.km])/(2*self.domain.dx) 
        trend_transport = -self.u*dstate - state*self.du
        
        
        #chemistry
        ROC, RP, NO, NO2, O3, SNGN = state 
        dtROC =  eROC*self.cadastre - l *ROC
        dtRP = k1(t)*ROC - RP*(k2*NO + 2*k6*NO2 + k5*RP) -l*RP
        dtNO = eNO*self.cadastre +k3(t)*NO2 -NO*(k2*RP + k4*O3) -l*NO
        dtNO2 = eNO2*self.car + k4*NO*O3 +k2*NO*RP - NO2*(k3(t)+2*k6*RP) -l*NO2
        dtO3 = k3(t)*NO2 - k4*NO*O3 -l*O3
        dtSNGN = 2*k6*NO2*RP -l*SNGN     
        trend_chem = np.array([dtROC,dtRP,dtNO,dtNO2,dtO3,dtSNGN])

        #GLOBAL :
        trend_state = trend_transport + trend_chem
        return trend_state
    
#     def _step(self, x, t, new_dt = None): 
#         x, t  = Model._step(self, x, t, new_dt)
#         #clipping negative values :
#         x = (x>0)*x  
#         return x, t

