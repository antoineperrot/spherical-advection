from ..tools.domain import Domain1D

from ..model import Model, Scheme

from .multivariate_1d_pkf import *
from .multivariate_1d_enkf import *
import numpy as np

from ..tools.domain import Domain1D

class OneDimLV(Model, MultivariateEnKF1D):
    def __init__(self, k1, k2, k3, Nx, activate_closure = False, dt = 1e-1, scheme = 'EE'):
        Model.__init__(self, dt, scheme)
        self.k1,self.k2,self.k3 = k1,k2,k3
        self.activate_closure = activate_closure
        self.domain = Domain1D(Nx)
        MultivariateEnKF1D.__init__(self, ['A','B'], self.domain)
    def trend(self, state,t):
        A,B = state
        trend_A = + self.k1*A - self.k2*A*B
        trend_B = - self.k3*B + self.k2*A*B
        trend_state = np.array([trend_A, trend_B])
        return trend_state

class PKFOneDimLV(Model, MultivariatePKF1D):
    def __init__(self, k1, k2, k3,Nx, activate_closure = False, dt = 1e-1, scheme = 'EE',name='None'):
        Model.__init__(self, dt, scheme, name)
        self.k1,self.k2,self.k3 = k1,k2,k3
        self.activate_closure = activate_closure*1.0
        self.domain = Domain1D(Nx)
        MultivariatePKF1D.__init__(self, ['A','B'], self.domain)

        
    def trend(self, state, t):
        A, B, Va, Vb, sa, sb, Vab = state
        dA, dB, dVa, dVb, dsa, dsb, dVab = self.domain.derive(state,axis=1)

        fA = + self.k1*A - self.k2*A*B - self.k2*Vab
        fB =  - self.k3*B + self.k2*A*B + self.k2*Vab
        fVab = Vab*(self.k1-self.k2*B-self.k3+self.k2*A) + self.k2*Va*B - self.k2*Vb*A 
        fVa = + 2* (Va * (self.k1 -self.k2*B) - self.k2*A*Vab ) 
        fVb = + 2* (Vb * (-self.k3+self.k2*A) + self.k2*B*Vab ) 

        closure1 = 2*Vab/(Va*Vb)**.5 / (sa+sb)
        closure2_3 = 0

        fsa = self.activate_closure*(- 2*self.k2*A*Vab*sa/Va + 2*self.k2*A*Vb**.5*sa**2*closure1/Va**.5 ) #+ k2**A*sa**2*closure2_3*dVb/(Va*Vb)**.5    - k2*A*Vb**.5*sa**2*closure2_3*dVa/Va**1.5 + 2*k2*Vb**.5*sa**2*closure2_3*dA/Va**.5
        fsb = self.activate_closure*(  2*self.k2*B*Vab*sb/Vb - 2*self.k2*B*Va**.5*sb**2*closure1/Vb**.5 )  #+ k2*B*Va**.5*sb**2*closure2_3*dVb /Vb**1.5 - k2*B*sb**2*closure2_3*dVa/(Va*Vb)**.5    - 2*k2*Va**.5*sb**2*closure2_3*dB/Vb**.5

        trend_state = np.array([fA,fB,fVa,fVb,fsa,fsb, fVab])
        return trend_state
