from ..model import Model, Scheme
from .multivariate_1d_enkf import *
import numpy as np

from ..tools.domain import Domain1D


from model import Model
class AdvectionLV(Model, MultivariateEnKF1D):
    def __init__(self, u, k1, k2, k3, dt = 1e-1, scheme = 'EE'):
        Model.__init__(self, dt, scheme)
        self.u,self.k1,self.k2,self.k3 = u,k1,k2,k3
        self.domain = Domain1D(len(u))
        self.du = self.domain.derive(u) 
        sup_u = np.linalg.norm(u,ord=np.inf)
        self.dt = self.domain.dx/sup_u
        
        MultivariateEnKF1D.__init__(self, ['A','B'], self.domain)

    def trend(self, state,t):
        A,B = state
        dA, dB = self.domain.derive(state,axis=2)
        trend_A = -self.u*dA - A*self.du + self.k1*A - self.k2*A*B
        trend_B = -self.u*dB - B*self.du - self.k3*B + self.k2*A*B
        trend_state = np.array([trend_A, trend_B])
        return trend_state


