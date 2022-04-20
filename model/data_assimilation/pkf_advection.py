from .multivariate_1d_pkf import MultivariatePKF1D
from ..model import Model

from ..tools.domain import Domain1D

import numpy as np

class PKFAdvection(Model, MultivariatePKF1D):
    def __init__(self, u, dt = 1e-1, scheme = 'EE'):
        Model.__init__(self, dt, scheme)
        self.u = u
        self.domain = Domain1D(len(u))
        self.du = self.domain.derive(u) 
        sup_u = np.linalg.norm(u,ord=np.inf)
        self.dt = self.domain.dx/sup_u
        
        MultivariatePKF1D.__init__(self, ['X'], self.domain)

        
    def trend(self, state,t):
        A,V,s = state
        dA, dV, ds = self.domain.derive(state,axis=1)
        trend_A = -self.u*dA - A * self.du
        trend_V = -self.u*dV -2*V *self.du
        trend_s = -self.u*ds +2*s *self.du
        trend_state = np.array([trend_A,trend_V,trend_s])
        return trend_state    
