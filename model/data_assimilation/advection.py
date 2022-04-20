from ..tools.domain import Domain1D

from ..model import Model

import numpy as np
from .multivariate_1d_enkf import *
#domain ?

class Advection(Model, MultivariateEnKF1D):
    def __init__(self, u, dt = 1e-1, scheme = 'EE'):
        Model.__init__(self, dt, scheme)
        self.domain = Domain1D(len(u))
        self.u = u
        self.domain = Domain1D(len(u))
        self.du = self.domain.derive(u) 
        sup_u = np.linalg.norm(u,ord=np.inf)
        self.dt = self.domain.dx/sup_u
        
        MultivariateEnKF1D.__init__(self, ['X'], self.domain)
        
    def trend(self, state,t):
        dx_state = self.domain.derive(state,axis=2)
        trend_state = -self.u * dx_state - state*self.du
        return trend_state

