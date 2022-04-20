


from .model import Model, Scheme
import numpy as np

## draw covariance matrix ?
import matplotlib.pyplot as plt

class LotkaVolterra_PKF(Model):
    def __init__(self, k1, k2, k3, dt = 1e-1, scheme = 'EE', name= 'Undefined model name'):
        self.k1, self.k2, self.k3 = k1, k2, k3
        super().__init__(dt, scheme, name)

        
    def trend(self, state, t):
        A, B, Va, Vb, Vab = state

        fA   =  self.k1*A - self.k2*A*B - self.k2*Vab
        fB   = -self.k3*B + self.k2*A*B + self.k2*Vab
        fVa  = 2* ( self.k1*Va - self.k2*Vab*A - self.k2*B*Va)
        fVb  = 2* (-self.k3*Vb + self.k2*Vab*B + self.k2*A*Vb)
        fVab = Vab*(self.k1-self.k2*B-self.k3+self.k2*A) +self.k2*Va*B -self.k2*Vb*A

        trend_state = np.array([fA,fB,fVa,fVb,fVab])
        return trend_state
    
    def plot(self, times, trajectory,**kwargs):
        A,B = np.array([trajectory[time][:2] for time in times]).T
        plt.plot(A,B,**kwargs)

            
            
def draw_covariance_matrix(mean, covariance_matrix, radius=1,**args):
    """
    Draw the ellispe corresponding to the covariance matrix M at location the mean and of radius radius
    """
    Lambda, P = np.linalg.eigh(covariance_matrix)
    theta = np.linspace(0,2*np.pi,200)
    y1 = np.sqrt(Lambda[0]) * np.cos(theta)
    y2 = np.sqrt(Lambda[1]) * np.sin(theta)
    y = np.array([y1,y2])
    x1, x2 = np.dot(P,y) * radius
    x1 += mean[0]
    x2 += mean[1]
    plt.plot(x1,x2,**args)
    plt.scatter(mean[0],mean[1],marker='o',zorder=3,**args)
    return
