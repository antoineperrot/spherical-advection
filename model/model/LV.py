from .model import Model, Scheme
import numpy as np
import matplotlib.pyplot as plt


class LotkaVolterra(Model):
    def __init__(self, k1, k2, k3, dt = 1e-1, scheme = 'EE', name= 'Undefined model name'):
        self.k1, self.k2, self.k3 = k1, k2, k3
        super().__init__(dt, scheme, name)
       
    
    def trend(self, state, t):
        A,B = state.T
        fA = self.k1*A - self.k2*A*B
        fB = self.k2*A*B - self.k3*B
        trend_state = np.vstack((fA,fB)).T
        return trend_state
    
    def plot(self, times, trajectory,**kwargs):
        x1,x2 = np.array([trajectory[time].mean(axis=0) for time in times]).T
        plt.plot(x1,x2,**kwargs)
        
    def assimilation(self, ensemble_f, y, H, R, Q): #obs, obs operator, obs cov matrix error, model cov matrix error
        Xf = ensemble_f.mean(axis=0)  #forecasts
        Ef = np.cov(ensemble_f.T) +  Q #error forecast
        K = Ef @ H.T @ np.linalg.inv(H.T @ Ef@ H + R)  # kalman gain
        
        Ne = len(ensemble_f)  #ensemble size
        Ea = (np.eye(2) - K @ H) @ Ef
        Xa = Xf + K@(y - H@Xf)
        
        ensemble_a = np.random.multivariate_normal(Xa, Ea, Ne)
        return ensemble_a
        

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
