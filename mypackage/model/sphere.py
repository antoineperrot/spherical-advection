from .model import Model
from ..tools import LatLonGrid
import numpy as np
class SemiLag(Model):
    def __init__(self, n_lat, wind, dt = 1e-1):
        self.grid = LatLonGrid(n_lat=n_lat)
        self.dt = dt
        self.wind = wind
        self.backtrack_points = self.compute_backtracking_points(wind)
        
    def compute_backtracking_points(self, wind):
        wind_x, wind_y, wind_z = wind
        
        # computing departures points in R^3 :
        Xd = self.grid.X - self.dt*wind_x
        Yd = self.grid.Y - self.dt*wind_y
        Zd = self.grid.Z - self.dt*wind_z
        
        # projecting them onto the sphere (approximation) :
        norm_d = np.linalg.norm(np.array([Xd,Yd,Zd]),axis=0)
        Xd = Xd/norm_d
        Yd = Yd/norm_d
        Zd = Zd/norm_d
        
        #converting them to (theta,phi) system of coordinates :
        phid = np.arctan2(Yd,Xd)
        thetad = np.arcsin( Zd/(Xd**2+Yd**2+Zd**2)**.5)

        #getting their indices :
        PHId = ((phid + np.pi)/(2*np.pi) * (self.grid.m))
        THETAd = (thetad + np.pi/2) / (np.pi/(self.grid.n+1))
        
        return (PHId, THETAd)
    
    def _step(self, x, t, new_dt = None):
        dt_ =  new_dt if new_dt else self.dt
        x = self.grid.BilinearSphericalInterpolation(self.backtrack_points, x)
        t += dt_
        
        return x, t    
