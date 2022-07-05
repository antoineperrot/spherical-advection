import numpy as np

from matplotlib import colors, cm
import matplotlib.pyplot as plt

from ..random import *

class LatLonGrid(object):
    def __init__(self,n_lat, n_long=None,R=1):
        self.n = n_lat
        self.m = 2*n_lat if n_long is None else n_long
        
        self.phi   = np.linspace(-np.pi, np.pi, self.m+1)[:-1]  #phi values
        self.theta = np.linspace(-np.pi/2,np.pi/2,self.n+2)[1:-1]#theta values
        self.R = R
        
        self.dtheta = np.abs(self.theta[1]-self.theta[0])
        self.dphi = np.abs(self.phi[1]-self.phi[0])
        
        self.THETA, self.PHI = np.meshgrid(range(self.n),range(self.m))  #indices meshgrid
        
        self.xtheta, self.xphi = np.meshgrid(self.theta, self.phi)  # values meshgrid
        
        self.ex, self.ey, self.ez = np.eye(3)  # basis vectors of R^3
        
        self.e_r = np.array([np.cos(self.xtheta)*np.cos(self.xphi),
               np.cos(self.xtheta)*np.sin(self.xphi),
               np.sin(self.xtheta)]) * self.R
        self.e_theta = np.array([-np.sin(self.xtheta)*np.cos(self.xphi),
                       -np.sin(self.xtheta)*np.sin(self.xphi),
                       np.cos(self.xtheta)]) * self.R

        self.e_phi = np.array([-np.cos(self.xtheta)*np.sin(self.xphi),
                         np.cos(self.xphi)*np.cos(self.xtheta),
                         np.cos(self.xtheta)*0]) * self.R
        
        self.X = self.R * np.cos(self.xtheta)*np.cos(self.xphi)  
        self.Y = self.R * np.cos(self.xtheta)*np.sin(self.xphi)
        self.Z = self.R * np.sin(self.xtheta)

        self.OM = self.X[:,:,np.newaxis]*self.ex\
                + self.Y[:,:,np.newaxis]*self.ey\
                + self.Z[:,:,np.newaxis]*self.ez  # all positions vectors in cartesian coordinates
        

    def projection(self, vecteur):
        return (np.dot(vecteur,self.ex),np.dot(vecteur,self.ey),np.dot(vecteur,self.ez))
    
    
   
    def sphere_plot(self, values,figsize=(16,16),cmap=cm.seismic,elev=20,azim=20,title='',**kwargs):
        
    
        fig = plt.figure(figsize=figsize)
        fig.suptitle(title, **kwargs)
        ax = fig.add_subplot(1,1,1,projection='3d')

        ax.view_init(elev=elev, azim=azim)

        strength = values
        norm=colors.Normalize(vmin = np.min(strength),
                              vmax = np.max(strength), clip = False)

        surface = ax.plot_surface(*self.projection(self.OM), rstride=1, cstride=1,
                               linewidth=0, antialiased=False,
                               facecolors=cmap(norm(strength)),zorder=1)
        return fig;
    
    def regular_points_generator_for_sphere(self, k):
    #thetas 
        x1s = (np.linspace(0,self.n-1, k)).astype(int)

        x1s_x = self.theta[x1s]

        #phis
        x2s = []
        for x1,x1_x in zip(x1s, x1s_x):
            k_ = int(np.cos(x1_x)* 2 * (k))
            vals = list((np.linspace(0,self.m-1, k_)[:-1]).astype(int))
            if k_ ==0 : vals= [0]

            x2s.append(vals)

        list_pos = []
        for i,pos0 in enumerate(x1s):
            for pos1 in x2s[i]:
                list_pos.append((pos1,pos0))

        return list_pos
    
    def dx_frompoint(self, tuple_index_pos):
        '''
        return a the grid of (d_theta, d_phi) value for a referent point of the grid, which indexes
        are indicated in tuple_index_pos
        '''
        Mstar = self.OM[tuple_index_pos]

        eR = np.cross(Mstar,self.ex) # rotation axis
        eR /= np.linalg.norm(eR)
        theta_star = -np.arccos(np.vdot(Mstar,self.ex))  # rotation angle

        e2 = np.cross(self.ex, eR)
        e2 /= np.linalg.norm(e2)

        P = np.array([self.ex, eR, e2]).T  # change of basis matrix from canonical basis to B' = {ex, eR, e2}

        R_BpBp = np.array([[ np.cos(theta_star), 0, np.sin(theta_star)], # R_PP of ex in B' coordinates
                         [0                  , 1, 0                 ], # R_PP of eR in B' coordinates
                         [-np.sin(theta_star), 0, np.cos(theta_star)]]) # R_PP of e2 in B' coordinates

        R = P @ R_BpBp @ np.linalg.inv(P)  # Rotation matrix in canonical basis

        # # verif ( equal [0,0,0]): 
        # np.round(R @ ex - Mstar, 3), np.round(  R @ eR - eR,3)


        outR = (self.OM @ R )

        # outR = outR / np.linalg.norm(outR, axis=2)[:,:,np.newaxis]

        xs = outR[:,:,0]; ys = outR[:,:,1] ; zs = outR[:,:,2]

        zs[zs >1] = 1
        zs[zs <-1] = -1

        dx1 = np.arccos(zs) - np.pi/2 # dtheta
        dx2 = np.arctan2(ys, xs)

        dx = np.stack([dx1,dx2],axis=2)
        
        return dx
    
    def isotropic_grf(self, ensemble_size=1, Lmax=20, alpha=0.7):
        return random_isotropic_ensemble(ensemble_size,
                                        self.m,self.n,
                                        Lmax=Lmax, alpha=alpha)

    def ensemble_diagnosis(self, ensemble, mode='aspect'):
        mean = ensemble.mean(axis=0)
        std = ensemble.std(axis=0)

        errors = (ensemble - mean)/std

        d_phi_errors =  ((np.roll(errors,-1,axis=1) - np.roll(errors,1,axis=1))/2/self.dphi ) / np.cos(self.xtheta)
        d_theta_errors =  ((np.roll(errors,-1,axis=2) - np.roll(errors,1,axis=2))/2/self.dtheta )
        del errors
        d_theta_errors[:,:,[0,self.n-1]] = 0
        d_phi_errors[:,:,[0,self.n-1]] = 0

        metric_tensor = np.zeros((2,2,mean.shape[0],mean.shape[1]))

        metric_tensor[0,0] = (d_theta_errors**2).mean(axis=0)
        metric_tensor[1,1] = (d_phi_errors**2).mean(axis=0)
        metric_tensor[0,1] = (d_phi_errors*d_theta_errors).mean(axis=0)
        metric_tensor[1,0] = (d_phi_errors*d_theta_errors).mean(axis=0)

        aspect_tensor =1/(metric_tensor[0,0] * metric_tensor[1,1] - metric_tensor[0,1]**2) * np.array([[metric_tensor[1,1],-metric_tensor[0,1]],
             [-metric_tensor[0,1], metric_tensor[0,0]]])
        aspect_tensor[:,:,:,[0,-1]] = 0

        if mode =='aspect':
            return mean, std, aspect_tensor[:,:,:,1:-1]
        else :
            return mean, std, metric_tensor[:,:,:,1:-1]


    def BilinearSphericalInterpolation(self, backtrack_points, state):
        PHId, THETAd = backtrack_points
        
        south_pole_mean = state[:,:,[0]]*0 + state[:,:,[0]].mean(axis=1)[:,np.newaxis]
        north_pole_mean = state[:,:,[-1]]*0 + state[:,:,[-1]].mean(axis=1)[:,np.newaxis]

        poles_extended_state = np.dstack([south_pole_mean, state, north_pole_mean])
    
     # bilinear interpolation on "rectangles" defined by the lat/lon grid.
        PHI_ = PHId.astype(int)
        PHIp_ = (PHId+1).astype(int) % self.m  # 
    
        THETA_ = THETAd.astype(int)
        THETAp_ = (THETAd+1).astype(int)
    
        interpolation = (1-THETAd%1)         *(  (1-PHId%1) * poles_extended_state[:,PHI_      , THETA_        ]
                                           + (PHId % 1) * poles_extended_state[:,PHIp_ , THETA_        ] ) \
                +  (THETAd % 1)          *(  (1-PHId%1) * poles_extended_state[:,PHI_      , THETAp_  ]
                                           + (PHId % 1) * poles_extended_state[:,PHIp_ , THETAp_  ] )
    
        return interpolation

        


