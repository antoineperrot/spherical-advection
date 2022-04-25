import numpy as np

from matplotlib import colors, cm
import matplotlib.pyplot as plt


class LatLonGrid(object):
    def __init__(self,n_lat, n_long=None,R=1):
        self.n = n_lat
        self.m = 2*n_lat if n_long is None else n_long
        
        self.phi   = np.linspace(-np.pi, np.pi, self.m+1)[:-1]  #phi values
        self.theta = np.linspace(-np.pi/2,np.pi/2,self.n+2)[1:-1]#theta values
        
        self.dtheta = np.abs(self.theta[1]-self.theta[0])
        self.dphi = np.abs(self.phi[1]-self.phi[0])
        
        self.THETA, self.PHI = np.meshgrid(range(self.n),range(self.m))  #indices meshgrid
        
        self.xtheta, self.xphi = np.meshgrid(self.theta, self.phi)  # values meshgrid
        
        self.ex, self.ey, self.ez = np.eye(3)  # basis vectors of R^3
        
        self.X = R * np.cos(self.xtheta)*np.cos(self.xphi)  
        self.Y = R * np.cos(self.xtheta)*np.sin(self.xphi)
        self.Z = R * np.sin(self.xtheta)

        self.OM = self.X[:,:,np.newaxis]*self.ex\
                + self.Y[:,:,np.newaxis]*self.ey\
                + self.Z[:,:,np.newaxis]*self.ez  # all positions vectors in cartesian coordinates
        

    def projection(self, vecteur):
        return (np.dot(vecteur,self.ex),np.dot(vecteur,self.ey),np.dot(vecteur,self.ez))
    
    
   
    def sphere_plot(self, values,figsize=(16,16),cmap=cm.seismic,title='',**kwargs):
        
    
        fig = plt.figure(figsize=figsize)
        fig.suptitle(title, **kwargs)
        ax = fig.add_subplot(1,1,1,projection='3d')

        ax.view_init(elev=20., azim=20)

        strength = values
        norm=colors.Normalize(vmin = 0,
                              vmax = np.max(strength), clip = False)

        surface = ax.plot_surface(*self.projection(self.OM), rstride=1, cstride=1,
                               linewidth=0, antialiased=False,
                               facecolors=cm.seismic(norm(strength)),zorder=1)
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
