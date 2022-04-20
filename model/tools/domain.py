import numpy as np

class Domain1D:
    def __init__(self, Nx):
        self.Nx = Nx
        self.x = np.linspace(0,1,Nx+1)[:-1]
        self.dx = self.x[1]
        self.shape = self.x.shape
        
    def derive(self, u,ord = 1,axis=0):
        if ord == 2:
            return ( np.roll(u,-1,axis=axis) - 2*u + np.roll(u,1,axis=axis) ) / self.dx**2
        else :
            return (np.roll(u,-1,axis=axis) - np.roll(u,1,axis=axis))/2/self.dx


class Domain2D:
    def __init__(self, Nx, Ny):
        self.Nx = Nx
        self.Ny = Ny
        self.dx = 1/Nx
        self.dy = 1/Ny
        self.dxs = [self.dx,self.dy]
        self.x = np.linspace(0,1,Nx)
        self.shape = self.x.shape
    def derive(self, u, axis=0):
        return (np.roll(u,-1,axis=axis) - np.roll(u,1,axis=axis))/self.dxs[axis]/2       
