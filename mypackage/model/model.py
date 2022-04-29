import numpy as np
import matplotlib.pyplot as plt


from tqdm import tqdm

class Scheme:
    "Available schemes : Explicit Euler, RK2, RK4"
    def __init__(self, scheme_name):
        "Available schemes : Explicit Euler, RK2, RK4"
        schemes  = ['EE','RK2','RK4']  #list of all available schemes
        if scheme_name not in schemes :
            raise ValueError("This scheme is not defined.")
        self.scheme_name = scheme_name
        self.scheme = self._make_scheme()

    def _make_scheme(self):
        def EE_scheme(f,x,t,dt):
            x_next = x + dt*f(x,t)
            return x_next
        
        def RK2_scheme(f,x,t,dt):
            k1 = x + dt/2 * f(x,t)
            k2 = f(k1, t + dt/2)
            x_next = x + dt * k2
            return x_next

        def RK4_scheme(f,x,t,dt):
            k1 = f(x,t)
            k2 = f(x + dt/2 * k1, t + dt/2)
            k3 = f(x + dt/2 * k2, t + dt/2)
            k4 = f(x + dt   * k3, t + dt  )
            x_next = x + dt/6 *(k1 + 2*k2 + 2*k3 + k4)
            return x_next
        

        if self.scheme_name =='EE': return EE_scheme
        if self.scheme_name =='RK2': return RK2_scheme
        if self.scheme_name =='RK4': return RK4_scheme


    def step(self,f,x,t,dt):
        return self.scheme(f,x,t,dt)


class Model:
    def __init__(self, dt = 1e-1, scheme = 'EE', name= 'Undefined model name'):
        self.dt = dt
        self.scheme = Scheme(scheme)
        self.name = name
        
    def trend(self, state, t):
        raise NotImplementedError
    
    def _step(self, x, t, new_dt = None):  
        dt_ =  new_dt if new_dt else self.dt
        x = self.scheme.step(self.trend, x, t, dt_)
        t += dt_
        
        return x, t
        
    def forecast(self, x0, t_end, t0=0, time_saving_interval=0.5, show_pbar=True):
        trajectory = {t0:x0}
        last_saved_time = t0
        t = t0
        x = x0
        
        
        if show_pbar : pbar = tqdm(initial=0,total=int(np.ceil((t_end - t)/self.dt)),position=0)
        while t < t_end :
            dt = self.dt if t + self.dt  <= t_end else t_end - t
            x, t = self._step(x, t, dt)
            if t - last_saved_time >= time_saving_interval :
                trajectory[t] = x;
                last_saved_time = t;
            if show_pbar : pbar.update(1)
        
        if show_pbar : pbar.close()
        
        if t != last_saved_time :
            trajectory[t] = x;
            last_saved_time = t;

        #could be usefull :
#        self.x = x
#        self.t = t
        
        return trajectory
