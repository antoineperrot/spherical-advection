from .multivariate_1d_pkf import MultivariatePKF1D
from ..model import Model

from ..tools.domain import Domain1D
import numpy as np

from scipy.interpolate import interp1d
min_to_hour = 60 
k3_values =np.array([0,0,0,0,0,0.00675528,.19,.3910734,.5074326,.5755002,.611526,.622824,.622824,.611526,.5755002,.5074326,.3910734,.1972314,.00675528,0,0,0,0,0,0])*min_to_hour

k3_temp = interp1d(range(25),k3_values+ [0],kind='quadratic')

sigma=0.02
# def smooth(x):
#     return np.maximum(x,0)**2 * np.exp(-x**2/sigma**2)+ np.maximum(x,0)* (1-np.exp(-x**2/sigma**2))


# def k3(t):
#     return smooth(k3_temp(t%24.0))


# approx for smoothing k3 func :
def k3(t):
    return np.exp( -np.abs((t%24.0)-12)**3/10**2) * 0.6239692319730302 * min_to_hour

def k1(t): return 0.00152 * k3(t)




k2 = 12.3* min_to_hour
k4 = 0.275* min_to_hour
k5 = 10.2* min_to_hour
k6 = 0.12* min_to_hour

k_2 = 12.3* min_to_hour
k_4 = 0.275* min_to_hour
k_5 = 10.2* min_to_hour
k_6 = 0.12* min_to_hour



day_to_hour = 1./24

eROC = 0.0235 * day_to_hour 
eNO = .243 * day_to_hour
eNO2 = .027 * day_to_hour
l = .02 * day_to_hour

class PKF_Advection_GRS(Model, MultivariatePKF1D):
    def __init__(self, u,cadastre, dt = 1e-1, scheme = 'RK4', name= 'Undefined model name'):
        Model.__init__(self, dt, scheme, name)
        self.u = u
        self.domain = Domain1D(len(u))
        self.du = self.domain.derive(u)
        sup_u = np.linalg.norm(u,ord=np.inf)
        self.dt = self.domain.dx/sup_u
        self.cadastre = cadastre
        MultivariatePKF1D.__init__(self, ['ROC', 'RP', 'NO', 'NO2', 'O3', 'SN(G)N'], self.domain)
        
        #rosenbrock constant :
        self.gamma = 1 + 1/np.sqrt(2)
        
        self.n_fields_ros = int(2*self.n_species + self.n_species*(self.n_species-1)/2)

    def forecast_implicit(self, x0, t_end, t0=0, time_saving_interval=0.5, show_pbar=True):
        trajectory = {t0:x0}
        last_saved_time = t0
        t = t0
        x = x0
        
        
        if show_pbar : pbar = tqdm(initial=0,total=int(np.ceil((t_end - t)/self.dt)),position=0)
        while t < t_end :
            dt = self.dt if t + self.dt  <= t_end else t_end - t
            x, t = self._step_implicit(x, t, dt)
            if t - last_saved_time >= time_saving_interval :
                trajectory[t] = x;
                last_saved_time = t;
            if show_pbar : pbar.update(1)
        
        if show_pbar : pbar.close()
        
        if t != last_saved_time :
            trajectory[t] = x;
            last_saved_time = t;

        #could be usefull :
        self.x = x
        self.t = t
        
        return trajectory
    
    def _step_ros2(self, state, t, dt =0.1):
        
        
        jacs = np.array([self._jacobian(state[:,i][:self.n_fields_ros], t) for i in range(self.domain.Nx) ])
        
        evals_for_k1 = self.chem_trend(state,t)
        k1s = np.array([np.linalg.solve(np.eye(self.n_fields_ros) - self.gamma * dt * jacs[i], evals_for_k1[:,i]) for i in range(self.domain.Nx)]).T
        
        evals_for_k2 = self.chem_trend(state + dt * evals_for_k1, t) - k1s
        k2s = np.array([np.linalg.solve(
            np.eye(self.n_fields_ros) - self.gamma * dt * jacs[i], evals_for_k2[:,i]) for i in range(self.domain.Nx)]).T
        
        out = state + 1.5*dt*k1s + 0.5*dt*k2s
        return out
    
    def _step_implicit(self, x, t, new_dt = None):  
        dt_ =  new_dt if new_dt else self.dt
        
        x_transprt = self.scheme.step(self.transport_trend, x, t, dt_)
        
        x_transprt_chem = x_transprt
        x_transprt_chem[:self.n_fields_ros] += self._step_ros2(x_transprt[:self.n_fields_ros], t ,dt_)
        
        t += dt_
        
        return x_transprt_chem, t
        
    def trend(self, state,t):
        
        out = self.transport_trend(state, t)
        out[:27] += self.chem_trend(state[:27], t)
        return out
    
    def transport_trend(self, state,t):
        #transport :
        dstate = self.domain.derive(state,axis = 1)
        trend_transport_temp = -self.u*dstate 

        state_mean = state[:6]
        state_var = state[6:12]
        state_aspect = state[27:]
        state_covar = state[12:27]

        trend_transport_mean     = -  self.du*state_mean
        trend_transport_var = -2*self.du*state_var
        trend_transport_aspect_tensor = 2*self.du*state_aspect
        trend_transport_covar = -2*self.du*state_covar

        out = trend_transport_temp + np.vstack((trend_transport_mean,
                                                trend_transport_var,
                                                trend_transport_covar,
                                                   trend_transport_aspect_tensor))

        
        return out
    
    
    def _jacobian(self, x, t):
        _Dummy_15162, _Dummy_15161, _Dummy_15165, _Dummy_15164, _Dummy_15163, _Dummy_15160, _Dummy_15141, _Dummy_15140, _Dummy_15143, _Dummy_15144, _Dummy_15142, _Dummy_15139, _Dummy_15150, _Dummy_15153, _Dummy_15152, _Dummy_15151, _Dummy_15149, _Dummy_15148, _Dummy_15147, _Dummy_15146, _Dummy_15145, _Dummy_15157, _Dummy_15156, _Dummy_15155, _Dummy_15159, _Dummy_15158, _Dummy_15154 = x
        k_3 = k3(t)
        k_1 = k1(t)
        return np.array([[-l, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [k_1, -2*_Dummy_15161*k_5 - 2*_Dummy_15164*k_6 - _Dummy_15165*k_2 - l, -_Dummy_15161*k_2, -2*_Dummy_15161*k_6, 0, 0, 0, -k_5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -k_2, -2*k_6, 0, 0, 0, 0, 0, 0, 0, 0], [0, -_Dummy_15165*k_2, -_Dummy_15161*k_2 - _Dummy_15163*k_4 - l, k_3, -_Dummy_15165*k_4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -k_2, 0, 0, 0, 0, -k_4, 0, 0, 0, 0], [0, -2*_Dummy_15164*k_6 + _Dummy_15165*k_2, _Dummy_15161*k_2 + _Dummy_15163*k_4, -2*_Dummy_15161*k_6 - k_3 - l, _Dummy_15165*k_4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, k_2, -2*k_6, 0, 0, 0, k_4, 0, 0, 0, 0], [0, 0, -_Dummy_15163*k_4, k_3, -_Dummy_15165*k_4 - l, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -k_4, 0, 0, 0, 0], [0, 2*_Dummy_15164*k_6, 0, 2*_Dummy_15161*k_6, 0, -l, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*k_6, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -2*l, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -4*_Dummy_15140*k_5 - 4*_Dummy_15147*k_6 - 2*_Dummy_15148*k_2, -2*_Dummy_15140*k_2, -4*_Dummy_15140*k_6, 0, 0, 0, -4*_Dummy_15161*k_5 - 4*_Dummy_15164*k_6 - 2*_Dummy_15165*k_2 - 2*l, 0, 0, 0, 0, 2*k_1, 0, 0, 0, 0, -2*_Dummy_15161*k_2, -4*_Dummy_15161*k_6, 0, 0, 0, 0, 0, 0, 0, 0], [0, -2*_Dummy_15143*k_2, -2*_Dummy_15148*k_2 - 2*_Dummy_15156*k_4, 0, -2*_Dummy_15143*k_4, 0, 0, 0, -2*_Dummy_15161*k_2 - 2*_Dummy_15163*k_4 - 2*l, 0, 0, 0, 0, 0, 0, 0, 0, -2*_Dummy_15165*k_2, 0, 0, 0, 2*k_3, -2*_Dummy_15165*k_4, 0, 0, 0, 0], [0, -4*_Dummy_15144*k_6 + 2*_Dummy_15157*k_2, 2*_Dummy_15147*k_2 + 2*_Dummy_15159*k_4, -4*_Dummy_15147*k_6, 2*_Dummy_15157*k_4, 0, 0, 0, 0, -4*_Dummy_15161*k_6 - 2*k_3 - 2*l, 0, 0, 0, 0, 0, 0, 0, 0, -4*_Dummy_15164*k_6 + 2*_Dummy_15165*k_2, 0, 0, 2*_Dummy_15161*k_2 + 2*_Dummy_15163*k_4, 0, 0, 2*_Dummy_15165*k_4, 0, 0], [0, 0, -2*_Dummy_15142*k_4, 0, -2*_Dummy_15156*k_4, 0, 0, 0, 0, 0, -2*_Dummy_15165*k_4 - 2*l, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*_Dummy_15163*k_4, 0, 2*k_3, 0, 0], [0, 4*_Dummy_15158*k_6, 0, 4*_Dummy_15145*k_6, 0, 0, 0, 0, 0, 0, 0, -2*l, 0, 0, 0, 0, 0, 0, 0, 0, 4*_Dummy_15164*k_6, 0, 0, 0, 0, 4*_Dummy_15161*k_6, 0], [0, -2*_Dummy_15150*k_5 - 2*_Dummy_15152*k_6 - _Dummy_15153*k_2, -_Dummy_15150*k_2, -2*_Dummy_15150*k_6, 0, 0, k_1, 0, 0, 0, 0, 0, -2*_Dummy_15161*k_5 - 2*_Dummy_15164*k_6 - _Dummy_15165*k_2 - 2*l, -_Dummy_15161*k_2, -2*_Dummy_15161*k_6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -_Dummy_15153*k_2, -_Dummy_15150*k_2 - _Dummy_15151*k_4, 0, -_Dummy_15153*k_4, 0, 0, 0, 0, 0, 0, 0, -_Dummy_15165*k_2, -_Dummy_15161*k_2 - _Dummy_15163*k_4 - 2*l, k_3, -_Dummy_15165*k_4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -2*_Dummy_15152*k_6 + _Dummy_15153*k_2, _Dummy_15150*k_2 + _Dummy_15151*k_4, -2*_Dummy_15150*k_6, _Dummy_15153*k_4, 0, 0, 0, 0, 0, 0, 0, -2*_Dummy_15164*k_6 + _Dummy_15165*k_2, _Dummy_15161*k_2 + _Dummy_15163*k_4, -2*_Dummy_15161*k_6 - k_3 - 2*l, _Dummy_15165*k_4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -_Dummy_15151*k_4, 0, -_Dummy_15153*k_4, 0, 0, 0, 0, 0, 0, 0, 0, -_Dummy_15163*k_4, k_3, -_Dummy_15165*k_4 - 2*l, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 2*_Dummy_15152*k_6, 0, 2*_Dummy_15150*k_6, 0, 0, 0, 0, 0, 0, 0, 0, 2*_Dummy_15164*k_6, 0, 2*_Dummy_15161*k_6, 0, -2*l, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -_Dummy_15143*k_2 - _Dummy_15148*k_2 - 2*_Dummy_15148*k_5 - 2*_Dummy_15157*k_6, -_Dummy_15140*k_2 - _Dummy_15146*k_4 - _Dummy_15148*k_2, -2*_Dummy_15148*k_6, -_Dummy_15148*k_4, 0, 0, -_Dummy_15165*k_2, -_Dummy_15161*k_2, 0, 0, 0, 0, k_1, 0, 0, 0, -_Dummy_15161*k_2 - 2*_Dummy_15161*k_5 - _Dummy_15163*k_4 - 2*_Dummy_15164*k_6 - _Dummy_15165*k_2 - 2*l, k_3, -_Dummy_15165*k_4, 0, -2*_Dummy_15161*k_6, 0, 0, 0, 0, 0], [0, -2*_Dummy_15144*k_6 - 2*_Dummy_15147*k_5 - 2*_Dummy_15147*k_6 + _Dummy_15148*k_2 - _Dummy_15157*k_2, _Dummy_15140*k_2 + _Dummy_15146*k_4 - _Dummy_15147*k_2, -2*_Dummy_15140*k_6 - 2*_Dummy_15147*k_6, _Dummy_15148*k_4, 0, 0, -2*_Dummy_15164*k_6 + _Dummy_15165*k_2, 0, -2*_Dummy_15161*k_6, 0, 0, 0, 0, k_1, 0, 0, _Dummy_15161*k_2 + _Dummy_15163*k_4, -2*_Dummy_15161*k_5 - 2*_Dummy_15161*k_6 - 2*_Dummy_15164*k_6 - _Dummy_15165*k_2 - k_3 - 2*l, _Dummy_15165*k_4, 0, -_Dummy_15161*k_2, 0, 0, 0, 0, 0], [0, -2*_Dummy_15146*k_5 - _Dummy_15156*k_2 - 2*_Dummy_15159*k_6, -_Dummy_15146*k_2 - _Dummy_15146*k_4, -2*_Dummy_15146*k_6, -_Dummy_15148*k_4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, k_1, 0, -_Dummy_15163*k_4, k_3, -2*_Dummy_15161*k_5 - 2*_Dummy_15164*k_6 - _Dummy_15165*k_2 - _Dummy_15165*k_4 - 2*l, 0, 0, -_Dummy_15161*k_2, 0, -2*_Dummy_15161*k_6, 0, 0], [0, -2*_Dummy_15145*k_5 + 2*_Dummy_15147*k_6 - _Dummy_15155*k_2 - 2*_Dummy_15158*k_6, -_Dummy_15145*k_2, 2*_Dummy_15140*k_6 - 2*_Dummy_15145*k_6, 0, 0, 0, 2*_Dummy_15164*k_6, 0, 0, 0, 0, 0, 0, 0, 0, k_1, 0, 2*_Dummy_15161*k_6, 0, -2*_Dummy_15161*k_5 - 2*_Dummy_15164*k_6 - _Dummy_15165*k_2 - 2*l, 0, 0, -_Dummy_15161*k_2, 0, -2*_Dummy_15161*k_6, 0], [0, _Dummy_15143*k_2 - _Dummy_15157*k_2 - 2*_Dummy_15157*k_6, -_Dummy_15147*k_2 + _Dummy_15148*k_2 + _Dummy_15156*k_4 - _Dummy_15159*k_4, -2*_Dummy_15148*k_6, _Dummy_15143*k_4 - _Dummy_15157*k_4, 0, 0, 0, _Dummy_15161*k_2 + _Dummy_15163*k_4, k_3, 0, 0, 0, 0, 0, 0, 0, -2*_Dummy_15164*k_6 + _Dummy_15165*k_2, -_Dummy_15165*k_2, 0, 0, -_Dummy_15161*k_2 - 2*_Dummy_15161*k_6 - _Dummy_15163*k_4 - k_3 - 2*l, _Dummy_15165*k_4, 0, -_Dummy_15165*k_4, 0, 0], [0, -_Dummy_15156*k_2, -_Dummy_15142*k_4 - _Dummy_15146*k_2 - _Dummy_15156*k_4, 0, -_Dummy_15143*k_4 - _Dummy_15156*k_4, 0, 0, 0, -_Dummy_15163*k_4, 0, -_Dummy_15165*k_4, 0, 0, 0, 0, 0, 0, 0, 0, -_Dummy_15165*k_2, 0, k_3, -_Dummy_15161*k_2 - _Dummy_15163*k_4 - _Dummy_15165*k_4 - 2*l, 0, k_3, 0, 0], [0, -_Dummy_15155*k_2 + 2*_Dummy_15157*k_6, -_Dummy_15145*k_2 - _Dummy_15154*k_4, 2*_Dummy_15148*k_6, -_Dummy_15155*k_4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*_Dummy_15164*k_6, 0, 0, -_Dummy_15165*k_2, 2*_Dummy_15161*k_6, 0, -_Dummy_15161*k_2 - _Dummy_15163*k_4 - 2*l, 0, k_3, -_Dummy_15165*k_4], [0, _Dummy_15156*k_2 - 2*_Dummy_15159*k_6, _Dummy_15142*k_4 + _Dummy_15146*k_2 - _Dummy_15159*k_4, -2*_Dummy_15146*k_6, _Dummy_15156*k_4 - _Dummy_15157*k_4, 0, 0, 0, 0, k_3, _Dummy_15165*k_4, 0, 0, 0, 0, 0, 0, 0, 0, -2*_Dummy_15164*k_6 + _Dummy_15165*k_2, 0, -_Dummy_15163*k_4, _Dummy_15161*k_2 + _Dummy_15163*k_4, 0, -2*_Dummy_15161*k_6 - _Dummy_15165*k_4 - k_3 - 2*l, 0, 0], [0, 2*_Dummy_15144*k_6 + _Dummy_15155*k_2 - 2*_Dummy_15158*k_6, _Dummy_15145*k_2 + _Dummy_15154*k_4, -2*_Dummy_15145*k_6 + 2*_Dummy_15147*k_6, _Dummy_15155*k_4, 0, 0, 0, 0, 2*_Dummy_15161*k_6, 0, 0, 0, 0, 0, 0, 0, 0, 2*_Dummy_15164*k_6, 0, -2*_Dummy_15164*k_6 + _Dummy_15165*k_2, 0, 0, _Dummy_15161*k_2 + _Dummy_15163*k_4, 0, -2*_Dummy_15161*k_6 - k_3 - 2*l, _Dummy_15165*k_4], [0, 2*_Dummy_15159*k_6, -_Dummy_15154*k_4, 2*_Dummy_15146*k_6, -_Dummy_15155*k_4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*_Dummy_15164*k_6, 0, 0, 0, -_Dummy_15163*k_4, 2*_Dummy_15161*k_6, k_3, -_Dummy_15165*k_4 - 2*l]])
    
    def chem_trend(self, state, t):
        
        out = np.zeros_like(state)
        
        state_mean = state[:6]
        state_var = state[6:12]
        state_aspect = state[27:]
        state_covar = state[12:27]

        ROC,RP,NO,NO2,O3,SNGN = state_mean
        V_ROC, V_RP,V_NO,V_NO2,V_O3,V_SNGN = state_var
        V_ROC_RP,V_ROC_NO,V_ROC_NO2,V_ROC_O3,V_ROC_SNGN,\
             V_RP_NO,V_RP_NO2,V_RP_O3,V_RP_SNGN,\
             V_NO_NO2,V_NO_O3,V_NO_SNGN,\
             V_NO2_O3,V_NO2_SNGN,\
             V_O3_SNGN = state_covar

        k3t = k3(t)
        k1t = k1(t)
        dtROC = eROC*self.cadastre -l*ROC
        dtRP = -k2*NO*RP-k2*V_RP_NO -k5*RP**2 - k5*V_RP -2*k6*NO2*RP-2*k6*V_RP_NO2-l*RP+ROC*k1t
        dtNO = eNO*self.cadastre-k2*NO*RP-k2*V_RP_NO-k4*NO*O3-k4*V_NO_O3-l*NO+NO2*k3t
        dtNO2 = eNO2*self.cadastre +k2*NO*RP+k2*V_RP_NO+k4*NO*O3+k4*V_NO_O3-2*k6*NO2*RP-2*k6*V_RP_NO2-l*NO2-NO2*k3t
        dtO3 = -k4*NO*O3-k4*V_NO_O3-l*O3+NO2*k3t
        dtSNGN = 2*k6*NO2*RP +2*k6*V_RP_NO2 -l*SNGN

        dtV_ROC = -2*l*V_ROC
        dtV_RP = -2*k2*NO*V_RP-2*k2*RP*V_RP_NO-4*k5*RP*V_RP -4*k6*NO2*V_RP-4*k6*RP*V_RP_NO2-2*l*V_RP+2*V_ROC_RP*k1t
        dtV_NO = -2*k2*NO*V_RP_NO-2*k2*RP*V_NO-2*k4*NO*V_NO_O3-2*k4*O3*V_NO-2*l*V_NO+2*V_NO_NO2*k3t 
        dtV_NO2 = 2*k2*NO*V_RP_NO2+2*k2*RP*V_NO_NO2+2*k4*NO*V_NO2_O3+2*k4*O3*V_NO_NO2 -4*k6*NO2*V_RP_NO2-4*k6*RP*V_NO2-2*l*V_NO2-2*k3t*V_NO2
        dtV_O3 = -2*k4*NO*V_O3-2*k4*O3*V_NO_O3-2*l*V_O3+2*V_NO2_O3*k3t
        dtV_SNGN = 4*k6*NO2*V_RP_SNGN + 4*k6*RP*V_NO2_SNGN -2*l*V_SNGN

        dtV_ROC_RP = -k2*NO*V_ROC_RP-k2*RP*V_ROC_NO-2*k5*RP*V_ROC_RP-2*k6*NO2*V_ROC_RP-2*k6*RP*V_ROC_NO2-2*l*V_ROC_RP+k1t*V_ROC
        dtV_ROC_NO = -k2*NO*V_ROC_RP-k2*RP*V_ROC_NO-k4*NO*V_ROC_O3-k4*O3*V_ROC_NO-2*l*V_ROC_NO+V_ROC_NO2*k3t
        dtV_ROC_NO2 = k2*NO*V_ROC_RP+k2*RP*V_ROC_NO+k4*NO*V_ROC_O3+k4*O3*V_ROC_NO-2*k6*NO2*V_ROC_RP -2*k6*RP*V_ROC_NO2-2*l*V_ROC_NO2-V_ROC_NO2*k3t
        dtV_ROC_O3 = -k4*NO*V_ROC_O3-k4*O3*V_ROC_NO-2*l*V_ROC_O3 +V_ROC_NO2*k3t
        dtV_ROC_SNGN = 2*k6*NO2*V_ROC_RP+2*k6*RP*V_ROC_NO2-2*l*V_ROC_SNGN

        dtV_RP_NO = -k2*NO*V_RP_NO-k2*NO*V_RP-k2*RP*V_RP_NO-k2*RP*V_NO-k4*NO*V_RP_O3-k4*O3*V_RP_NO-2*k5*RP*V_RP_NO-2*k6*NO2*V_RP_NO-2*k6*RP*V_NO_NO2-2*l*V_RP_NO+V_ROC_NO*k1t + V_RP_NO2*k3t
        dtV_RP_NO2 = -k2*NO*V_RP_NO2+k2*NO*V_RP-k2*RP*V_NO_NO2+k2*RP*V_RP_NO+k4*V_RP_O3+k4*O3*V_RP_NO-2*k5*RP*V_RP_NO2-2*k6*NO2*V_RP_NO2-2*k6*NO2*V_RP-2*k6*RP*V_RP_NO2-2*k6*RP*V_NO2-2*l*V_RP_NO2+V_ROC_NO2*k1t-V_RP_NO2*k3t
        dtV_RP_O3 = -k2*NO*V_RP_O3-k2*RP*V_NO_O3-k4*NO*V_RP_O3-k4*O3*V_RP_NO-2*k5*RP*V_RP_O3-2*k6*NO2*V_RP_O3-2*k6*RP*V_NO2_O3-2*l*V_RP_O3 + V_ROC_O3*k1t+V_RP_NO2*k3t
        dtV_RP_SNGN = -k2*NO*V_RP_SNGN-k2*RP*V_NO_SNGN-2*k5*RP*V_RP_SNGN-2*k6*NO2*V_RP_SNGN+2*k6*NO2*V_RP-2*k6*RP*V_NO2_SNGN+2*k6*RP*V_RP_NO2-2*l*V_RP_SNGN+V_ROC_SNGN*k1t

        dtV_NO_NO2 = k2*NO*V_RP_NO-k2*NO*V_RP_NO2-k2*RP*V_NO_NO2+k2*RP*V_NO-k4*NO*V_NO2_O3+k4*NO*V_NO_O3-k4*O3*V_NO_NO2+k4*O3*V_NO-2*k6*NO2*V_RP_NO-2*k6*RP*V_NO_NO2-2*l*V_NO_NO2-V_NO_NO2*k3t +k3t*V_NO2
        dtV_NO_O3 = -k2*NO*V_RP_O3-k2*RP*V_NO_O3-k4*NO*V_NO_O3-k4*NO*V_O3-k4*O3*V_NO_O3-k4*O3*V_NO-2*l*V_NO_O3+V_NO2_O3*k3t +V_NO_NO2*k3t
        dtV_NO_SNGN = -k2*NO*V_RP_SNGN-k2*RP*V_NO_SNGN-k4*NO*V_O3_SNGN-k4*O3*V_NO_SNGN +2*k6*NO2*V_RP_NO+2*k6*RP*V_NO_NO2-2*l*V_NO_SNGN+V_NO2_SNGN*k3t

        dtV_NO2_O3 = k2*NO*V_RP_O3+k2*RP*V_NO_O3-k4*NO*V_NO2_O3+k4*NO*V_O3-k4*O3*V_NO_NO2+k4*O3*V_NO_O3-2*k6*NO2*V_RP_O3-2*k6*RP*V_NO2_O3-2*l*V_NO2_O3-V_NO2_O3*k3t+k3t*V_NO2
        dtV_NO2_SNGN = k2*NO*V_RP_SNGN+k2*RP*V_NO_SNGN+k4*NO*V_O3_SNGN+k4*O3*V_NO_SNGN+2*k6*NO2*V_RP_NO2-2*k6*NO2*V_RP_SNGN-2*k6*RP*V_NO2_SNGN+2*k6*RP*V_NO2-2*l*V_NO2_SNGN-V_NO2_SNGN*k3t

        dtV_O3_SNGN = -k4*NO*V_O3_SNGN-k4*O3*V_NO_SNGN+2*k6*NO2*V_RP_O3+2*k6*RP*V_NO2_O3-2*l*V_O3_SNGN+V_NO2_SNGN*k3t

        chem_trend_mean = np.array([dtROC,dtRP,dtNO,dtNO2,dtO3,dtSNGN])
        chem_trend_var = np.array([ dtV_ROC,dtV_RP,dtV_NO,dtV_NO2,dtV_O3,dtV_SNGN])
        chem_trend_covar = np.array([dtV_ROC_RP,dtV_ROC_NO,dtV_ROC_NO2,dtV_ROC_O3,dtV_ROC_SNGN,\
           dtV_RP_NO,dtV_RP_NO2,dtV_RP_O3,dtV_RP_SNGN,\
           dtV_NO_NO2,dtV_NO_O3,dtV_NO_SNGN,\
           dtV_NO2_O3,dtV_NO2_SNGN,\
           dtV_O3_SNGN])
        
        out[:6] += chem_trend_mean
        out[6:12] += chem_trend_var
        out[12:27] += chem_trend_covar
        
        return out
        


