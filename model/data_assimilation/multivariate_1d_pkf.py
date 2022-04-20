import numpy as np

class MultivariatePKF1D(object):   
    def __init__(self, species_name, domain):
        self.species_name = species_name
        self.n_species = len(species_name)
        self.domain = domain
        
        self.mapping = {}
        self.mapping['Mean concentration'] = {}
        self.mapping['Std'] = {}
        self.mapping['Length-scale'] = {}
        self.mapping['Covariance'] = {}
        total_n_fields = int(self.n_species * 3 + (self.n_species-1)*self.n_species / 2)
        cov_index = list(range(int((self.n_species-1)*self.n_species / 2)))
        for i, specie in enumerate(species_name):
            self.mapping['Mean concentration'][specie] = i
            self.mapping['Std'][specie] = self.n_species + i
            self.mapping['Length-scale'][specie] = int(2*self.n_species + (self.n_species-1)*self.n_species / 2  + i)
            self.mapping['Covariance'][specie]= {}
        for i, specie in enumerate(species_name):
            for j, specie2 in enumerate(self.species_name):
                if j>i:
                    self.mapping['Covariance'][specie][specie2] = 2*self.n_species + cov_index[0]
                    self.mapping['Covariance'][specie2][specie] = 2*self.n_species + cov_index[0]
                    cov_index.pop(0)

    def _autocorrelation_model(self, s, loc):
        dist = np.abs(self.domain.x - self.domain.x[loc])
        dist = np.minimum(dist, 0.5- (dist -0.5))
        output = (s[loc]*s)**.25/(0.5*(s[loc]+s))**.5 * \
            np.exp(- dist**2/ (s[loc]+s))
        return output
    
    def _autocovariance_model(self, V, s, loc):
        dist = np.abs(self.domain.x - self.domain.x[loc])
        dist = np.minimum(dist, 0.5- (dist -0.5))
        output = (V[loc]*V)**.5 *\
            (s[loc]*s)**.25/(0.5*(s[loc]+s))**.5 * \
            np.exp(- dist**2/ (s[loc]+s))
        return output
    
    def _crosscorrelation_model(self, Vab,Va,Vb,sA,sB, loc):
        dist = np.abs(self.domain.x - self.domain.x[loc])
        dist = np.minimum(dist, 0.5- (dist -0.5))
        output = 0.5 * (Vab[loc]/(Va[loc]*Vb[loc])**.5 + Vab/(Va*Vb)**.5) * \
                np.exp(-dist**2/( 0.5*(sA[loc]+sA +sB[loc]+sB) ))
        return output
        

    def _crosscovariance_model(self, Vab,Va,Vb,sA,sB, loc):
        dist = np.abs(self.domain.x - self.domain.x[loc])
        dist = np.minimum(dist, 0.5- (dist -0.5))
        output = 0.5 * (Vab[loc]/(Va[loc]*Vb[loc])**.5 + Vab/(Va*Vb)**.5) * \
                np.exp(-dist**2/( 0.5*(sA[loc]+sA+sB[loc]+sB) )) * \
                np.sqrt(Va[loc]*Vb)
        return output
    
    def _get_covariance_func(self, pkf_state, loc, specie):
        funcs = {}
        for i, specie2 in enumerate(self.species_name):
            if specie2 == specie :
                func = self._autocovariance_model(pkf_state[self.mapping["Std"][specie]],
                                                  pkf_state[self.mapping["Length-scale"][specie]], loc)
            else :
                func = self._crosscovariance_model(
                pkf_state[self.mapping['Covariance'][specie][specie2]],
                pkf_state[self.mapping['Std'][specie]],
                pkf_state[self.mapping['Std'][specie2]],
                pkf_state[self.mapping['Length-scale'][specie]],
                pkf_state[self.mapping['Length-scale'][specie2]],
                loc)
            funcs[specie + '/'+specie2] = func
        return funcs
    
    def _get_correlation_func(self, pkf_state, loc, specie):
        funcs = {}
        for i, specie2 in enumerate(self.species_name):
            if specie2 == specie :
                func = self._autocorrelation_model(pkf_state[self.mapping["Length-scale"][specie]], loc)
            else :
                func = self._crosscorrelation_model(
                pkf_state[self.mapping['Covariance'][specie][specie2]],
                pkf_state[self.mapping['Std'][specie]],
                pkf_state[self.mapping['Std'][specie2]],
                pkf_state[self.mapping['Length-scale'][specie]],
                pkf_state[self.mapping['Length-scale'][specie2]],
                loc)
            funcs[specie + '/'+specie2] = func
        return funcs
            
                        
    def render_traj(self, traj, cor_func_locations = [], cor_func_species=[]):
        rendered_traj = {}
        times = traj.keys()
        for time in times :
            rendered_traj[time]= {'Mean concentration':{}, 'Std':{}, 'Length-scale':{}, 'Correlation': {} }
            for i, specie in enumerate(self.species_name) :       
                rendered_traj[time]['Mean concentration'][specie] = traj[time][self.mapping['Mean concentration'][specie]]
                rendered_traj[time]['Std'][specie] = traj[time][self.mapping['Std'][specie]]**.5
                rendered_traj[time]['Length-scale'][specie] = traj[time][self.mapping['Length-scale'][specie]]**.5   
                rendered_traj[time]['Correlation'][specie] = {}
                
            for i, specie in enumerate(self.species_name) :
                for j, specie2 in enumerate(self.species_name):
                    if j > i :
                        rendered_traj[time]['Correlation'][specie][specie2] = traj[time][self.mapping['Covariance'][specie][specie2]]/traj[time][self.mapping['Std'][specie]]**.5/traj[time][self.mapping['Std'][specie2]]**.5
                        rendered_traj[time]['Correlation'][specie2][specie] = traj[time][self.mapping['Covariance'][specie][specie2]]/traj[time][self.mapping['Std'][specie]]**.5/traj[time][self.mapping['Std'][specie2]]**.5
                            
                rendered_traj[time]['Correlation functions'] = {}
                rendered_traj[time]['Covariance functions'] = {}
                for loc, specie in zip(cor_func_locations, cor_func_species):
                    rendered_traj[time]['Correlation functions'][loc] = self._get_correlation_func(traj[time], loc, specie)
                    rendered_traj[time]['Covariance functions'][loc] = self._get_covariance_func(traj[time], loc, specie)
        
        return rendered_traj
    
    
    def multivariate_assimilation(self, pkf_state, obs_values, obs_locations,
                                  obs_species, obs_variances, method='O2'):
        Nx = self.domain.Nx
        
        for obs_value, obs_loc, obs_specie, obs_var in zip(obs_values, obs_locations,
                                  obs_species, obs_variances):
            
            pkf_state_a = np.zeros_like(pkf_state)
            
            x_f = np.hstack([pkf_state[i] for i in range(self.n_species)])
            v_f = np.hstack([pkf_state[i+self.n_species] for i in range(self.n_species)])
            g_f = 1/np.hstack([pkf_state[i+2*self.n_species] for i in range(self.n_species)])
            d_v_f = np.hstack([self.domain.derive(pkf_state[i+self.n_species]) for i in range(self.n_species)])
                              
            cor_funcs = []
            d_cor_funcs = []

            x_f_loc = pkf_state[self.mapping['Mean concentration'][obs_specie]][obs_loc]
            v_f_loc = pkf_state[self.mapping['Std'][obs_specie]][obs_loc]
            for k, specie in enumerate(self.species_name):
                if specie == obs_specie:
                    func = self._autocorrelation_model(pkf_state[self.mapping['Length-scale'][specie]],
                                                       obs_loc)
                    auto_var_func = pkf_state[k+self.n_species]
                    auto_cov_func = self._autocovariance_model(pkf_state[self.mapping['Std'][specie]],
                                                                pkf_state[self.mapping['Length-scale'][specie]],
                                           obs_loc)
                else :
                    func = self._crosscorrelation_model(
                        pkf_state[self.mapping['Covariance'][obs_specie][specie]],
                        pkf_state[self.mapping['Std'][obs_specie]],
                        pkf_state[self.mapping['Std'][specie]],
                        pkf_state[self.mapping['Length-scale'][obs_specie]],
                        pkf_state[self.mapping['Length-scale'][specie]],
                        obs_loc)
                    
                cor_funcs.append(func)
                d_cor_funcs.append(self.domain.derive(func))
            cor_funcs = np.hstack(cor_funcs);
            d_cor_funcs = np.hstack(d_cor_funcs)
                              
            # state update :
            x_a = x_f + v_f**.5 * cor_funcs * v_f_loc**.5/ (v_f_loc + obs_var) *(obs_value - x_f_loc)
            
            #variance update :
            v_a = v_f * (1 - cor_funcs**2 * v_f_loc / (v_f_loc + obs_var))
            
            # metric tensor update :
            if method == 'O1':
                g_a = v_f/v_a * g_f 
            
            if method == 'O2':
                tmp  = np.hstack([self.domain.derive((v_f**.5 * cor_funcs)[i*Nx:(i+1)*Nx]) for i in range(self.n_species)])
                tmp2 = np.hstack([self.domain.derive((v_a)[i*Nx:(i+1)*Nx]) for i in range(self.n_species)])
                g_a = v_f/v_a * g_f + 1/(4*v_f*v_a) * (d_v_f)**2 \
                    - 1/v_a *tmp**2* v_f_loc/(v_f_loc + obs_var) \
                    - 1/(4*v_a**2) * tmp2**2
            s_a = 1/g_a
            
            # covariance update :
            for k1, sp1 in enumerate(self.species_name):
                for k2, sp2 in enumerate(self.species_name) :
                    if k2 > k1 :
                        v_ab_f = pkf_state[self.mapping['Covariance'][sp1][sp2]]
                        if sp1 == obs_specie :
                            cov_func_a_l = auto_cov_func
                        else :                    	
                            cov_func_a_l = self._crosscovariance_model(
                                    pkf_state[self.mapping['Covariance'][obs_specie][sp1]],
                                    pkf_state[self.mapping['Std'][obs_specie]],
                                    pkf_state[self.mapping['Std'][sp1]],
                                    pkf_state[self.mapping['Length-scale'][obs_specie]],
                                    pkf_state[self.mapping['Length-scale'][sp1]],
                                    obs_loc)
                        if sp2 == obs_specie :
                            cov_func_b_l = auto_cov_func
                        else :
                            cov_func_b_l = self._crosscovariance_model(
                                    pkf_state[self.mapping['Covariance'][obs_specie][sp2]],
                                    pkf_state[self.mapping['Std'][obs_specie]],
                                    pkf_state[self.mapping['Std'][sp2]],
                                    pkf_state[self.mapping['Length-scale'][obs_specie]],
                                    pkf_state[self.mapping['Length-scale'][sp2]],
                                    obs_loc)
                                    
                        v_ab_a = v_ab_f - 1/(v_f_loc + obs_var) * cov_func_a_l * cov_func_b_l
                        pkf_state_a[self.mapping['Covariance'][sp1][sp2]] = v_ab_a

            for k, specie in enumerate(self.species_name):
                pkf_state_a[self.mapping['Mean concentration'][specie]] = x_a[k*Nx:(k+1)*Nx]
                pkf_state_a[self.mapping['Std'][specie]] = v_a[k*Nx:(k+1)*Nx]
                pkf_state_a[self.mapping['Length-scale'][specie]] = s_a[k*Nx:(k+1)*Nx]
                
            
            pkf_state = pkf_state_a
            
        return pkf_state
