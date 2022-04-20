import numpy as np

class MultivariateEnKF1D(object):   
    def __init__(self, species_name, domain):
        self.species_name = species_name
        self.n_species = len(species_name)
        self.domain = domain
        
        self.mapping = {}
        for i, specie in enumerate(self.species_name) : self.mapping[specie] = i
        
    def _correlation_matrix(self, ensemble):
        tmp = np.hstack([ensemble[i] for i in range(ensemble.shape[0])]).T
        return np.corrcoef(tmp)

    def _covariance_matrix(self, ensemble):
	    tmp = np.hstack([ensemble[i] for i in range(ensemble.shape[0])]).T
	    return np.cov(tmp)
        
    def assimilate(self, enkf_state, obs_network, obs_species, obs_values, obs_variances):
        
        state_f = np.hstack([enkf_state[i].mean(axis=0) for i in range(self.n_species) ])
        error_f = self._covariance_matrix(enkf_state)
        

        #observation operator : 
        H = np.zeros((len(obs_network),self.domain.Nx*self.n_species))
        
        
        for i, (obs_loc, obs_specie) in enumerate(zip(obs_network, obs_species)):
        	index_specie = self.mapping[obs_specie]
        	
        	H[i, obs_loc + index_specie * self.domain.Nx] = 1
        	

        # observation error covariance matrix
        R = np.diag(obs_variances)

        # kalman gain
        K = error_f @ H.T @ np.linalg.inv(H @ error_f @ H.T + R)

        #analysis step :
        error_a = (np.eye(self.n_species*self.domain.Nx) - K @ H) @ error_f
        state_a = state_f + K @ (obs_values - H @ state_f)
        
        
        #new analysis ensemble :
        ensemble_size = enkf_state.shape[1]
        samples_a = np.random.multivariate_normal(state_a, error_a, ensemble_size)
        ensemble_a = []
        for i in range(self.n_species):
            ensemble_a.append(samples_a[:,i*self.domain.Nx:(i+1)*self.domain.Nx])
	      
        ensemble_a = np.array(ensemble_a)

        return ensemble_a
        
    def _formated_correlation_matrix(self, ensemble):
        cor_mat = self._correlation_matrix(ensemble)
        out = {}
        Nx = self.domain.Nx
        for i, specie in enumerate(self.species_name):
            out[specie] = {}
            for j, specie2 in enumerate(self.species_name):
                out[specie][specie2] = cor_mat[i*Nx:(i+1)*Nx, j*Nx:(j+1)*Nx]
        return out
    
    def _formated_covariance_matrix(self, ensemble):
        cov_mat = self._covariance_matrix(ensemble)
        out = {}
        Nx = self.domain.Nx
        for i, specie in enumerate(self.species_name):
            out[specie] = {}
            for j, specie2 in enumerate(self.species_name):
                out[specie][specie2] = cov_mat[i*Nx:(i+1)*Nx, j*Nx:(j+1)*Nx]
        return out
    
    
    def diagnosis(self, traj, cor_func_locations= [], cor_func_species=[]):
        diagnosis = {}
        times = traj.keys()
        for time in times:
            diagnosis[time]= {'Mean concentration':{}, 'Std':{}, 'Length-scale':{}, 'Correlation': {} }
            for i,species in enumerate(self.species_name):
                diagnosis[time]['Mean concentration'][species] = traj[time][i].mean(axis=0)
                diagnosis[time]['Std'][species] = traj[time][i].std(axis=0)
                errors = traj[time][i] - traj[time][i].mean(axis=0)
                normalized_errors = errors/errors.std(axis=0)

                dx_normalized_errors = self.domain.derive(normalized_errors,axis=1)
                diagnosis[time]['Length-scale'][species] = (1/np.mean(dx_normalized_errors**2,axis=0))**.5
                diagnosis[time]['Correlation'][species] = {}                
                
            for i,species in enumerate(self.species_name):
                for j in range(i+1, self.n_species):
                    errors = traj[time][i] - traj[time][i].mean(axis=0)
                    errors_j = traj[time][j] - traj[time][j].mean(axis=0)
                    species_j = self.species_name[j]
                    diagnosis[time]['Correlation'][species][species_j] = np.mean(errors*errors_j,axis=0)/errors.std(axis=0)/errors_j.std(axis=0)
                    diagnosis[time]['Correlation'][species_j][species] = np.mean(errors*errors_j,axis=0)/errors.std(axis=0)/errors_j.std(axis=0)
                    
            diagnosis[time]['Correlation functions'] = {}
            diagnosis[time]['Covariance functions'] = {}
            cor_mat = self._formated_correlation_matrix(traj[time])
            cov_mat = self._formated_covariance_matrix(traj[time])
            for loc, specie in zip(cor_func_locations, cor_func_species):
                diagnosis[time]['Correlation functions'][loc] = {}
                diagnosis[time]['Covariance functions'][loc] = {}
                for specie2 in self.species_name:
                    diagnosis[time]['Correlation functions'][loc][specie+'/'+specie2] = cor_mat[specie][specie2][loc]
                    diagnosis[time]['Covariance functions'][loc][specie+'/'+specie2] = cov_mat[specie][specie2][loc]
                
                
        return diagnosis

