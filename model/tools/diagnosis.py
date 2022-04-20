
import numpy as np

def diag_with_closures(traj_enkf):
    diag = {}
    
    for t, _ensemble in traj_enkf.items() : 
        ensemble = np.swapaxes(np.array(_ensemble),0,1)
        a = ensemble[0].mean(axis=0)
        b = ensemble[1].mean(axis=0)
        dx = 1/ ensemble.shape[-1]
        

        
        std_a = ensemble[0].std(axis=0)
        std_b = ensemble[1].std(axis=0)
        
        errors_A = ensemble[0] - a 
        errors_B = ensemble[1] - b 
        
        vab = (errors_A*errors_B).mean(axis=0)
        norm_errors_A = errors_A/std_a
        norm_errors_B = errors_B/std_b

        dx_norm_errors_A =   (np.roll(norm_errors_A,-1,axis=1) - np.roll(norm_errors_A,1,axis=1))/2/dx
        dx_norm_errors_B =   (np.roll(norm_errors_B,-1,axis=1) - np.roll(norm_errors_B,1,axis=1))/2/dx        
        
        sa = 1/(dx_norm_errors_A**2).mean(axis=0)
        sb = 1/(dx_norm_errors_B**2).mean(axis=0)
        
        eb_dx_ea = (norm_errors_B*dx_norm_errors_A).mean(axis=0)
        ea_dx_eb = (dx_norm_errors_B*norm_errors_A).mean(axis=0)
        dx_ea_dx_eb = (dx_norm_errors_B*dx_norm_errors_A).mean(axis=0)
        
        diag[t] = {}
        diag[t]['A'] = a; diag[t]['B'] = b;
        diag[t]['Va'] = std_a**2; diag[t]['Vb'] = std_b**2;
        diag[t]['Vab'] = vab
        diag[t]['sA'] = sa; diag[t]['sB'] = sb;
        diag[t]['eb_dx_ea'] = eb_dx_ea
        diag[t]['ea_dx_eb'] = ea_dx_eb
        diag[t]['dx_ea_dx_eb'] = dx_ea_dx_eb
    
    return diag
