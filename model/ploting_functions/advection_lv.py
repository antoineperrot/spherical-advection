import numpy as np
import matplotlib.pyplot as plt

def Comparison_Anisotropies_LV_HO(diag_enkf, initial_settings,k,title='', c1 = 'fuchsia', c2 = "lawngreen"):
    times = np.array(list(diag_enkf.keys()))
    
    stdA = np.array([diag_enkf[t]['Va'] for t in times]).mean(axis=1)**.5
    stdB = np.array([diag_enkf[t]['Vb'] for t in times]).mean(axis=1)**.5
    la = np.array([diag_enkf[t]['sA'] for t in times]).mean(axis=1)**.5
    lb = np.array([diag_enkf[t]['sB'] for t in times]).mean(axis=1)**.5
    Vab = np.array([diag_enkf[t]['Vab'] for t in times]).mean(axis=1)
    dx_ea_dx_eb = np.array([diag_enkf[t]['dx_ea_dx_eb'] for t in times]).mean(axis=1)
    cAB = Vab/stdA/stdB
    
    fig, ax = plt.subplots(4,2,figsize=(16,8),sharex='all')

    fig.suptitle(title,fontsize=15)
    
    ax[0,0].plot(times,2*stdA/(stdA[0]+stdB[0]),label='$2 \sigma_A/(\sigma_A^0+ \sigma_B^0)$',c=c1)
    ax[0,0].plot(times,2*stdB/(stdA[0]+stdB[0]),label='$2 \sigma_B/(\sigma_A^0+ \sigma_B^0)$',c=c2)
    
    ax[1,0].plot(times,cAB,label=r'$corr(A, B)$',c='gray')
    ax[1,0].set_ylim(-1,1)

    ax[2,0].plot(times,2*la/(la[0]+lb[0]),label='$2 L_A/(L_A^0+ L_B^0)$',c=c1)
    ax[2,0].plot(times,2*lb/(la[0]+lb[0]),label='$2 L_B/(L_A^0+ L_B^0)$',c=c2)
    
    ax[3,0].plot(times,dx_ea_dx_eb,label=r'$\overline{\partial_x \tilde{\varepsilon}_A\partial_x \tilde{\varepsilon}_B}$',c='gray')
    ax[3,0].set_xlabel("$t$",fontsize=14)
    
    time_window = np.array(times)
    
    LA0 = initial_settings['Length-scale A']
    LB0 = initial_settings['Length-scale B']
    
    VA0 = initial_settings['std A']**2
    VB0 = initial_settings['std B']**2
    
    gA0, gB0 = 1/LA0**2, 1/LB0**2
    
    # Calcul des séries temporelles des statistiques
    ct, st = np.cos(k*time_window), np.sin(k*time_window)
    
    VA = ct**2 *VA0 + st**2*VB0
    VB = ct**2 *VB0 + st**2*VA0
    VAB = ct*st*(VA0-VB0) 
    
    stdA, stdB = np.sqrt(VA), np.sqrt(VB) 
    
    gA = (ct**2 *VA0 *gA0 + st**2*VB0 *gB0)/VA
    gB = (st**2 *VA0 *gA0 + ct**2*VB0 *gB0)/VB
    
    LA = 1/np.sqrt(gA)
    LB = 1/np.sqrt(gB)

    EAB = ct*st/(stdA*stdB)*(VA0*gA0-VB0*gB0)
    
    
    ax[0,1].plot(times,2*stdA/(stdA[0]+stdB[0]),label='$2 \sigma_A/(\sigma_A^0+ \sigma_B^0)$',c=c1,linestyle='--')
    ax[0,1].plot(times,2*stdB/(stdA[0]+stdB[0]),label='$2 \sigma_B/(\sigma_A^0+ \sigma_B^0)$',c=c2,linestyle='--')
    
    
    ax[1,1].plot(times,VAB/stdA/stdB,label=r'$corr(A, B)$',c='gray',linestyle='--')
    ax[1,1].set_ylim(-1,1)
    
    ax[2,1].plot(times,2*LA/(LA[0]+LB[0]),label='$2 L_A/(L_A^0+ L_B^0)$',c=c1,linestyle='--')
    ax[2,1].plot(times,2*LB/(LA[0]+LB[0]),label='$2 L_B/(L_A^0+ L_B^0)$',c=c2,linestyle='--')

    ax[2,0].set_ylim(0.8,1.2)
    ax[2,1].set_ylim(0.8,1.2)
    
    ax[3,1].plot(times,EAB,label=r'$\overline{\partial_x \tilde{\varepsilon}_A\partial_x \tilde{\varepsilon}_B}$',c='gray',linestyle='--')
    ax[3,1].set_xlabel("$t$",fontsize=14)
    ax[3,1].axhline(0,c='grey',alpha=0.3)
    ax[3,0].axhline(0,c='grey',alpha=0.3)
    ax[1,0].axhline(0,c='grey',alpha=0.3)
    ax[1,1].axhline(0,c='grey',alpha=0.3)
    
    fields = ['normalized std','correlation', 'normalized length-scales',r'$\overline{\partial_x \tilde{\varepsilon}_A\partial_x \tilde{\varepsilon}_B}$']
    system = ['in the LV dyn.','in the HO dyn.']
    
    for i, (panel, index) in enumerate(zip(ax.flat, 'abcdefghij')):
        panel.legend();
        panel.set_title('('+index+') '+ fields[i//2] +' ' + system[i%2],fontsize=14);
    
    fig.tight_layout();
    return fig;