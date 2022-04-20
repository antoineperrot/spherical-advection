import matplotlib.pyplot as plt
import numpy as np

def univariate_plot(times,x, enkf_traj, pkf_traj, normalizations, ensemble_size):
    fig, ax = plt.subplots(1,3,figsize=(16,5))
    fig.subplots_adjust(hspace=.25,wspace=.25)
    panel_index = [f'({c}) ' for c in 'abc']
    keys = ['Mean concentration','Std','Length-scale']
    fields = ['Mean concentration','Standard-deviation','Length-scale']
    normalization_labels=[r'$A/A^0$', r'$\sigma_x/\sigma^0$',r'$L_x/l_h^0$']

    alphas = np.linspace(0,1,len(times))

    for k, field in enumerate(fields):
        
        ax[k].set_title(panel_index[k] + field,fontsize=15)
        for i, time in enumerate(times):
            alpha = alphas[i]**2
            linewidth = 2 if time in [times[0],times[-1]] else 1
            label = f'{ensemble_size}-Ensemble' if time == times[-1] else None
            c = 'b' 
            ax[k].set_ylabel(normalization_labels[k],fontsize=15)
            ax[k].plot(x, enkf_traj[time][keys[k]]['X']/normalizations[k],c=c,
                     linewidth=linewidth,label=label,linestyle='--',alpha=alpha)

            c = 'r' #if i in [0,n_stops-1] else 'lightsalmon'
            label = 'PKF-NN' if time == times[-1] else None
            ax[k].plot(x, pkf_traj[time][keys[k]]['X']/normalizations[k],c=c,
                     linewidth=linewidth,linestyle='-',label=label,zorder=3,alpha=alpha)
            ax[k].set_xlabel(r'$\mathbf{x}/D$',fontsize=15)
            if k==0 and time == times[-1] : ax[k].legend(); ax[k].text(.73,1.53,r'$t=3 \tau_{adv}$',fontsize=14)
                
    return fig
