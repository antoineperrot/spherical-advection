import numpy as np
import matplotlib.pyplot as plt

def plot_univariate_statistics(xx, times_stop, traj, species_name, normalization,
                               colors, format_ = (6,3), fig =  None, ax=None, suptitle=None, **kwargs):
    n_species = len(species_name)
    import string
    panels_index = [f'({c}) ' for c in string.ascii_lowercase]
    
    if fig is None : 
        figsize=(18,9*3)
        fig, ax = plt.subplots(*format_,figsize=figsize)
    if suptitle is not None : fig.suptitle(suptitle, fontsize=16)
    n_stops = len(times_stop)
    for k, type_field in enumerate(['Mean concentration','Std','Length-scale']):
        for i, field in enumerate(species_name):
            index_plot = i+k*n_species
            ax.flat[index_plot].set_title(panels_index[index_plot] + type_field + ' '+field,fontsize=14)  
            if k%2 == 1 : ax.flat[i+k*n_species].set_facecolor(3*(0.96,))
            for j,time in enumerate(times_stop):
                alpha = ((j+1)/(n_stops+1))**2
                ax.flat[index_plot].plot(xx, traj[time][type_field +' '+ field]/\
                                normalization[type_field][field],c=colors[field],alpha=alpha,
                                        **kwargs)
                
    fig.tight_layout()
    return fig, ax
    
    
def plot_multivariate_statistics(xx, times, traj, species_name, normalization,
                                 colors, figsize=26, fig=None, ax=None, **kwargs):
    n_species = len(species_name) ; n_stops = len(times);
    if fig is None:
        fig, ax = plt.subplots(n_species, n_species, figsize =2*(figsize,))
        
    import copy
    kwargs2 = copy.deepcopy(kwargs)
    del kwargs2['c']
    
    for i, field1 in enumerate(species_name):
        for j,field2 in enumerate(species_name):
            if j>i:
                ax[i,j].set_title(field1+"/"+field2,fontsize=14)    
                ax[j,i].set_title(field2+"/"+field1,fontsize=14)

                for k,time  in enumerate(times):
                    alpha = ((k+1)/(n_stops+1))**2
                    normalization_ = (traj[time]['Std '+field1]*traj[time]['Std '+field2])

                    ax[i,j].plot(domain1d.x,traj[time]['Covariance ' + field1+'/'+field2]/normalization_,
                                    **kwargs)
                    ax[j,i].plot(domain1d.x,traj[time]['Covariance ' + field1+'/'+field2]/normalization_,
                                    **kwargs)

                    ax[i,j].set_ylim(-1.1,1.1);ax[j,i].set_ylim(-1.1,1.1)
            if j==i :
                ax[i,j].set_title('Std '+field1, fontsize=14)    

                for k,time  in enumerate(times):
                    alpha = ((k+1)/(n_stops+1))**2
                    normalization_ = normalization['Std'][field1]
                    ax[i,j].plot(domain1d.x,diag[time]['Std ' + field1]/normalization_,
                                    c=colors[field],alpha=alpha,**kwargs2)
    fig.tight_layout();
    
    return fig, ax
