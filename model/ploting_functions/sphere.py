import matplotlib.pyplot as plt
from matplotlib import colors, cm
def plot_sphere(X,Y,Z, values, cmap=cm.PuRd, ds=1,cmin=None, cmax=None,elev=0,azim=0,figsize=(16,16),title=None):
    fig = plt.figure(figsize=figsize)
    fig.suptitle(title,fontsize=17)
    ax = fig.add_subplot(1,1,1,projection='3d')
    ax.view_init(elev=elev, azim=azim)

    strength = values[::ds, ::ds]
    norm=colors.Normalize(vmin = np.min(strength) if cmin is None else cmin ,
                          vmax = np.max(strength) if cmax is None else cmax,
                          clip = False)

    surface = ax.plot_surface(X[::ds, ::ds],Y[::ds, ::ds],Z[::ds, ::ds], rstride=1, cstride=1,cmap = cm.PuRd,
                           linewidth=0, antialiased=False,
                           facecolors=cmap(norm(strength)),zorder=1)

    ax.set_xlabel("$\mathbf{x}$",fontsize=15)
    ax.set_ylabel("$\mathbf{y}$",fontsize=15)
    ax.set_zlabel("$\mathbf{z}$",fontsize=15);
    return fig
