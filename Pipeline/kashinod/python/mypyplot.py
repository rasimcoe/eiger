import numpy as np

def imshow_ptiles(ax, img, ptiles=[0,100], origin='lower', **kwargs):
    vmin,vmax = np.nanpercentile(img, ptiles)
    ax.imshow(img, vmin=vmin, vmax=vmax, **kwargs, origin=origin)
    
