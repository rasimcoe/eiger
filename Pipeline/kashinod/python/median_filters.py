import numpy as np


def median_filter1d(array, kernel, hole, nan_preserve=False):

    array = np.array(array)
    array_ext = np.append(np.append(np.nan, array),np.nan)

    n_x = array_ext.size
    
    ## NaN value will be ignored
    kx = kernel
    idx_x = np.fromfunction(lambda  i, j: i+j, (n_x, kx), dtype=np.int64) - kx // 2

    # Exclude the center gap
    if hole > 0:
        kx_gap = hole
        #print('Gap size: ', kx_gap, flush=True)
        cols_remain=np.append(np.arange((kx - kx_gap)/2., dtype=np.int64),
                              np.flip(kx-1-np.arange((kx - kx_gap)/2., dtype=np.int64)))
        idx_x = idx_x[:, cols_remain]

    idx_x[idx_x < 0]=0
    idx_x[idx_x > n_x-1]=n_x-1

    med_arr = np.nanmedian(array_ext[idx_x], axis=1)
    if nan_preserve==True:
        med_arr[np.isnan(array_ext)]=np.nan
        
    return med_arr[1:-1]

def median_filter2d(array, kx, ky, kx_gap, nan_preserve=False):
    
    array = np.array(array)
    n_y0, n_x0 = array.shape
    n_y = n_y0 + 2
    n_x = n_x0 + 2

    
    array_ext = np.zeros((n_y,n_x))+np.nan
    array_ext[1:-1, 1:-1] = array

    ## NaN value will be ignored
    idx_x = np.fromfunction(lambda  i, j: i+j, (n_x, kx), dtype=np.int64) - kx // 2
    idx_y = np.fromfunction(lambda  i, j: i+j, (n_y, ky), dtype=np.int64) - ky // 2

    # Exclude the center gap
    if kx_gap > 0:
        #print('Gap size: ', kx_gap, flush=True)
        cols_remain=np.append(np.arange((kx - kx_gap)/2., dtype=np.int64),
                              np.flip(kx-1-np.arange((kx - kx_gap)/2., dtype=np.int64)))
        idx_x = idx_x[:, cols_remain]


    idx_x[idx_x < 0]=0
    idx_x[idx_x > n_x-1]=n_x-1
    idx_y[idx_y < 0]=0
    idx_y[idx_y > n_y-1]=n_y-1

    out_img = np.zeros(array_ext.shape)
    for iy in np.arange(n_y):
        med_arr = np.nanmedian(array_ext[idx_y[iy,0]:idx_y[iy,-1]+1,idx_x], axis=[0,2])
        out_img[iy,:] = med_arr
        
    return out_img[1:-1,1:-1]  


def median_filter_2dKernel(array, 
                           kernel, 
                           xcen=None,
                           ycen=None,
                           preserve_nan=False):

    if xcen==None:
        xcen = kernel.shape[1]//2
    if ycen==None:
        ycen = kernel.shape[0]//2
        
    array = np.array(array)
    ny, nx = array.shape
    ky, kx = kernel.shape
    nxe = nx + kx*2
    nye = ny + ky*2
    array_ext = np.zeros((nye,nxe))+np.nan
    array_ext[ky:ky+ny, kx:kx+nx] = array
    
    outimg = np.zeros(array.shape)+np.nan
    print('outimg.shape: ', outimg.shape)
    for ix in range(nx):
        for iy in range(ny):
            arr = array_ext[ky+iy-ycen:ky+iy-ycen+ky,
                            kx+ix-xcen:kx+ix-xcen+kx]
            
            outimg[iy,ix] = np.nanmedian(arr[kernel])
            
    if preserve_nan:
        outimg[np.isnan(array)]=np.nan
    return outimg
        
