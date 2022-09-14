import os
from astropy.nddata import Cutout2D
from astropy.visualization import simple_norm
from astropy.wcs import WCS, utils
from astropy.io import fits
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
import grismconf
from jwst import datamodels


def xy_drizzled_to_xy_rate(emfile,drizzlefile,x_drizzled,y_drizzled):
    #the emfile is the single emission-line file (from a _cal.fits) corresponding to a rate.fits 
    #the outout are the pixel positions in the original rate file
    hdul = fits.open(emfile)
    wcs_em=WCS(hdul['SCI'].header) 

    hdu2=fits.open(drizzlefile)
    wcs_drizzle=WCS(hdu2['SCI'].header)

    ra_pix,dec_pix=wcs_drizzle.all_pix2world(x_drizzled,y_drizzled,0)
    x_rate,y_rate=wcs_em.all_world2pix(ra_pix,dec_pix,0)

    return x_rate,y_rate

def xylamb_dispersed_to_radec_wcs(inpwcs, 
                                  x_line, 
                                  y_line, 
                                  obs_lamb_line, 
                                  filt='F356W', module='A', version='V4'):

    
    #configfile='/Users/kashino/Studies/JWST/mirage2/my_jwreftools_execution_V2/NIRCAM_F356W_modA_R.conf'):
    #delta_lamb_int =0.0144 O3_5008 to Hbeta
    #ratefile='jw01243001001_01101_00001_nrca5_flatfieldstep.fits'
    #grism_with_wcs = datamodels.open(ratefile)

    if module.upper()!='A' and module.upper()!='B':
        raise ValueError('module must be A or B.')

    # if version==2:
    #     configfile='/Users/kashino/Studies/JWST/mirage2/referenceFiles/mirage_data/nircam/GRISM_NIRCAM/GRISM_NIRCAM/V2/NIRCAM_'+filter+'_mod'+module.upper()+'_R.conf'
    # elif version==4:
    #     configfile='/Users/kashino/Studies/JWST/GRISM_NIRCAM/V4/NIRCAM_'+filter+'_mod'+module.upper()+'_R.conf'
    # else:
    #     raise ValueError('GRISM version musht be either 2 or 4.')
    condir = '/scratch/kashinod/JWST/GRISM_NIRCAM/'+version
    configfile = os.path.join(condir, 'NIRCAM_'+filt+'_mod'+module.upper()+'_R.conf')
    
    pix_to_world = inpwcs.get_transform('detector','world')

    C = grismconf.Config(configfile)

    xref=x_line
    yref=y_line
    
    t= C.INVDISPL('+1',xref,yref,obs_lamb_line)
    dx=C.DISPX('+1',xref,yref,t) 
    dy=C.DISPY('+1',xref,yref,t) 
    x0=x_line-dx
    y0=y_line-dy
    #RA,DEC,a,b=pix_to_world(x0,y0,0,0)
    RA,DEC=pix_to_world(x0,y0)
    
    return RA,DEC

def radeclamb_to_xy_dispersed(inpwcs, ra, dec, obs_lamb_line, filter='F356W', module='A', version='V4'):
    #configfile='/Users/kashino/Studies/JWST/mirage2/my_jwreftools_execution_V2/NIRCAM_F356W_modA_R.conf'):
    #delta_lamb_int =0.0144 O3_5008 to Hbeta
    #ratefile='jw01243001001_01101_00001_nrca5_flatfieldstep.fits'
    #grism_with_wcs = datamodels.open(ratefile)

    if module.upper()!='A' and module.upper()!='B':
        raise ValueError('module must be A or B.')
    
    # if version==2:
    #     configfile='/Users/kashino/Studies/JWST/mirage2/referenceFiles/mirage_data/nircam/GRISM_NIRCAM/GRISM_NIRCAM/V2/NIRCAM_'+filter+'_mod'+module.upper()+'_R.conf'
    # elif version==4:
    #     configfile='/Users/kashino/Studies/JWST/GRISM_NIRCAM/V4/NIRCAM_'+filter+'_mod'+module.upper()+'_R.conf'
    # else:
    #     raise ValueError('GRISM version musht be either 2 or 4.')
    #condir = '/Users/kashino/Studies/JWST/GRISM_NIRCAM/'+version
    condir = '/scratch/kashinod/JWST/GRISM_NIRCAM/'+version
    configfile = os.path.join(condir, 'NIRCAM_'+filter+'_mod'+module.upper()+'_R.conf')
    
    C = grismconf.Config(configfile)    

    use_WCS_world_to_pixel=False
    if use_WCS_world_to_pixel:
        print('use world_to_pixel.', flush=True)
        world_to_pix = inpwcs.world_to_pixel
        x0,y0=world_to_pix(SkyCoord(ra=ra, dec=dec, unit='deg'))
    else:
        world_to_pix = inpwcs.get_transform('world','detector')
        x0,y0=world_to_pix(ra, dec)
            
    xref=x0 #2048.
    yref=y0 #2048.

    t= C.INVDISPL('+1',xref,yref,obs_lamb_line)
    dx=C.DISPX('+1',xref,yref,t) 
    dy=C.DISPY('+1',xref,yref,t) 
    #RA,DEC,a,b=pix_to_world(x0,y0,0,0)

    x_line = x0 + dx
    y_line = y0 + dy
    return x_line,y_line


def radeclamb_to_xy_coadded_dispersed(inpwcs, 
                                      ra,
                                      dec,
                                      obs_lamb_line,
                                      cal_wcs,
                                      i2d_wcs,
                                      **kwargs):

    x_cal, y_cal = radeclamb_to_xy_dispersed(
        inpwcs, ra, dec, obs_lamb_line, **kwargs)
    sky = cal_wcs.pixel_to_world(x_cal, y_cal)
    x, y = i2d_wcs.world_to_pixel(sky)
    return x, y


def wavelengths():
    w_4341  = 4341.692     * u.angstrom
    w_4862  = 4862.69      * u.angstrom
    w_4960  = 4960.295     * u.angstrom
    w_5008  = 5008.240     * u.angstrom
    w_6564  = 6564.61      * u.angstrom
    w_6585  = 6585.27      * u.angstrom
    w_6718  = 6718.294     * u.angstrom
    w_6732  = 6732.673     * u.angstrom
    w_9071  = 9071.1       * u.angstrom
    w_9533  = 9533.2       * u.angstrom
    w_10052 = 10052.123    * u.angstrom         # P-delta
    w_10830 = 10833.306444 * u.angstrom         # HeI
    w_10941 = 10941.082    * u.angstrom         # P-gamma
    w_12570 = 12570.2068   * u.angstrom         # [FeII] 12570
    w_12821 = 12821.576    * u.angstrom         # P-beta
    w_18756 = 18756.096    * u.angstrom         # P-alpha
    w_19450 = 19450.89     * u.angstrom         # Br-delta
    w_21661 = 21661.21     * u.angstrom         # Br-gamma
    w_26258 = 26258.68     * u.angstrom         # Br-beta

    w_lines = {'4341': w_4341,
               '4862' :w_4862, 
               '4960' :w_4960, 
               '5008' :w_5008,
               '6564' :w_6564, 
               '6585' :w_6585, 
               '6718' :w_6718, 
               '6732' :w_6732, 
               '9071' :w_9071, 
               '9533' :w_9533, 
               '10052':w_10052,
               '10830':w_10830,
               '10941':w_10941,
               '12570':w_12570,
               '12821':w_12821,
               '18756':w_18756,
               '19450':w_19450,
               '19450':w_19450,
               '21661':w_21661,
               '26258':w_26258}
    return w_lines





import numpy as np
from scipy.interpolate import interpn
from astropy.convolution import convolve, Tophat2DKernel

def get_tadpole_dxdy(src_x, src_y, tp_model):

    fill_value=None
    src_xy = np.array([src_y, src_x]).T
    tp1_dx = interpn((tp_model['GRISM_Y_IND'],tp_model['GRISM_X_IND']), 
                     tp_model['DELTA_X_L'],
                     src_xy,
                     method='linear', bounds_error=False,
                     fill_value=fill_value)
    tp1_dy = interpn((tp_model['GRISM_Y_IND'],tp_model['GRISM_X_IND']), 
                     tp_model['DELTA_Y_L'], 
                     src_xy,
                     method='linear', bounds_error=False, 
                     fill_value=fill_value)
    tp2_dx = interpn((tp_model['GRISM_Y_IND'],tp_model['GRISM_X_IND']), 
                     tp_model['DELTA_X_R'], 
                     src_xy,
                     method='linear', bounds_error=False, 
                     fill_value=fill_value)
    tp2_dy = interpn((tp_model['GRISM_Y_IND'],tp_model['GRISM_X_IND']), 
                     tp_model['DELTA_Y_R'], 
                     src_xy,
                     method='linear', bounds_error=False, 
                     fill_value=fill_value)
    return tp1_dx, tp1_dy, tp2_dx, tp2_dy

def get_tadpole_mask_for_xy(src_x, src_y, wcs, tp_model, rpix=10., src_x_min = -220, src_y_min=-200):

    nx = 2048 ## for LW cal.fits
    ny = 2048 ## for LW cal.fits

    mask_nx = nx+100
    mask_ny = ny+100

    tp1_dx, tp1_dy, tp2_dx, tp2_dy = get_tadpole_dxdy(src_x, src_y, tp_model)

    tp1_x = src_x + tp1_dx
    tp1_y = src_y + tp1_dy
    tp2_x = src_x + tp2_dx
    tp2_y = src_y + tp2_dy

    idx_mask = np.where((src_x>src_x_min)&
                        (src_y>src_y_min)&
                        (tp1_x<mask_nx)&
                        (tp2_y<mask_ny))[0]
    if idx_mask.size==0:
        return np.ones((ny, nx), dtype=int)



    ### Tophat kernel
    tophatk = Tophat2DKernel(rpix)
    
    mask_arr = np.ones((idx_mask.size, mask_ny, mask_nx))

    print('idx_mask.size: ', idx_mask.size, flush=True)
    for i in range(idx_mask.size):
        i_src = idx_mask[i]

        xx = np.arange(np.int32(tp1_x[i_src]),
                       np.int32(tp2_x[i_src]+1)+1)
        yy = np.int32(np.interp(xx, 
                                np.array([tp1_x[i_src], tp2_x[i_src]]),
                                np.array([tp1_y[i_src], tp2_y[i_src]])))
        tmp = np.where((xx>=0)&(xx<mask_nx)&(yy>=0)&(yy<mask_ny))[0]
        if tmp.size>0:
            mask_arr[i, yy[tmp],xx[tmp]]=0            
            x0 = np.int(tp1_x[i_src] - 2*rpix)
            x1 = np.int(tp2_x[i_src] + 2*rpix)
            y0 = np.int(tp2_y[i_src] - 2*rpix)
            y1 = np.int(tp1_y[i_src] + 2*rpix)
            x0 = np.max([0, x0])
            x1 = np.min([mask_nx-1, x1])
            y0 = np.max([0, y0])
            y1 = np.min([mask_ny-1, y1])
            print(x0,x1,y0,y1, flush=True)
            mask_arr[i, y0:y1,x0:x1] = convolve(mask_arr[i, y0:y1,x0:x1], tophatk, boundary='fill', fill_value=1)
            mask_arr[i, y0:y1,x0:x1][np.where(mask_arr[i, y0:y1,x0:x1]<0.999999)]=np.nan
            mask_arr[i, y0:y1,x0:x1][np.where(mask_arr[i, y0:y1,x0:x1]>0.999999)]=1.0

    mask = np.prod(mask_arr[:, 0:ny, 0:nx], axis=0)
    return mask

    
    
def get_tadpole_mask_for_sky(skycoord, wcs, tp_model, rpix=10., src_x_min = -220, src_y_min=-200, src_xoff=0., src_yoff=0.):

    src_x, src_y = wcs.world_to_pixel(skycoord)

    ### Arbitrary offset (optional)
    src_x = src_x + src_xoff
    src_y = src_y + src_yoff

    mask = get_tadpole_mask_for_xy(src_x, src_y, wcs, tp_model, rpix=rpix, src_x_min=src_x_min, src_y_min=src_y_min)
    return mask


    



    


#===========================

def myextract(ra0, dec0, hdul, grism_with_wcs, yoffset=0., yhsize=21, use_wcs=False, wcs=None):
    
    # Get header
    h = hdul[0].header
    
    # Find some information about the observing mode of this rate file.
    filt = h["FILTER"] # Filter name, e.g. F410M
    grism = h["PUPIL"][-1] # R or C
    module = h["MODULE"] # Which NIRCAM module, A or B
    print("Filter:",filt)
    print("grism:",grism)
    print("Module:",module)
    
    import grismconf
    condir = '/scratch/kashinod/JWST/GRISM_NIRCAM/V4'
    configfile = os.path.join(condir, 'NIRCAM_'+filt+'_mod'+module.upper()+'_'+grism+'.conf')
    C = grismconf.Config(configfile)

    # We initialize a datamodel WCS to the data
    # grism_with_wcs = datamodels.open(dataname)
    world_to_pix = grism_with_wcs.meta.wcs.get_transform('world','detector')
    #pix_to_world = grism_with_wcs.meta.wcs.get_transform('detector','world')
    x0, y0 = world_to_pix(ra0,dec0)
    if use_wcs:
        x0, y0 = wcs.world_to_pixel(SkyCoord(ra=ra0, dec=dec0, unit='degree'))
        
        
    # Compute the position of the source in the image in  pixel coordinates
    
    y0 = y0 + yoffset
    #print(dataname,"Target is at ",x0,y0)

    xref=x0 #2048.
    yref=y0 #2048.

    t= C.INVDISPL('+1',xref,yref,np.array([2.99, 4.21]))
    dx=C.DISPX('+1',xref,yref,t) 
    dy=C.DISPY('+1',xref,yref,t) 
    x_line = x0 + dx
    y_line = y0 + dy
    minx0 = np.max([0, np.int32(np.min(x_line))])
    maxx0 = np.min([2047, np.int32(np.max(x_line))])
    miny0 = np.max([0, np.int32(np.min(y_line-yhsize))])
    maxy0 = np.min([2047, np.int32(np.max(y_line+yhsize))])
    
    if minx0>2047 or maxx0 < 0 or miny0>2047 or maxy0 < 0:
        print('ERROR: Spectrum is out of fiels.')
        print('--- minx0, maxx0: ', minx0, maxx0)
        print('--- miny0, maxy0: ', miny0, maxy0)
        return None
        
    # We load the observation, extimate and subtract the data, masking sources using the Mirage simulation
    data = hdul["SCI"].data  ## MJy/sr
    err = hdul["ERR"].data 
    dq = hdul["DQ"].data
  
    # We trim our data to be the stamp containing the spectrum we want to extract
    # Since we use the minx0,maxx0,miny0,maxy0 coordinates from our simulation, this stamp should
    # be just like our simulated data but containing the actual observation.
    data = data[miny0:maxy0+1,minx0:maxx0+1]
    err = err[miny0:maxy0+1,minx0:maxx0+1]
    dq = dq[miny0:maxy0+1,minx0:maxx0+1]
    
    # We estimate of the contamination, from our Mirage simulation of the whole field, as noted
    # above, contamination is just models-model0
    #contam = models[miny0:maxy0+1,minx0:maxx0+1] - model0

    # We load trhe grismconf module and load the proper configuration for this NIRCAM mode
    # Description of what grismconf does and how grism calibration  is done is available in 
    # WFC3 ISR 2017-01 from Pirzkal et al. (http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2017-01.pdf)
    
    # These are the coordinates of all the pixels in our 2D stamp, but in the full image (wrt calibration is known)
    ys,xs = np.indices((maxy0-miny0+1,maxx0-minx0+1))
    # xs and ys are now the relatice dx and dy offsets from the position of our source
    # They are both 2D arrays of x and y coordinates.
    xs = xs + minx0 - x0
    ys = ys + miny0 - y0

    # Depending on whether the grism disperse in the x or y direction, we use the INVDISPX or INVDISPY functions
    # to compute the value for t for every pixel in our 2D stamps
    if grism=="R":
        ts = C.INVDISPX("+1",x0,y0,xs)
        dys = C.DISPY("+1",x0,y0,ts) + ys -  2*C.DISPY("+1",x0,y0,ts)
        
    if grism=="C":
        ts = C.INVDISPY("+1",x0,y0,ys)
        dys = C.DISPX("+1",x0,y0,ts) + xs
    
    # Now compute the wavelength of every pixel in our 2D stamp
    ws = C.DISPL("+1",x0,y0,ts)

    # Now, depending of whether things are in the row or col, we transpose things so that we can look at them 
    # properly (i.e. row direction)
    if grism=="C":
        m = np.transpose(model0) # The model counts in each pixel
        l = np.transpose(ws)     # The wavelength of each pixel
        d = np.transpose(data)   # The data counts in each pixel
        c = np.transpose(contam) # The contamination estimated counts in each pixel
        e = np.transpose(err)    # THe data error estimates in each pixel
        q = np.transpose(dq)     # The data DQ in each pixel
        y = np.transpose(dys)    # The cross-dispersion distance of each pixel from the trace

    if grism=="R":
#        m = model0
        l = ws
        d = data
#        c = contam
        e = err
        q = dq
        y = dys

    # We now create optimal extraction weights, using our simulation (which accounts for the object profile etc)
    # The extraxction weights are based on the normalized simulated spectral profile at each wavelength
    # Sum up the model in the y-direction and replicate that 
 #   ysum = np.sum(m,axis=0) 
 #   w = np.repeat([ysum], np.shape(m)[0], axis=0)
 #   weight = m/w
    
    # We make sure our original data only has valid data and return everything we computed
    # What we return are no longer 2D stamps but 1D vectors
    #ok = np.isfinite(d)
    
#    return  d[ok], e[ok], q[ok], l[ok], c[ok], m[ok], weight[ok], y[ok], C
    #return  d[ok], e[ok], q[ok], l[ok], 0,      0,     0, y[ok], C
    return  d, e, q, l, 0,      0,     0, y, C, filt, grism, module

def shift_columns(y,d,var,lamb,q,ysize=31):
    
    ly,lx=np.shape(d)
    shifted_y=np.zeros((ysize,lx))
    shifted_y_var=np.zeros((ysize,lx))
    shifted_y_lamb=np.zeros((ysize,lx))

    for ix in range(lx):
        flx_col = d[:, ix]
        var_col = var[:,ix]
        lamb_col = lamb[:,ix]
        q_col=q[:,ix]

        bad_col = var_col*0.
        bad_col[np.where((var_col<=0.0)|(np.isnan(var_col))|(np.isnan(flx_col)))[0]]=1.0
        bad_col[np.mod(q_col,2)==1]=1.0

        ### Cumulative sum
        flx_col_cum = np.nancumsum(flx_col)
        bad_col_cum = np.nancumsum(bad_col)
        var_col_cum = np.nancumsum(var_col)
        lamb_col_cum = np.nancumsum(lamb_col)


        ### New axis
        yaxis=y[:,ix]
        #yaxis_new = np.arange(ysize)-0.5-ysize/2.  #y_trace_xaxis[ix]-15+np.arange(31)
        yaxis_new = np.arange(ysize) - np.int(ysize/2.)

        flx_col_new0 = np.interp(yaxis_new - 0.5, yaxis+0.5, flx_col_cum, left=np.nan, right=np.nan)
        flx_col_new1 = np.interp(yaxis_new + 0.5, yaxis+0.5, flx_col_cum, left=np.nan, right=np.nan)
        var_col_new0 = np.interp(yaxis_new - 0.5, yaxis+0.5, var_col_cum, left=np.nan, right=np.nan)
        var_col_new1 = np.interp(yaxis_new + 0.5, yaxis+0.5, var_col_cum, left=np.nan, right=np.nan)
        bad_col_new0 = np.interp(yaxis_new - 0.5, yaxis+0.5, bad_col_cum, left=np.nan, right=np.nan)
        bad_col_new1 = np.interp(yaxis_new + 0.5, yaxis+0.5, bad_col_cum, left=np.nan, right=np.nan)

        lamb_col_new0 = np.interp(yaxis_new - 0.5, yaxis+0.5, lamb_col_cum, left=np.nan, right=np.nan)
        lamb_col_new1 = np.interp(yaxis_new + 0.5, yaxis+0.5, lamb_col_cum, left=np.nan, right=np.nan)


        lamb_col_new = lamb_col_new1-lamb_col_new0

        bad_col_new = bad_col_new1-bad_col_new0
        flx_col_new = bad_col_new * 0.
        var_col_new = bad_col_new * 0.

        idx_bad  = np.where(bad_col_new>=0.5)[0]
        idx_good = np.where(bad_col_new<0.5)[0]

        bad_col_new[idx_bad]=1.0


        flx_col_new[idx_good] = (flx_col_new1-flx_col_new0)[idx_good]/(1-bad_col_new[idx_good])
        flx_col_new[idx_bad] = np.nan
        var_col_new[idx_good] = (var_col_new1-var_col_new0)[idx_good]/(1-bad_col_new[idx_good])
        var_col_new[idx_bad] = np.nan

        #flx_col_new = (flx_col_new1-flx_col_new0)#[idx_good]/(1-bad_col_new[idx_good])


        shifted_y[:,ix] = flx_col_new
        shifted_y_var[:,ix]=var_col_new
        shifted_y_lamb[:,ix]=lamb_col_new

    return shifted_y,shifted_y_var,shifted_y_lamb


def scrunch_columns(lambda_array,l,shifted_y,shifted_y_var,module,flambda=False):

    ## shifted_y in the unit of DN/s (count rate)
    
    ly,lx=np.shape(shifted_y)
    len_lamb=len(lambda_array)
    scrunched_x=np.zeros((ly,len_lamb))
    scrunched_x_var=np.zeros((ly,len_lamb))
    scrunched_lamb=np.zeros((ly,len_lamb))

    dwav=lambda_array[1]-lambda_array[0]


    for iy in range(ly):
        lambaxis = l[iy,:]
        flx_col = shifted_y[iy,:]
        var_col = shifted_y_var[iy,:]


        if module.upper()=='B':
            lambaxis=lambaxis[::-1]
            axis=lambaxis[::-1]
            flx_col=flx_col[::-1]
            var_col=var_col[::-1]

        lxi = np.arange(lambaxis.size)
        lambgrid = np.interp(lxi-0.5, lxi, lambaxis)
        lambgrid[0]=lambgrid[1]-(lambgrid[2]-lambgrid[1])
        lambgrid = np.append(lambgrid,lambgrid[-1]+(lambgrid[-1]-lambgrid[-2]))
        

        if flambda==False:
            dlambdadx = (lambgrid[1:]-lambgrid[0:-1])*1e4
            flx_col = flx_col * dlambdadx ## (DN/s) --> (DN/s)/angstrom
            var_col = var_col * dlambdadx ## (DN/s)^2 --> (DN/s)^2/angstrom
        
        bad_col = var_col*0.
        bad_col[np.where((var_col<=0.0)|(np.isnan(var_col))|(np.isnan(flx_col)))[0]]=1.0

        ### Cumulative sum
        flx_row_cum = np.append(0,np.nancumsum(flx_col))
        bad_row_cum = np.append(0,np.nancumsum(bad_col))
        var_row_cum = np.append(0,np.nancumsum(var_col))

  
        #print(np.flip(lambaxis)+0.5)
        #print((lambda_array-dwav/2.),np.flip(lambaxis)+0.5)
        flx_row_new0 = np.interp(lambda_array-dwav/2., lambgrid, flx_row_cum, left=np.nan, right=np.nan)
        flx_row_new1 = np.interp(lambda_array+dwav/2., lambgrid, flx_row_cum, left=np.nan, right=np.nan)
        var_row_new0 = np.interp(lambda_array-dwav/2., lambgrid, var_row_cum, left=np.nan, right=np.nan)
        var_row_new1 = np.interp(lambda_array+dwav/2., lambgrid, var_row_cum, left=np.nan, right=np.nan)
        bad_row_new0 = np.interp(lambda_array-dwav/2., lambgrid, bad_row_cum, left=np.nan, right=np.nan)
        bad_row_new1 = np.interp(lambda_array+dwav/2., lambgrid, bad_row_cum, left=np.nan, right=np.nan)

        bad_row_new = bad_row_new1-bad_row_new0
        flx_row_new = bad_row_new * 0.
        var_row_new = bad_row_new * 0.

        idx_bad  = np.where(bad_row_new>=0.5)[0]
        idx_good = np.where(bad_row_new<0.5)[0]

        bad_row_new[idx_bad]=1.0

        flx_row_new[idx_good] = (flx_row_new1-flx_row_new0)[idx_good]/(dwav*1e4)/(1-bad_row_new[idx_good])
        flx_row_new[idx_bad] = np.nan

        var_row_new[idx_good] = (var_row_new1-var_row_new0)[idx_good]/(dwav*1e4)/(1-bad_row_new[idx_good])
        var_row_new[idx_bad] = np.nan


        scrunched_x[iy,:] = flx_row_new
        scrunched_x_var[iy,:] = var_row_new
        scrunched_lamb[iy,:]=lambda_array
    return scrunched_x,scrunched_x_var,scrunched_lamb



def mygauss(x, p0,p1,p2,p3,p4):
    xx = np.power((x-p1)/p2,2)
    res = p0 * np.exp(-xx/2.) + (p3 + p4*(x-p1))
    return res





