import numpy as np
from scipy import interpolate
from scipy.signal import find_peaks, medfilt
from sklearn import linear_model
import matplotlib.pyplot as plt
from eiger.QuasarSpec.loadQsoSpec import loadQsoSpec
import os
from astropy.table import Table

def redsideContin(specdict, instrument, writefile=False):

    spec = specdict[instrument]
    redshift = specdict['z_em']
    
    wave     = spec['wave']
    restwave = wave / (1+redshift)
    flux     = spec['flux']
    ivar     = spec['ivar']
    naxis1   = len(flux)

    fitmask  = np.ones(len(wave),dtype=np.bool8)

    #############
    #
    # Mask out known problem areas
    fitmask[restwave < 1210] = False
    fitmask[np.logical_and(wave > 13590, wave < 14100)] = False
    fitmask[np.logical_and(wave > 17900, wave < 19400)] = False

    if (instrument=='MOSFIRE_Y'):
        fitmask[np.logical_or(wave < 9700, wave > 10800)] = False
    if (instrument=='MOSFIRE_J'):
        fitmask[np.logical_or(wave < 11500, wave > 13550)] = False
    if (instrument=='MOSFIRE_H'):
        fitmask[np.logical_or(wave < 14700, wave > 17900)] = False
    if (instrument=='MOSFIRE_K'):
        fitmask[wave < 19520] = False
    #
    #####################
        
    # 0. Divide by a smoothed spectrum to take out large departures from
    #    the sensitivity functions

    binwidth  = 150
    halfwidth = binwidth // 2
    nbins     = (naxis1 // binwidth) - 1
    indx = halfwidth + np.arange(nbins) * binwidth
    gd = (fitmask[indx] == True)
    xbin = restwave[indx[gd]]

    # Note this intentionally fits high (70 %-ile)
    # to eliminate dipping into absorption lines
    ybin = np.array([np.quantile(flux[ii-halfwidth:ii+halfwidth],0.70) for ii in indx[gd]])

    # Screen bad areas and save points for later
    ybin[xbin < 1216] = ybin[xbin > 1216][0]
    ybin[ybin < 1e-20] = np.median(ybin)
    xbin1 = xbin
    ybin1 = ybin

    # Fit the knots with a bspline and normalize
    t, c, k = interpolate.splrep(xbin, ybin, k=3, s=0)
    flux_filt = interpolate.BSpline(t,c,k)
    flux = flux / flux_filt(restwave)

    ##################################

    binwidth  = 200
    halfwidth = binwidth // 2
    nbins     = (naxis1 // binwidth) - 1
    indx = halfwidth + np.arange(nbins) * binwidth
    gd = (fitmask[indx] == True)
    xbin = restwave[indx[gd]]
    
    # 1. Calculate an upper envelope

    ybin = [np.quantile(flux[ii-halfwidth:ii+halfwidth],0.75) for ii in indx[gd]]
    flux_filt2 = interpolate.interp1d(xbin, ybin, bounds_error=False, \
                                     fill_value='extrapolate')

    pk_mask, _ = find_peaks(flux,height=(flux_filt2(restwave),None))

    # With successive iterations, tighten the criteria to keep points
    if (instrument == 'MOSFIRE'):
        niter=3
        discard = [5,4,3]
    elif (instrument == 'FIRE' or instrument == 'XSH_NIR'):
        niter=3
        discard = [4,3,2]
    elif (instrument == 'XSH_VIS'):
        niter=3
        discard = [4,3,2]
    else:
        niter=2
        discard = [4,3]

    bspline = True

    # Iteratively reject outliers while fitting a bspline continuum
    for i in range(niter):
    
        if (i == 0):
            f = interpolate.interp1d(restwave[pk_mask],flux[pk_mask],\
                                     bounds_error=False, \
                                     fill_value='extrapolate')
        else:
            f = c
        
        # 3. Apply RANSAC to detect outlying features (esp. absorption) and mask

        mad = np.median(np.abs(np.median(flux)-flux))

        ransac = linear_model.RANSACRegressor(random_state=0, \
                                              loss='absolute_loss',\
                                              residual_threshold=discard[i]*mad)

        ransac.fit(restwave.reshape(len(restwave),1), flux, \
                   sample_weight=np.abs(ivar))

        inlier_mask  = ransac.inlier_mask_
        inlier_mask[np.logical_and(restwave > 1215, restwave < 1230)] = True
        inlier_mask[np.logical_and(wave > 13590, wave < 14100)] = False
        inlier_mask[np.logical_and(wave > 17900, wave < 19400)] = False
        inlier_mask[restwave < 1215] = False

        if (instrument=='MOSFIRE_Y'):
            inlier_mask[np.logical_or(wave < 9700, wave > 10800)] = False
        if (instrument=='MOSFIRE_J'):
            inlier_mask[np.logical_or(wave < 11500, wave > 13550)] = False
        if (instrument=='MOSFIRE_H'):
            inlier_mask[np.logical_or(wave < 14700, wave > 17900)] = False
        if (instrument=='MOSFIRE_K'):
            inlier_mask[wave < 19520] = False
        
        # 4. Smooth the inlier data points to build the continuum

        if (bspline):
            gd = (fitmask[indx] == True)
            ybin = [np.median(flux[ii-halfwidth:ii+halfwidth]) for ii in indx[gd]]
            t, c, k = interpolate.splrep(xbin, \
                                         ybin, \
                                         k=3, s=0)
            c = interpolate.BSpline(t,c,k)
            
        else:
            c = interpolate.interp1d(restwave[inlier_mask],\
                                     medfilt(flux[inlier_mask],kernel_size=51), \
                                     fill_value='extrapolate')

        cont = c(restwave)

    cont = cont * flux_filt(restwave)
    cont[restwave < 1216] = (cont[restwave > 1216])[0]
    flux = flux * flux_filt(restwave)
        
    if (True):
        plt.step(wave,flux,where='mid')
        plt.plot(wave,cont,color='r')
#        plt.plot(wave,flux_filt(restwave))
#        plt.plot(xbin1 * (1+redshift),ybin1,'+')
#        plt.plot(wave[inlier_mask==False],flux[inlier_mask==False], '+')
        plt.ylim(-1.0*np.median(cont),3*np.median(cont))
        plt.show()

    # Write the output to disk as a fits binary table
    # in the same directory as the data file
    if (writefile):
        infile   = loadQsoSpec(specdict['objid'], files=True)
        path     = os.path.dirname(infile[instrument])
        contname = os.path.basename(infile[instrument])[:-5]+'_contin.fits'
        outfile  = path+'/'+contname
        t = Table([wave,cont,inlier_mask],names=('wave','cont','mask'))
        t.write(outfile, format='fits', overwrite=True)
        print(f"Writing: {outfile}")
        
    return(cont)


def allContin(writefile=False):

    # 1 = SDSS1030
    qso       = loadQsoSpec(1)
    firecontin      = redsideContin(qso,'FIRE',writefile=writefile)
    xsh_nir_contin  = redsideContin(qso,'XSH_NIR',writefile=writefile)
    xsh_vis_contin  = redsideContin(qso,'XSH_VIS',writefile=writefile)
    hirescontin     = redsideContin(qso,'HIRES',writefile=writefile)

    # 2 = P159-02
    qso = loadQsoSpec(2)
    firecontin      = redsideContin(qso,'FIRE',writefile=writefile)
    xsh_nir_contin  = redsideContin(qso,'XSH_NIR',writefile=writefile)
    xsh_vis_contin  = redsideContin(qso,'XSH_VIS',writefile=writefile)
    
    # 3 = 1120+0641
    qso = loadQsoSpec(3)
    firecontin      = redsideContin(qso,'FIRE',writefile=writefile)
    xsh_nir_contin  = redsideContin(qso,'XSH_NIR',writefile=writefile)
    xsh_vis_contin  = redsideContin(qso,'XSH_VIS',writefile=writefile)

    # 4 = J0100+2802
    qso = loadQsoSpec(4)
    firecontin  = redsideContin(qso,'FIRE',writefile=writefile)
    xsh_nir_contin  = redsideContin(qso,'XSH_NIR',writefile=writefile)
    xsh_vis_contin  = redsideContin(qso,'XSH_VIS',writefile=writefile)
    # HIRES spectrum already normalized
    # hires_contin  = redsideContin(qso,'HIRES')
    hires_contin = np.ones(len(qso['HIRES']['flux']))
    hires_mask   = np.ones(len(qso['HIRES']['flux']),dtype=np.bool8)
    hires_mask[qso['HIRES']['wave'] < 1215] = False
    if (writefile):
        infile   = loadQsoSpec(4,files=True)['HIRES']
        path     = os.path.dirname(infile)
        contname = os.path.basename(infile)[:-5]+'_contin.fits'
        outfile  = path+'/'+contname
        t = Table([qso['HIRES']['wave'],hires_contin,hires_mask],\
                  names=('wave','cont','mask'))
        t.write(outfile,format='fits',overwrite=True)
        print(f"writing: {outfile}")
                   
    # 5 = 1148+5251
    qso = loadQsoSpec(5)
    mosfire_y  = redsideContin(qso,'MOSFIRE_Y',writefile=writefile)
    mosfire_j  = redsideContin(qso,'MOSFIRE_J',writefile=writefile)
    mosfire_h  = redsideContin(qso,'MOSFIRE_H',writefile=writefile)
    mosfire_k  = redsideContin(qso,'MOSFIRE_K',writefile=writefile)
    # HIRES spectrum already normalized
    # hires_contin  = redsideContin(qso,'HIRES')
    hires_contin = np.ones(len(qso['HIRES']['flux']))
    hires_mask   = np.ones(len(qso['HIRES']['flux']),dtype=np.bool8)
    hires_mask[qso['HIRES']['wave'] < 1215] = False
    if (writefile):
        infile   = loadQsoSpec(4,files=True)['HIRES']
        path     = os.path.dirname(infile)
        contname = os.path.basename(infile)[:-5]+'_contin.fits'
        outfile  = path+'/'+contname
        t = Table([qso['HIRES']['wave'],hires_contin,hires_mask],\
                  names=('wave','cont','mask'))
        t.write(outfile,format='fits',overwrite=True)
        print(f"writing: {outfile}")
    
    # 6 = J0148+0600
    qso = loadQsoSpec(6)
    xsh_nir_contin  = redsideContin(qso,'XSH_NIR',writefile=writefile)
    xsh_vis_contin  = redsideContin(qso,'XSH_VIS',writefile=writefile)

