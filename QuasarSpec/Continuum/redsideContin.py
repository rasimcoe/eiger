import numpy as np
from scipy import interpolate
from scipy.signal import find_peaks, medfilt
from sklearn import linear_model
import matplotlib.pyplot as plt
from eiger.QuasarSpec.loadQsoSpec import loadQsoSpec
import os
from astropy.table import Table

def redsideContin(specdict, instrument, writefile=False, bspline=True):

    spec = specdict[instrument]
    redshift = specdict['z_em']
    
    wave     = spec['wave']
    restwave = wave / (1+redshift)
    all_wave      = spec['wave']
    all_restwave  = all_wave / (1+redshift)
    flux     = spec['flux']
    ivar     = spec['ivar']
    naxis1   = len(flux)

    fitmask  = np.ones(len(wave),dtype=np.bool8)

    #############
    #
    # Mask out known problem areas
    #
    fitmask[restwave < 1210] = False
    fitmask[np.logical_and(wave > 13490, wave < 14300)] = False
    fitmask[np.logical_and(wave > 17940, wave < 19500)] = False

    if (instrument=='MOSFIRE_Y'):
        fitmask[np.logical_or(wave < 9580, wave > 11400)] = False
    if (instrument=='MOSFIRE_J'):
        fitmask[np.logical_or(wave < 11350, wave > 13750)] = False
    if (instrument=='MOSFIRE_H'):
        fitmask[np.logical_or(wave < 14750, wave > 17950)] = False
    if (instrument=='MOSFIRE_K'):
        fitmask[wave < 19520] = False

    if (specdict['objname']=='J1030+0524' or specdict['objid']==1):
        fitmask[np.logical_and(wave > 8990,wave < 9024)] = False
        fitmask[np.logical_and(wave > 15600,wave < 15700)] = False
        fitmask[np.logical_and(wave > 10540,wave < 10620)] = False
    if (specdict['objname']=='P159-02' or specdict['objid']==2):
        fitmask[np.logical_and(wave > 13040,wave < 13100)] = False
    if (specdict['objname']=='J0100+2802' or specdict['objid']==4):
        fitmask[np.logical_and(wave > 14600,wave < 14660)] = False
        fitmask[np.logical_and(wave > 9450,wave < 9500)] = False
        fitmask[np.logical_and(wave > 9430,wave < 9440)] = False
        fitmask[np.logical_and(wave > 9215,wave < 9225)] = False
        fitmask[np.logical_and(wave > 9257,wave < 9270)] = False
        fitmask[np.logical_and(wave > 9800,wave < 9850)] = False
    if (specdict['objname']=='J0148+0600' or specdict['objid']==6):
        fitmask[np.logical_and(wave > 9715,wave < 9765)] = False
        
    #
    #####################
        
    # 0. Divide by a bspline spaced at every "binwidth" pixels to
    # smooth out any large scale variations caused by narrow emission
    # lines or dodgy sensitivity functions in pypeit reductions.
    # Rather than a median, this fits intentionally a bit high (upper
    # 75th %-ile) because that screens out broad absorption lines a bit better.
    # We will fix that systematic difference further below in iterative fit.

    if ('MOSFIRE' in instrument):
        binwidth  = 75
    else:
        binwidth  = 250
    halfwidth = binwidth // 2
    nbins     = (naxis1 // binwidth) - 1
    indx = halfwidth + np.arange(nbins) * binwidth
    gd = (fitmask[indx] == True)
    xbin = restwave[indx[gd]]

    # Note this intentionally fits high (75 %-ile)
    # to eliminate dipping into absorption lines
    ybin = np.array([np.quantile(flux[ii-halfwidth:ii+halfwidth],0.5) for ii in indx[gd]])

    # Tidy up some areas we won't use, and save points for later
    # ybin[xbin < 1216] = ybin[xbin > 1216][0]
    ybin[ybin < 1e-20] = np.median(ybin)
    xbin1 = xbin
    ybin1 = ybin

    # Fit the knots with a bspline and normalize it out
    t, c, k = interpolate.splrep(xbin, ybin, k=3, s=0)
    flux_filt = interpolate.BSpline(t,c,k)

    if (False):
        plt.plot(wave,flux)
        plt.plot(xbin*(1+redshift),ybin,'*',color='c')
        plt.plot(wave,flux_filt(restwave))
        plt.show()

    flux = flux / flux_filt(restwave)
    ivar = ivar * (flux_filt(restwave))**2
    fitmask[flux*np.sqrt(ivar) < 3] = False

    
    ##################################
    #
    # 1. Instatiate a mask for points to include in the fit by
    # refitting an interpolated function to the normalized spectrum
    # As before, bias the fit a bit high becuase we want to avoid dipping
    # into real absorption lines.

    if (instrument=='HIRES'):
        binwidth=250
    elif (instrument == 'FIRE_XSH'):
        binwidth = 15
    elif ('MOSFIRE' in instrument):
        binwidth = 50
    else:
        binwidth  = 60

    halfwidth = binwidth // 2
    nbins     = (naxis1 // binwidth) - 1
    indx = halfwidth + np.arange(nbins) * binwidth
    gd = (fitmask[indx] == True)

    xbin = np.array(restwave[indx[gd]])
    ybin = np.array([np.quantile(flux[ii-halfwidth:ii+halfwidth],0.75) for ii in indx[gd]])

    flux_filt2 = interpolate.interp1d(xbin, ybin, bounds_error=False, \
                                     fill_value='extrapolate')

    pk_mask, _ = find_peaks(flux,height=(flux_filt2(restwave),None))

    
    ######################################
    #
    # 2. Now we have pk_mask, which contains the points we will use.
    # Next is to iteratively fit and reject outliers using RANSAC
    
    # Set the number of iterations and rejection criteria, which
    # vary according to which spectrograph is used.
    
    if ('MOSFIRE' in instrument):
        niter=1
        discard = [4,3]
    elif (instrument == 'FIRE' or instrument == 'XSH_NIR'):
        niter=3
        discard = [4,3,2]
    elif (instrument == 'XSH_VIS'):
        niter=3
        discard = [4,3,2]
    else:
        niter=1
        discard = [10,5]

    cont = np.ones(len(restwave))
        
    for i in range(niter):
    
        if (i == 0):
            # f is the new continuum model, on first pass it is a fit
            # to the peaks found above
            f = interpolate.interp1d(restwave[pk_mask],flux[pk_mask],\
                                     bounds_error=False, \
                                     fill_value='extrapolate')
            inmask = fitmask
        else:
            # On 2nd and later passes, f is taken from the previous fit.
            f = c
            inmask = inlier_mask
            
        # 3. Apply RANSAC to detect outlying features (esp. absorption) and mask

        inmask_ind = np.where(inmask)[0]
        mad = np.median(np.abs(np.median(flux[inmask_ind]/cont[inmask_ind])-flux[inmask_ind]/cont[inmask_ind]))
        # print("MAD = {}".format(mad))
        
        ransac = linear_model.RANSACRegressor(random_state=0, \
                                              loss='absolute_loss',\
                                              residual_threshold=discard[i]*mad)

        fitwv = restwave[inmask_ind]
        fitflux = flux[inmask_ind]/cont[inmask_ind]
        fitivar = ivar[inmask_ind]/cont[inmask_ind]
        ransac.fit(fitwv.reshape(len(fitwv),1), fitflux, \
                   sample_weight=fitivar)

        inlier_mask  = fitmask
        
        # Recall that fitmask contains only the regions we masked by hand
        # We need to be sure to re-eliminate them at each iteration
        inlier_mask[fitmask == False] = False
        inlier_mask[np.logical_and(restwave > 1215, restwave < 1230)] = True
        
        # 4. Smooth the inlier data points to build the continuum

        gd = (fitmask[indx] == True)
        xbin = np.zeros(len(indx[gd]))
        ybin = np.zeros(len(indx[gd]))
        for i in range(len(indx[gd])):
            ii = (indx[gd])[i]
            xbin[i] = restwave[ii]
            medpoints = flux[ii-halfwidth:ii+halfwidth]
            medmask   = inlier_mask[ii-halfwidth:ii+halfwidth]
            ybin[i] = np.median(medpoints[medmask])
        # print(xbin,ybin)
            
        use = ~np.isnan(ybin)
        if (bspline):
            t, c, k = interpolate.splrep(xbin[use], ybin[use], k=3, s=0)
            c = interpolate.BSpline(t,c,k)
        else:
            c = interpolate.interp1d(restwave[inlier_mask],\
                                     medfilt(flux[inlier_mask],kernel_size=51), \
                                     fill_value='extrapolate')


        cont = c(restwave)

    # Go back and put back in the crude normalization that we did at the
    # start of the process, so that the output is the same as the input flux
    cont = cont * flux_filt(restwave)
    cont[restwave < 1216] = (cont[restwave > 1216])[0]
    flux = flux * flux_filt(restwave)

    fullcontin = c(all_restwave) * flux_filt(all_restwave)
    fullcontin[all_restwave < 1216] = (fullcontin[all_restwave > 1216])[0]
    
    if (True):
        plt.step(wave,flux,where='mid')
        plt.plot(wave,cont,color='r')
        plt.plot(all_wave,fullcontin,color='r')
        plt.plot(xbin[use]*(1+redshift),ybin[use]*flux_filt(xbin[use]), 'o',color='y',markersize=3)
        plt.ylim(-1.0*np.median(cont),3*np.median(cont))
        plt.show()

    # Write the output to disk as a fits binary table
    # in the same directory as the data file
    if (writefile):
        #infile   = loadQsoSpec(specdict['objid'],revision=specdict['revision'],files=True)
        #path     = os.path.dirname(infile[instrument])
        #contname = os.path.basename(infile[instrument])[:-5]+'_contin.fits'
        # outfile  = path+'/'+contname
        outfile  = writefile
        t = Table([all_wave,fullcontin],names=('wave','cont'))
        t.write(outfile, format='fits', overwrite=True)
        print(f"Writing: {outfile}")
        
    # return(xbin[use]*(1+redshift),ybin[use])


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

    # The VIS reduction of 1120 looks correupted somehow
    # xsh_vis_contin  = redsideContin(qso,'XSH_VIS',writefile=writefile)

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

def j1148():

    qso = loadQsoSpec(5)

    path='/Users/simcoe/Science/eiger/Cache/Keck/MOSFIRE/'

    output = path+'J1148_coadd_Y_1013_tellcorr_contin.fits'
    mosfire_y  = redsideContin(qso,'MOSFIRE_Y',writefile=output)
    output = path+'J1148_coadd_J_1013_tellcorr_contin.fits'
    mosfire_j  = redsideContin(qso,'MOSFIRE_J',writefile=output)
    output = path+'J1148_coadd_H_1011_tellcorr_contin.fits'
    mosfire_h  = redsideContin(qso,'MOSFIRE_H',writefile=output)
    output = path+'J1148_coadd_K_1011_tellcorr_contin.fits'
    mosfire_k  = redsideContin(qso,'MOSFIRE_K',writefile=output)
