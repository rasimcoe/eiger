import numpy
import matplotlib.pyplot as plt
from jwst import datamodels
from scipy.stats import binned_statistic_2d
import os
from astropy.visualization import simple_norm
import h5py
import numpy as np
from astropy.io import fits
import grismconf
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, Angle

import astropy.units as u


import numpy as np
from scipy import special
import warnings
import numpy
from scipy import interpolate
import astropy.constants


from lmfit import Model

def gaussian(x,totflux,c,x0,sigma):
    return totflux*((sigma)**-1 * (2*np.pi)**-0.5 *np.exp(-(x-x0)**2/(2*sigma**2)))+c   



def extract(ID,ra0,dec0,dataname,h5name,h5name_2): #typically h5name is ptsrcstamps, h5name_2 is galstamps
    with fits.open(dataname) as fin:
        h = fin[0].header
        
    # Find some information about the observing mode of this rate file.
    filt = h["FILTER"] # Filter name, e.g. F410M
    grism = h["PUPIL"][-1] # R or C
    module = h["MODULE"] # Which NIRCAM module, A or B
    print("Filter:",filt)
    print("grism:",grism)
    print("Module:",module)
 
    # We load a model of the background for this combination of Module, Filter and grism
    bck = fits.open("../simulate_mirage/MIRAGE_DATA/mirage_data/nircam/GRISM_NIRCAM/current/{}_mod{}_{}_back.norm.fits".format(filt,module,grism))[1].data

    # We initialize a datamodel WCS to the data
    grism_with_wcs = datamodels.open(dataname)

    world_to_pix = grism_with_wcs.meta.wcs.get_transform('world','detector')
    pix_to_world = grism_with_wcs.meta.wcs.get_transform('detector','world')

    # Compute the position of the source in the image in  pixel coordinates
    x0, y0, t, tt = world_to_pix(ra0,dec0,0,0)
    #x0, y0 = world_to_pix(ra0,dec0)

    print(dataname,"Target is at ",x0,y0)

    # We will need the gain file since data are not in e-/s
    if module=='A':
        gainfile = '/scratch/EIGER/simulate_mirage/MIRAGE_DATA/crds_cache/references/jwst/nircam/jwst_nircam_gain_0056.fits'
    if module=='B':
        gainfile = '/scratch/EIGER/simulate_mirage/MIRAGE_DATA/crds_cache/references/jwst/nircam/jwst_nircam_gain_0054.fits'
    gain = fits.open(gainfile)[1].data

    
    # We use our simulation to figure our where our spectrum is, approximately.
    # models will be an array  containinig ALL of the spectra from all of the sources
    # model0 is an image containing only a simulation of our target
    # Note that models-model0 will therefore be what we think contamination is.
    #print(np.shape(bck))
    models = np.zeros(np.shape(bck),np.float)
    #print(np.shape(models))
    stp = h5py.File(h5name,"r")
    minx0 = None
    for k in list(stp.keys()):
        minx = stp[k].attrs["minx"]
        maxx = stp[k].attrs["maxx"]
        miny = stp[k].attrs["miny"]
        maxy = stp[k].attrs["maxy"]
        m = stp[k][:] 
        models[miny:maxy+1,minx:maxx+1] += m 

        if k=="{}_+1".format(ID):
            model0 = m
            minx0 = stp[k].attrs["minx"]
            maxx0 = stp[k].attrs["maxx"]
            miny0 = stp[k].attrs["miny"]
            maxy0 = stp[k].attrs["maxy"]

    stp.close()

    stp = h5py.File(h5name_2,"r") #this is to also have galstamps

    for k in list(stp.keys()):
        minx = stp[k].attrs["minx"]
        maxx = stp[k].attrs["maxx"]
        miny = stp[k].attrs["miny"]
        maxy = stp[k].attrs["maxy"]
        m = stp[k][:] 
        models[miny:maxy+1,minx:maxx+1] += m 
        #print(k,np.nanmax(m))

        if k=="{}_+1".format(ID):
            model0 = m
            minx0 = stp[k].attrs["minx"]
            maxx0 = stp[k].attrs["maxx"]
            miny0 = stp[k].attrs["miny"]
            maxy0 = stp[k].attrs["maxy"]

    stp.close()
    if minx0==None:
        return None

    print('minx,maxx,miny,maxy',minx0,maxx0,miny0,maxy0)


    #MAKE EXTRACTION BOX BIGGER
    boxheight=maxy0-miny0
    delta_boxheight=50-boxheight
    oldmaxy0=maxy0
    oldminy0=miny0
    if boxheight < 50:
        maxy0=maxy0+int(delta_boxheight/2.)
        miny0=miny0-int(delta_boxheight/2.)
        newheight=1+maxy0-miny0
    else:
        newheight=boxheight

    modelold=model0
    model0=np.zeros((newheight,np.shape(model0)[1]))
    print(np.shape(model0))
    model0[int(delta_boxheight/2.):-int(delta_boxheight/2.),:]+=modelold
    model0[:int(delta_boxheight/2.),:]=0.
    #model0[-int(delta_boxheight/2.):,:]=0.




    # We load the observation, extimate and subtract the data, masking sources using the Mirage simulation
    data = fits.open(dataname)["SCI"].data * gain
    # Create a mask using of model of the whole thing. We mask out any objects with a flux greater  than 0.001 e-/s
    # We also mask part of the image which would be where we know the background is getting very small (10% level of maximum)
    #ok = (models<0.001) & (bck>0.1*np.max(bck)) 
    #ok = (bck>0.1*np.max(bck)) 
#
    ok = (models<0.001) & (bck>0.999)
    print('max bck',np.max(bck))

    # The sky level is now just the median of the ratio of our observation and model, when applying the mask we computed
    # to avoid any spectra
    sky_level = np.nanmedian((data/bck)[ok])
    back_use=np.zeros(np.shape(bck))
    back_use[ok]=1.
    fits.writeto('use_for_back.fits',back_use,overwrite=True)
    #fits.writeto('back.fits',bck,overwrite=True)

    print("Estimated background level:",sky_level)
    
    # We substract the bacvkground from our data
    data = data - sky_level*bck
    fits.writeto('background_subtracted.fits',data,overwrite=True)

    # We trim our data to be the stamp containing the spectrum we want to extract
    # Since we use the minx0,maxx0,miny0,maxy0 coordinates from our simulation, this stamp should
    # be just like our simulated data but containing the actual observation.
    data = data[miny0:maxy0+1,minx0:maxx0+1]




    # Same for error array
    err = fits.open(dataname)["ERR"].data*gain
    fits.writeto('err.fits',err,overwrite=True)

    err = err[miny0:maxy0+1,minx0:maxx0+1]

    # Same for DQ array
    dq = fits.open(dataname)["DQ"].data
    dq = dq[miny0:maxy0+1,minx0:maxx0+1]

    # We estimate of the contamination, from our Mirage simulation of the whole field, as noted
    # above, contamination is just models-model0
    contam = models[miny0:maxy0+1,minx0:maxx0+1] - model0



    # We load trhe grismconf module and load the proper configuration for this NIRCAM mode
    # Description of what grismconf does and how grism calibration  is done is available in 
    # WFC3 ISR 2017-01 from Pirzkal et al. (http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2017-01.pdf)
    C = grismconf.Config("../simulate_mirage/MIRAGE_DATA/mirage_data/nircam/GRISM_NIRCAM/current/NIRCAM_{}_mod{}_{}.conf".format(filt,module,grism))

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
        dys = C.DISPY("+1",x0,y0,ts) + ys  - 2*C.DISPY("+1",x0,y0,ts) ## I Dont know why  I need to do the -2*C.. JM 15 marh 2022
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
        m = model0
        l = ws
        d = data
        c = contam
        e = err
        q = dq
        y = dys

    # We now create optimal extraction weights, using our simulation (which accounts for the object profile etc)
    # The extraxction weights are based on the normalized simulated spectral profile at each wavelength
    # Sum up the model in the y-direction and replicate that 
    ysum = np.sum(m,axis=0) 
    w = np.repeat([ysum], np.shape(m)[0], axis=0)
    weight = m/w
    
    # We make sure our original data only has valid data and return everything we computed
    # What we return are no longer 2D stamps but 1D vectors
    ok = np.isfinite(d)

    print(len(d[ok]),'<<<<------')

    
    return  d[ok], e[ok], q[ok], l[ok], c[ok], m[ok], weight[ok], y[ok], C
    #return  d, e, q, l, c, m, weight, y, C


def extract_with_emline(ID,ra0,dec0,dataname,SBE_method,DIRECT_SBE,h5name,h5name_2,to_replace): #typically h5name is ptsrcstamps, h5name_2 is galstamps
    with fits.open(dataname) as fin:
        h = fin[0].header
        
    # Find some information about the observing mode of this rate file.
    filt = h["FILTER"] # Filter name, e.g. F410M
    grism = h["PUPIL"][-1] # R or C
    module = h["MODULE"] # Which NIRCAM module, A or B
    print("Filter:",filt)
    print("grism:",grism)
    print("Module:",module)
 
     # We load trhe grismconf module and load the proper configuration for this NIRCAM mode
    # Description of what grismconf does and how grism calibration  is done is available in 
    # WFC3 ISR 2017-01 from Pirzkal et al. (http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2017-01.pdf)
    C = grismconf.Config("../simulate_mirage/MIRAGE_DATA/mirage_data/nircam/GRISM_NIRCAM/current/NIRCAM_{}_mod{}_{}.conf".format(filt,module,grism))

    # We load a model of the background for this combination of Module, Filter and grism
    bck = fits.open("../simulate_mirage/MIRAGE_DATA/mirage_data/nircam/GRISM_NIRCAM/current/{}_mod{}_{}_back.norm.fits".format(filt,module,grism))[1].data


    #USE FLATFIELDSTEP:
    # We initialize a datamodel WCS to the data
    grism_with_wcs = datamodels.open(dataname)

    world_to_pix = grism_with_wcs.meta.wcs.get_transform('world','detector')
    pix_to_world = grism_with_wcs.meta.wcs.get_transform('detector','world')

    # Compute the position of the source in the image in  pixel coordinates
    x0, y0, t, tt = world_to_pix(ra0,dec0,0,0)
    #x0, y0 = world_to_pix(ra0,dec0)

    print(dataname,"Target is at (old) ",x0,y0)

    ###x0=x0+1.
    #USE RATE UHEAD_CAL
    # hdul=fits.open(dataname)
    # wcs=WCS(hdul['SCI'].header)
    # x0,y0=wcs.all_world2pix(ra0,dec0,0.) #lamb_line=

    # print('New:',x0,y0)



    # We will need the gain file since data are not in e-/s
    if module=='A':
        gainfile = '/scratch/EIGER/simulate_mirage/MIRAGE_DATA/crds_cache/references/jwst/nircam/jwst_nircam_gain_0056.fits'
    if module=='B':
        gainfile = '/scratch/EIGER/simulate_mirage/MIRAGE_DATA/crds_cache/references/jwst/nircam/jwst_nircam_gain_0054.fits'
    gain = fits.open(gainfile)[1].data
    if SBE_method==True:
    # We use our simulation to figure our where our spectrum is, approximately.
    # models will be an array  containinig ALL of the spectra from all of the sources
    # model0 is an image containing only a simulation of our target
    # Note that models-model0 will therefore be what we think contamination is.
    #print(np.shape(bck))
        models = np.zeros(np.shape(bck),np.float)

        models_id = np.zeros(np.shape(bck),np.float)

        #print(np.shape(models))
        stp = h5py.File(h5name,"r")

        minx0 = None
        ii=0
        for k in list(stp.keys()):
            #print(k)
            ii+=1
            minx = stp[k].attrs["minx"]
            maxx = stp[k].attrs["maxx"]
            miny = stp[k].attrs["miny"]
            maxy = stp[k].attrs["maxy"]
            m = stp[k][:] 
            models[miny:maxy+1,minx:maxx+1] += m 
            thism=models[miny:maxy+1,minx:maxx+1] 

            models_id[miny:maxy+1,minx:maxx+1][thism>0.05] +=ii#float(k.replace('_+1','')) #still wrong


            if k=="{}_+1".format(ID):
                model0 = m
                minx0 = stp[k].attrs["minx"]
                maxx0 = stp[k].attrs["maxx"]
                miny0 = stp[k].attrs["miny"]
                maxy0 = stp[k].attrs["maxy"]


        stp.close()

        if DIRECT_SBE==False:
            stp = h5py.File(h5name_2,"r") #this is to also have galstamps

            for k in list(stp.keys()):
                ii+=1
                minx = stp[k].attrs["minx"]
                maxx = stp[k].attrs["maxx"]
                miny = stp[k].attrs["miny"]
                maxy = stp[k].attrs["maxy"]
                m = stp[k][:] 
                models[miny:maxy+1,minx:maxx+1] += m 

                thism=models[miny:maxy+1,minx:maxx+1]
                models_id[miny:maxy+1,minx:maxx+1][thism>0.05] +=ii#float(k.replace('_+1',''))

                #print(k,np.nanmax(m))

                if k=="{}_+1".format(ID):
                    model0 = m
                    minx0 = stp[k].attrs["minx"]
                    maxx0 = stp[k].attrs["maxx"]
                    miny0 = stp[k].attrs["miny"]
                    maxy0 = stp[k].attrs["maxy"]

            stp.close()
        if minx0==None:
            return None
        if minx0<0:
            minx0=0.


    if SBE_method==False:
        ysize=10
        ts0=C.INVDISPX("+1",x0,y0,0)
        ts1=C.INVDISPX("+1",x0,y0,1)

        minx0,maxx0=int(C.DISPX("+1",x0,y0,0)+x0),int(C.DISPX("+1",x0,y0,1)+x0)
        miny0,maxy0=int(y0-ysize+C.DISPY("+1",x0,y0,0.5)),int(y0+ysize+C.DISPY("+1",x0,y0,0.5))

        print(minx0,maxx0,miny0,maxy0,'STOPHERE')
        
        if module=='B':
            oldmaxx0=np.copy(maxx0)
            oldminx0=np.copy(minx0)
            minx0=oldmaxx0
            maxx0=oldminx0

        print(maxy0,miny0)
        if miny0<0:
            miny0=0
        if maxy0>2048:
            maxy0=2047
        #if minx0>=2047:
           # just_crash


        if maxx0>2048:
            maxx0=2047
        if minx0<0:
            minx0=0

        models = np.zeros(np.shape(bck),np.float)
        models_id = np.zeros(np.shape(bck),np.float)
        model0=models[miny0:maxy0+1,minx0:maxx0+1]


    print('minx,maxx,miny,maxy',minx0,maxx0,miny0,maxy0)

    #MAKE EXTRACTION BOX BIGGER
    boxheight=maxy0-miny0
    delta_boxheight=50-boxheight
    oldmaxy0=maxy0
    oldminy0=miny0
    if boxheight < 50:
        maxy0=maxy0+int(delta_boxheight/2.)
        miny0=miny0-int(delta_boxheight/2.)
            
        newheight=1+maxy0-miny0
        modelold=model0

        ly=np.shape(modelold)

        model0=np.zeros((newheight,np.shape(model0)[1]))
        print(np.shape(model0))
        model0[int(delta_boxheight/2.):-int(delta_boxheight/2.),:]+=modelold
        model0[:int(delta_boxheight/2.),:]=0.
    else:
        newheight=boxheight

    #model0[-int(delta_boxheight/2.):,:]=0.

    print(np.shape(model0),'<<<<')


    # We load the observation, extimate and subtract the data, masking sources using the Mirage simulation
    data = fits.open(dataname)["SCI"].data * gain
    bck=bck*gain
    # Create a mask using of model of the whole thing. We mask out any objects with a flux greater  than 0.001 e-/s
    # We also mask part of the image which would be where we know the background is getting very small (10% level of maximum)
    #ok = (models<0.001) & (bck>0.1*np.max(bck)) 
    #ok = (bck>0.1*np.max(bck)) 
#
#    ok = (models<0.001) & (bck>0.999)
    ok = (models<0.001) & (bck>0.1*np.max(bck))

    print('max bck',np.max(bck))

    # The sky level is now just the median of the ratio of our observation and model, when applying the mask we computed
    # to avoid any spectra
    sky_level = np.nanmedian((data/bck)[ok])
    back_use=np.zeros(np.shape(bck))
    back_use[ok]=1.
    #fits.writeto('use_for_back.fits',back_use,overwrite=True)
    #fits.writeto('back.fits',bck,overwrite=True)

    print("Estimated background level:",sky_level)
    
    # We substract the bacvkground from our data
    data = data - sky_level*bck
    #fits.writeto('background_subtracted.fits',data,overwrite=True)

    # We trim our data to be the stamp containing the spectrum we want to extract
    # Since we use the minx0,maxx0,miny0,maxy0 coordinates from our simulation, this stamp should
    # be just like our simulated data but containing the actual observation.
    data = data[miny0:maxy0+1,minx0:maxx0+1]



    # Same for error array
    err = fits.open(dataname)["ERR"].data*gain
    #fits.writeto('err.fits',err,overwrite=True)

    err=err**2 #Change to variance

    err = err[miny0:maxy0+1,minx0:maxx0+1]

    # Same for DQ array
    dq = fits.open(dataname)["DQ"].data
    dq = dq[miny0:maxy0+1,minx0:maxx0+1]


    # Same for Emission-line map
    emname=dataname.replace(to_replace,'emline_oldhead')
    contname=dataname.replace(to_replace,'continua_oldhead')

    em = fits.open(emname)["SCI"].data * gain
    em = em[miny0:maxy0+1,minx0:maxx0+1]

    cont = fits.open(contname)["SCI"].data * gain
    cont = cont[miny0:maxy0+1,minx0:maxx0+1]

    # We estimate of the contamination, from our Mirage simulation of the whole field, as noted
    # above, contamination is just models-model0
    contam = models[miny0:maxy0+1,minx0:maxx0+1] - model0

    contam_id=models_id[miny0:maxy0+1,minx0:maxx0+1]


    # These are the coordinates of all the pixels in our 2D stamp, but in the full image (wrt calibration is known)
    ys,xs = np.indices((maxy0-miny0+1,maxx0-minx0+1))
    # xs and ys are now the relatice dx and dy offsets from the position of our source
    # They are both 2D arrays of x and y coordinates.
    xs = xs + minx0 - x0 
    ys = ys + miny0 - y0

    # Depending on whether the grism disperse in the x or y direction, we use the INVDISPX or INVDISPY functions
    # to compute the value for t for every pixel in our 2D stamps
    if grism=="R":
        print('-----',xs)
        ts = C.INVDISPX("+1",x0,y0,xs)
        dys = C.DISPY("+1",x0,y0,ts) + ys  - 2*C.DISPY("+1",x0,y0,ts) ## I Dont know why  I need to do the -2*C.. JM 15 marh 2022
    if grism=="C":
        ts = C.INVDISPY("+1",x0,y0,ys)
        dys = C.DISPX("+1",x0,y0,ts) + xs
    

    # Now compute the wavelength of every pixel in our 2D stamp
    ws = C.DISPL("+1",x0,y0,ts)


    print(np.shape(ws),np.shape(ts),np.shape(xs),'<<<')


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

        cont=np.transpose(cont) #continum map Daichi-method
        em=np.transpose(em) #emission-line map Daichi-method
        contam_id=np.transpose(contam_id)

    if grism=="R":
        m = model0
        l = ws
        d = data
        c = contam
        e = err
        q = dq
        y = dys


    # We now create optimal extraction weights, using our simulation (which accounts for the object profile etc)
    # The extraxction weights are based on the normalized simulated spectral profile at each wavelength
    # Sum up the model in the y-direction and replicate that 
    ysum = np.sum(m,axis=0) 
    w = np.repeat([ysum], np.shape(m)[0], axis=0)
    weight = m/w
    
    # We make sure our original data only has valid data and return everything we computed
    # What we return are no longer 2D stamps but 1D vectors
    ok = np.isfinite(d)
    print(np.shape(d),d)
    senscorrect=True
    if senscorrect==True:
        lam = np.nanmean(l,axis=0)
        C = grismconf.Config("../simulate_mirage/MIRAGE_DATA/mirage_data/nircam/GRISM_NIRCAM/current/NIRCAM_%s_mod%s_R.conf"%(filt,module))
        s = 1E-18 * C.SENS["+1"](lam) #ALways use the one for module A because Module B data has been rescaled.

        ffilter=np.zeros(np.shape(d)) 
        for jj in range(len(ffilter[:,0])):
            ffilter[jj,:]=s


        d=d/s
        e=e/s
        m=m/s
        c=c/s
        em=em/s
        cont=cont/s

    print(np.shape(d),d)
    print('^^^^')
  
    #return  d[ok], e[ok], q[ok], l[ok], c[ok], m[ok], weight[ok], y[ok], cont[ok], em[ok], contam_id[ok], C
    return  d, e, q, l, c, m, weight, y, cont, em, contam_id, C




def extract_emline_only(ID,ra0,dec0,dataname,to_replace): 
    with fits.open(dataname) as fin:
        h = fin[0].header
        
    # Find some information about the observing mode of this rate file.
    filt = h["FILTER"] # Filter name, e.g. F410M
    grism = h["PUPIL"][-1] # R or C
    module = h["MODULE"] # Which NIRCAM module, A or B
    #print("Filter:",filt)
    #print("grism:",grism)
    #print("Module:",module)
 
     # We load trhe grismconf module and load the proper configuration for this NIRCAM mode
    # Description of what grismconf does and how grism calibration  is done is available in 
    # WFC3 ISR 2017-01 from Pirzkal et al. (http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2017-01.pdf)
    C = grismconf.Config("../simulate_mirage/MIRAGE_DATA/mirage_data/nircam/GRISM_NIRCAM/current/NIRCAM_{}_mod{}_{}.conf".format(filt,module,grism))

    # We load a model of the background for this combination of Module, Filter and grism
    bck = fits.open("../simulate_mirage/MIRAGE_DATA/mirage_data/nircam/GRISM_NIRCAM/current/{}_mod{}_{}_back.norm.fits".format(filt,module,grism))[1].data

    # We initialize a datamodel WCS to the data
    grism_with_wcs = datamodels.open(dataname)

    world_to_pix = grism_with_wcs.meta.wcs.get_transform('world','detector')
    pix_to_world = grism_with_wcs.meta.wcs.get_transform('detector','world')

    # Compute the position of the source in the image in  pixel coordinates
    x0, y0, t, tt = world_to_pix(ra0,dec0,0,0)

    print(dataname,"Target is at ",x0,y0)

    # We will need the gain file since data are not in e-/s
    if module=='A':
        gainfile = '/scratch/EIGER/simulate_mirage/MIRAGE_DATA/crds_cache/references/jwst/nircam/jwst_nircam_gain_0056.fits'
    if module=='B':
        gainfile = '/scratch/EIGER/simulate_mirage/MIRAGE_DATA/crds_cache/references/jwst/nircam/jwst_nircam_gain_0054.fits'
    gain = fits.open(gainfile)[1].data

    gain=np.zeros(np.shape(gain))+1. ######ZZZZ CHECK

    ysize=20
    minx0,maxx0=int(C.DISPX("+1",x0,y0,0)+x0),int(C.DISPX("+1",x0,y0,1)+x0)
    miny0,maxy0=int(y0-ysize+C.DISPY("+1",x0,y0,0.5)),int(y0+ysize+C.DISPY("+1",x0,y0,0.5))

    if module=='B':
        oldmaxx0=np.copy(maxx0)
        oldminx0=np.copy(minx0)
        minx0=oldmaxx0
        maxx0=oldminx0

    #print(maxy0,miny0)
    if miny0<0:
        miny0=0
    if maxy0>2048:
        maxy0=2047

    if maxx0>2048:
        maxx0=2047
    if minx0<0:
        minx0=0

    models = np.zeros(np.shape(bck),np.float)
    models_id = np.zeros(np.shape(bck),np.float)
    model0=models[miny0:maxy0+1,minx0:maxx0+1]


    print('minx,maxx,miny,maxy',minx0,maxx0,miny0,maxy0)


    #MAKE EXTRACTION BOX BIGGER
    boxheight=maxy0-miny0
    delta_boxheight=50-boxheight
    oldmaxy0=maxy0
    oldminy0=miny0
    if boxheight < 50:
        maxy0=maxy0+int(delta_boxheight/2.)
        miny0=miny0-int(delta_boxheight/2.)
            
        newheight=1+maxy0-miny0
        modelold=model0

        ly=np.shape(modelold)

        model0=np.zeros((newheight,np.shape(model0)[1]))
        #print(np.shape(model0))
        model0[int(delta_boxheight/2.):-int(delta_boxheight/2.),:]+=modelold
        model0[:int(delta_boxheight/2.),:]=0.
    else:
        newheight=boxheight


    # We load the observation, extimate and subtract the data, masking sources using the Mirage simulation
    data = fits.open(dataname)["SCI"].data * gain
    bck=bck*gain
    # We also mask part of the image which would be where we know the background is getting very small (10% level of maximum)
    ok = (bck>0.1*np.max(bck))

    sky_level = np.nanmedian((data/bck)[ok])
    back_use=np.zeros(np.shape(bck))
    back_use[ok]=1.

    print("Estimated background level:",sky_level)
    # We substract the bacvkground from our data
    data = data - sky_level*bck

    # We trim our data to be the stamp containing the spectrum we want to extract
    # Since we use the minx0,maxx0,miny0,maxy0 coordinates from our simulation, this stamp should
    # be just like our simulated data but containing the actual observation.
    data = data[miny0:maxy0+1,minx0:maxx0+1]




    # Same for error array
    err = fits.open(dataname)["ERR"].data*gain
    #fits.writeto('err.fits',err,overwrite=True)

    err=err**2 #Change to variance

    err = err[miny0:maxy0+1,minx0:maxx0+1]

    # Same for DQ array
    dq = fits.open(dataname)["DQ"].data
    dq = dq[miny0:maxy0+1,minx0:maxx0+1]

    # Same for Emission-line map
    emname=dataname.replace(to_replace,'emline_oldhead')

    em = fits.open(emname)["SCI"].data * gain
    em = em[miny0:maxy0+1,minx0:maxx0+1] #* gain


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
        dys = C.DISPY("+1",x0,y0,ts) + ys  - 2*C.DISPY("+1",x0,y0,ts) ## I Dont know why  I need to do the -2*C.. JM 15 marh 2022
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
        e = np.transpose(err)    # THe data error estimates in each pixel
        y = np.transpose(dys)    # The cross-dispersion distance of each pixel from the trace

        em=np.transpose(em) #emission-line map Daichi-method
        q=np.transpose(dq)
    if grism=="R":
        m = model0
        l = ws
        d = data
        e = err
        y = dys
        q = dq

    # We make sure our original data only has valid data and return everything we computed
    # What we return are no longer 2D stamps but 1D vectors
   # fits.writeto('test2d_d.fits',d,overwrite=True)
   # fits.writeto('test2d_y.fits',y,overwrite=True)
   # fits.writeto('test2d_l.fits',l,overwrite=True)


    ok = np.isfinite(d)
    
    return  d[ok], e[ok], q[ok], l[ok], y[ok], em[ok], C



def create_cutout(directimage,RA,DEC,lx,ly):
    #directimage='reduced/test_10bright_mag26_z4to7/imaging_F356W/step_i2d.fits'
    hdu = fits.open(directimage)
    wcs = WCS(hdu['SCI'].header)
    positions=SkyCoord(RA,DEC,unit="deg")
    size = (ly,lx)
    cutout = Cutout2D(hdu['SCI'].data, position=positions, size=size, wcs=wcs)
    hdu['SCI'].data = cutout.data
    hdu['SCI'].header.update(cutout.wcs.to_header())

    return cutout.data,hdu['SCI'].header


def xy_drizzled_to_xy_rate(emfile,drizzlefile,x_drizzled,y_drizzled):
    #the emfile is the single emission-line fine (from a _cal.fits) corresponding to a rate.fits 
    #the outout are the pixel positions in the original rate file
    hdul = fits.open(emfile)
    wcs_em=WCS(hdul['SCI'].header) 

    hdu2=fits.open(drizzlefile)
    wcs_drizzle=WCS(hdu2['SCI'].header)

    ra_pix,dec_pix=wcs_drizzle.all_pix2world(x_drizzled,y_drizzled,0)
    x_rate,y_rate=wcs_em.all_world2pix(ra_pix,dec_pix,0)

    return x_rate,y_rate
    #The rate file needs to have a WCS assigned.
    #If you want to do this:
    #from jwst.assign_wcs.assign_wcs_step import AssignWcsStep          
    #step = AssignWcsStep()
    #step.save_results = False#True
    #step.output_dir = output_dir
    #result = step.run(rate_file)



def xylamb_dispersed_to_radec(ratefile,x_line,y_line,obs_lamb_line,configfile='../simulate_mirage/MIRAGE_DATA/mirage_data/nircam/GRISM_NIRCAM/current/NIRCAM_F356W_modA_R.conf'):
    #delta_lamb_int =0.0144 O3_5008 to Hbeta
    #ratefile='jw01243001001_01101_00001_nrca5_flatfieldstep.fits'
    C = grismconf.Config(configfile)
    t= C.INVDISPL('+1',x_line,y_line,obs_lamb_line)
    dx=C.DISPX('+1',x_line,y_line,t) 
    dy=C.DISPY('+1',x_line,y_line,t) 
    x0=x_line-dx
    y0=y_line-dy

    hdul=fits.open(ratefile)
    wcs=WCS(hdul['SCI'].header)
    RA,DEC=wcs.all_pix2world(x0,y0,0)
    #RA,DEC,a,b=pix_to_world(x0,y0,0,0)
    return RA,DEC



def radeclambda_to_dispersed(emline_primdit,RA,DEC,obs_lamb_line,configfile='../simulate_mirage/MIRAGE_DATA/mirage_data/nircam/GRISM_NIRCAM/current/NIRCAM_F356W_modA_R.conf'):
    #delta_lamb_int =0.0144 O3_5008 to Hbeta
    #ratefile='jw01243001001_01101_00001_nrca5_flatfieldstep.fits'
    hdul=fits.open(emline_primdit)
    wcs=WCS(hdul['SCI'].header)
    x0,y0=wcs.all_world2pix(RA,DEC,0)

    C = grismconf.Config(configfile)
    print(x0,y0,obs_lamb_line)

    tlist=np.zeros(len(x0))
    for q in range(len(x0)):
        tlist[q]=C.INVDISPL('+1',x0[q],y0[q],obs_lamb_line[q])
   # t= C.INVDISPL('+1',x0[0],y0[0],obs_lamb_line[0])
    #print('::::',t,C.INVDISPL('+1',0,0,obs_lamb_line))
    dx=C.DISPX('+1',x0,y0,tlist) 
    dy=C.DISPY('+1',x0,y0,tlist) 
    x_grism=x0+dx
    y_grism=y0+dy


    #RA,DEC,a,b=pix_to_world(x0,y0,0,0)
    return x_grism,y_grism


def shift_columns(y,d,var,lamb,q,ysize=31,yoffset=-2):
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
        yaxis_new = np.arange(ysize)+0.5-ysize/2.  + yoffset  #y_trace_xaxis[ix]-15+np.arange(31)

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


def scrunch_columns(lambda_array,l,shifted_y,shifted_y_var,module):
    ly,lx=np.shape(shifted_y)
    len_lamb=len(lambda_array)
    scrunchd_x=np.zeros((ly,len_lamb))
    scrunchd_x_var=np.zeros((ly,len_lamb))
    scrunchd_lamb=np.zeros((ly,len_lamb))

    dwav=lambda_array[1]-lambda_array[0]



    for iy in range(ly):
        lambaxis=l[iy,:]
        flx_col = shifted_y[iy,:]
        var_col = shifted_y_var[iy,:]


        if module=='b':
            lambaxis=lambaxis[::-1]
            flx_col=flx_col[::-1]
            var_col=var_col[::-1]

        bad_col = var_col*0.
        bad_col[np.where((var_col<=0.0)|(np.isnan(var_col))|(np.isnan(flx_col)))[0]]=1.0

        ### Cumulative sum
       # flx_row_cum = np.append([0],np.nancumsum(flx_col)) #Update 22 aug, [0] append
        #bad_row_cum = np.append([0],np.nancumsum(bad_col)) #Update 22 aug, [0] append
        #var_row_cum = np.append([0],np.nancumsum(var_col)) #Update 22 aug, [0] append

        flx_row_cum = np.nancumsum(flx_col) #Update 22 aug, [0] append
        bad_row_cum = np.nancumsum(bad_col) #Update 22 aug, [0] append
        var_row_cum = np.nancumsum(var_col) #Update 22 aug, [0] append

  
        #print(np.flip(lambaxis)+0.5)
        #print((lambda_array-dwav/2.),np.flip(lambaxis)+0.5)
        flx_row_new0 = np.interp(lambda_array-dwav/2., lambaxis, flx_row_cum, left=np.nan, right=np.nan)
        flx_row_new1 = np.interp(lambda_array+dwav/2., lambaxis, flx_row_cum, left=np.nan, right=np.nan)
        var_row_new0 = np.interp(lambda_array-dwav/2., lambaxis, var_row_cum, left=np.nan, right=np.nan)
        var_row_new1 = np.interp(lambda_array+dwav/2., lambaxis, var_row_cum, left=np.nan, right=np.nan)
        bad_row_new0 = np.interp(lambda_array-dwav/2., lambaxis, bad_row_cum, left=np.nan, right=np.nan)
        bad_row_new1 = np.interp(lambda_array+dwav/2., lambaxis, bad_row_cum, left=np.nan, right=np.nan)

        bad_row_new = bad_row_new1-bad_row_new0
        flx_row_new = bad_row_new * 0.
        var_row_new = bad_row_new * 0.

        idx_bad  = np.where(bad_row_new>=0.5)[0]
        idx_good = np.where(bad_row_new<0.5)[0]

        bad_row_new[idx_bad]=1.0

        flx_row_new[idx_good] = (flx_row_new1-flx_row_new0)[idx_good]/(1-bad_row_new[idx_good])
        flx_row_new[idx_bad] = np.nan

        var_row_new[idx_good] = (var_row_new1-var_row_new0)[idx_good]/(1-bad_row_new[idx_good])
        var_row_new[idx_bad] = np.nan


        scrunchd_x[iy,:] = flx_row_new/(dwav*1E4) #this is bc of the /Angstrom term
        scrunchd_x_var[iy,:] = var_row_new/(dwav*1E4) #this is bc of the /Angstrom term


        scrunchd_lamb[iy,:]=lambda_array
    return scrunchd_x,scrunchd_x_var,scrunchd_lamb



    
def stack_with_reject_outliers(data,err,nobs, m = 5., method='mean',iterations=1):
    for i in range(iterations):
    	median=np.nanmedian(data,axis=0)
    	std=np.nanstd(data,axis=0)
    	if i>1:
    		mask_outliers+=np.abs((data-median))/std > m
    	else:
    		mask_outliers=np.abs((data-median))/std > m
    	data[mask_outliers]=np.nan
    	
    var=err**2 
    var[mask_outliers]=np.nan

    if method=='mean':
        fdata = np.nanmean(data,axis=0)
    elif method=='median':
        fdata = np.nanmedian(data,axis=0)
    
    return fdata,(np.nanmean(var,axis=0))**0.5, mask_outliers

    
def stack_with_mask_outliers(data, mask_outliers, method='mean'):
    data[mask_outliers] = np.nan
    if method=='mean':
        fdata = np.nanmean(data,axis=0)
    elif method=='median':
        fdata = np.nanmedian(data,axis=0)
    
    return fdata


def create_simple_optweight(wavelength,scrunched_data,kernel):
    yrange=np.array(range(len(scrunched_data[:,0])))
    print(len(scrunched_data[:,0]))
    opt_weight=np.zeros(np.shape(scrunched_data))

    midy=0.5+len(scrunched_data[:,0])/2.
    model = Model(gaussian) #

    params = model.make_params(totflux=1000.,c=0.,x0=midy,sigma=kernel) 
    opti=model.eval(params,x=yrange)
    opti=opti/np.nansum(opti)


    for bb in range(len(opt_weight[0,:])):
        opt_weight[:,bb]=opti
    return opt_weight




def extract_dk(ra0, dec0,
               scifile,
               emfile,
               contfile,
               ratefile, 
               yoffset=0., 
               yhsize=20, 
               use_wcs=False, 
               senscorrect=False):
    
    hdul = fits.open(scifile) 
    wcs= WCS(hdul[1].header)

    grism_with_wcs = datamodels.open(ratefile)

    ### contfile=emfile.replace('emline.fits','continua.fits')
    hdul_em=fits.open(emfile)
    hdul_cont=fits.open(contfile)

    # Get header
    h = hdul[0].header
    
    # Find some information about the observing mode of this rate file.
    filt = h["FILTER"] # Filter name, e.g. F410M
    grism = h["PUPIL"][-1] # R or C
    module = h["MODULE"] # Which NIRCAM module, A or B
    print("Filter:",filt)
    print("grism:",grism)
    print("Module:",module)
    
    condir='/scratch/EIGER/simulate_mirage/MIRAGE_DATA/mirage_data/nircam/GRISM_NIRCAM/current/'
    configfile = os.path.join(condir, 'NIRCAM_'+filt+'_mod'+module.upper()+'_'+grism+'.conf')
    C = grismconf.Config(configfile)

    # We initialize a datamodel WCS to the data
    # grism_with_wcs = datamodels.open(dataname)
    world_to_pix = grism_with_wcs.meta.wcs.get_transform('world','detector')
    #pix_to_world = grism_with_wcs.meta.wcs.get_transform('detector','world')
    x0, y0,empty,empty2 = world_to_pix(ra0,dec0,0,0)
    if use_wcs:
        x0, y0 = wcs.world_to_pixel(SkyCoord(ra=ra0, dec=dec0, unit='degree'))
        
        
    # Compute the position of the source in the image in  pixel coordinates
    
    # Offset correction if given
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
    
    print('minx0,maxx0,miny0,maxy0',minx0,maxx0,miny0,maxy0)
    if minx0>2047 or maxx0 < 0 or miny0>2047 or maxy0 < 0:
        print('ERROR: Spectrum is out of fiels.')
        print('--- minx0, maxx0: ', minx0, maxx0)
        print('--- miny0, maxy0: ', miny0, maxy0)
        return None
        
    # We load the observation, extimate and subtract the data, masking sources using the Mirage simulation
    data = hdul["SCI"].data  ## MJy/sr
    #ratefile=ratefile.replace('flatfieldstep','rate')
    #ratedata=fits.getdata(ratefile)

    emdata = hdul_em["SCI"].data  ## MJy/sr
    contdata = hdul_cont["SCI"].data  ## MJy/sr

    err = hdul["ERR"].data 
    dq = hdul["DQ"].data
  
    # We trim our data to be the stamp containing the spectrum we want to extract
    # Since we use the minx0,maxx0,miny0,maxy0 coordinates from our simulation, this stamp should
    # be just like our simulated data but containing the actual observation.
    data = data[miny0:maxy0+1,minx0:maxx0+1]
    err = err[miny0:maxy0+1,minx0:maxx0+1]
    dq = dq[miny0:maxy0+1,minx0:maxx0+1]
    
    emdata = emdata[miny0:maxy0+1,minx0:maxx0+1]
    contdata = contdata[miny0:maxy0+1,minx0:maxx0+1]

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

        em = np.transpose(emdata)    # The cross-dispersion distance of each pixel from the trace
        con = np.transpose(contdata)    # The cross-dispersion distance of each pixel from the trace

    if grism=="R":
#        m = model0
        l = ws
        d = data
#        c = contam
        e = err
        q = dq
        y = dys
        em = emdata
        cont= contdata

    # We now create optimal extraction weights, using our simulation (which accounts for the object profile etc)
    # The extraxction weights are based on the normalized simulated spectral profile at each wavelength
    # Sum up the model in the y-direction and replicate that 
 #   ysum = np.sum(m,axis=0) 
 #   w = np.repeat([ysum], np.shape(m)[0], axis=0)
 #   weight = m/w
    
    # We make sure our original data only has valid data and return everything we computed
    # What we return are no longer 2D stamps but 1D vectors
    #ok = np.isfinite(d)

    ## SENSITIVITY CORRECTION WILL BE EXPLICITLY APPLIED AFTER "EXTRACT" 
    if senscorrect==True:
        lam = np.nanmean(l,axis=0)
        s = 1E-18 * C.SENS["+1"](lam) #ALways use the one for module A because Module B data has been rescaled.

        ffilter=np.zeros(np.shape(d)) 
        for jj in range(len(ffilter[:,0])):
            ffilter[jj,:]=s

        d=d/s
        e=e/s
        em=em/s
        cont=cont/s
    
#    return  d[ok], e[ok], q[ok], l[ok], c[ok], m[ok], weight[ok], y[ok], C
    #return  d[ok], e[ok], q[ok], l[ok], 0,      0,     0, y[ok], C
    return  d, e, q, l, 0,      0,     0, y, em, cont, C#, filt, grism, module




def shift_columns_dk(y,d,var,lamb,q,ysize=31):
    
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


def scrunch_columns_dk(lambda_array,l,shifted_y,shifted_y_var,module,flambda=False):

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
        ##-dwav in lambrgrid-dwav added JM 30 Nov 22 ==> DK deleted 
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





def extract_fresco(ra0, dec0,scifile,emfile,ratefile, yoffset=0., yhsize=20, use_wcs=False):
    hdul = fits.open(scifile)
    wcs= WCS(hdul[1].header)


    grism_with_wcs = datamodels.open(ratefile)

    contfile=emfile.replace('emline.fits','continua.fits')
    hdul_em=fits.open(emfile)
    hdul_cont=fits.open(contfile)

    # Get header
    h = hdul_em[0].header
    
    # Find some information about the observing mode of this rate file.
    filt = h["FILTER"] # Filter name, e.g. F410M
    grism = h["PUPIL"][-1] # R or C
    module = h["MODULE"] # Which NIRCAM module, A or B
    print("Filter:",filt)
    print("grism:",grism)
    print("Module:",module)
    
    condir='/scratch/EIGER/simulate_mirage/MIRAGE_DATA/mirage_data/nircam/GRISM_NIRCAM/current/'
    configfile = os.path.join(condir, 'NIRCAM_'+filt+'_mod'+module.upper()+'_'+grism+'.conf')
    C = grismconf.Config(configfile)

    # We initialize a datamodel WCS to the data
    # grism_with_wcs = datamodels.open(dataname)
    world_to_pix = grism_with_wcs.meta.wcs.get_transform('world','detector')
    #pix_to_world = grism_with_wcs.meta.wcs.get_transform('detector','world')
    x0, y0,empty,empty2 = world_to_pix(ra0,dec0,0,0)
    if use_wcs:
        x0, y0 = wcs.world_to_pixel(SkyCoord(ra=ra0, dec=dec0, unit='degree'))
        
        
    # Compute the position of the source in the image in  pixel coordinates
    
    y0 = y0 + yoffset
    #print(dataname,"Target is at ",x0,y0)

    xref=x0 #2048.
    yref=y0 #2048.

    t= C.INVDISPL('+1',xref,yref,np.array([3.79, 5.15]))
    dx=C.DISPX('+1',xref,yref,t) 
    dy=C.DISPY('+1',xref,yref,t) 
    x_line = x0 + dx #CHANGED
    y_line = y0 + dy
    minx0 = np.max([0, np.int32(np.min(x_line))])
    maxx0 = np.min([2047, np.int32(np.max(x_line))])
    miny0 = np.max([0, np.int32(np.min(y_line-yhsize))])
    maxy0 = np.min([2047, np.int32(np.max(y_line+yhsize))])
    
    print('minx0,maxx0,miny0,maxy0',minx0,maxx0,miny0,maxy0)
    if minx0>2047 or maxx0 < 0 or miny0>2047 or maxy0 < 0:
        print('ERROR: Spectrum is out of fiels.')
        print('--- minx0, maxx0: ', minx0, maxx0)
        print('--- miny0, maxy0: ', miny0, maxy0)
        return None
        
    # We load the observation, extimate and subtract the data, masking sources using the Mirage simulation
    data = hdul["SCI"].data  ## MJy/sr
    emdata = hdul_em["SCI"].data  ## MJy/sr
    contdata = hdul_cont["SCI"].data  ## MJy/sr



    err = hdul_em["ERR"].data 
    dq = hdul_em["DQ"].data
  
    # We trim our data to be the stamp containing the spectrum we want to extract
    # Since we use the minx0,maxx0,miny0,maxy0 coordinates from our simulation, this stamp should
    # be just like our simulated data but containing the actual observation.
    data = data[miny0:maxy0+1,minx0:maxx0+1]
    err = err[miny0:maxy0+1,minx0:maxx0+1]
    dq = dq[miny0:maxy0+1,minx0:maxx0+1]
    
    emdata = emdata[miny0:maxy0+1,minx0:maxx0+1]
    contdata = contdata[miny0:maxy0+1,minx0:maxx0+1]


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

        em = np.transpose(emdata)    # The cross-dispersion distance of each pixel from the trace
        con = np.transpose(contdata)    # The cross-dispersion distance of each pixel from the trace

    if grism=="R":
#        m = model0
        l = ws
        d = data
#        c = contam
        e = err
        q = dq
        y = dys
        em = emdata
        cont= contdata

    # We now create optimal extraction weights, using our simulation (which accounts for the object profile etc)
    # The extraxction weights are based on the normalized simulated spectral profile at each wavelength
    # Sum up the model in the y-direction and replicate that 
 #   ysum = np.sum(m,axis=0) 
 #   w = np.repeat([ysum], np.shape(m)[0], axis=0)
 #   weight = m/w
    
    # We make sure our original data only has valid data and return everything we computed
    # What we return are no longer 2D stamps but 1D vectors
    #ok = np.isfinite(d)
    senscorrect=True
    if senscorrect==True:
        lam = np.nanmean(l,axis=0)
        s = 1E-18 * C.SENS["+1"](lam) #ALways use the one for module A because Module B data has been rescaled.

        ffilter=np.zeros(np.shape(d)) 
        for jj in range(len(ffilter[:,0])):
            ffilter[jj,:]=s

        d=d/s
        e=e/s
        em=em/s
        cont=cont/s
    
#    return  d[ok], e[ok], q[ok], l[ok], c[ok], m[ok], weight[ok], y[ok], C
    #return  d[ok], e[ok], q[ok], l[ok], 0,      0,     0, y[ok], C
    return  d, e, q, l, 0,      0,     0, y, em, cont, C#, filt, grism, module



#### ======================================================================
def get_yoffset(x, y, module, yoffset_map=None):

    if yoffset_map is None:
        yoffset_map = np.load(eiger_reference_directory,
                              'yoffset_polyreg_F356W.R.Mod'+module.upper()+'_jw01076101.npy', 
                              allow_pickle=True)[()]
    else:
        try:
            yoffset_map = yoffset_map[module]
        except:
            pass

    if (x >= yoffset_map['xind'][0] and 
        x <= yoffset_map['xind'][-1] and
        y >= yoffset_map['yind'][0] and 
        y <= yoffset_map['yind'][-1]):
        x0 = int(x) - yoffset_map['xind'][0]
        y0 = int(y) - yoffset_map['yind'][0]
        yoffset = yoffset_map['yoffset_polyreg'][y0,x0]
        return yoffset
    else:
        return np.nan




def offset_wcs_in_sky(ra, dec, p):

    n = np.array(ra).size
    
    skycen = SkyCoord(ra=p[2], dec=p[3], unit='deg')
    sky = SkyCoord(ra=ra, dec=dec, unit='deg')
    sky_off = sky.spherical_offsets_by(np.repeat(p[0],n)*u.mas, 
                                       np.repeat(p[1],n)*u.mas)
    
    sky_off_t = skycen.spherical_offsets_to(sky_off)               

    c, s = np.cos(p[4]), np.sin(p[4])  ## p[4] in radian
    
    sky_new_t = SkyCoord(ra = (c*sky_off_t[0].deg - s*sky_off_t[1].deg),
                         dec = (s*sky_off_t[0].deg + c*sky_off_t[1].deg),
                         unit='deg')
    
    sky_new = skycen.spherical_offsets_by(sky_new_t.ra.deg*u.deg, 
                                          sky_new_t.dec.deg*u.deg)
    return sky_new.ra.deg, sky_new.dec.deg

def offset_wcs(x, y, p):
    x_off = x + p[0]
    y_off = y + p[1]
    
    x_off_t = x_off - p[2]
    y_off_t = y_off - p[3]
    
    c, s = np.cos(p[4]), np.sin(p[4])
    #R = np.array(((c,-s),(s,c)))
    
    x_new = (c*x_off_t - s*y_off_t) + p[2]
    y_new = (s*x_off_t + c*y_off_t) + p[3]
    return x_new, y_new


   


def offset_wcs_residual(pars, x, y, x2, y2):
    parvals = pars.valuesdict()
    p = [parvals['dx'],
         parvals['dy'],
         parvals['rot_xcen'],
         parvals['rot_ycen'],
         parvals['theta']]
    x_new, y_new = offset_wcs(x, y, p)
    residual = np.sqrt((x_new - x2)**2 + (y_new - y2)**2)
    return residual


def get_wcs_offset_params_dict(field, mode='sky'):

    if field=='J0100':
        #     fil = '/scratch/kashinod/EIGER/J0100/reduction_imaging/checkWCS_v2/wcs_offset_params_visit'+str(visit)+module.lower()+'.txt'
        #     params = np.loadtxt(fil)
        if mode=='sky':
            params_dict = {'1a':[-1.92284861e+00,  1.88044700e-02,  1.50258035e+01,  2.81031424e+01, -5.27552483e-04],
                           '1b':[-1.33417926e+00,  4.50286320e+00,  1.50127890e+01,  2.81313436e+01, -3.94051962e-04],
                           '2a':[-1.17131648e+02,  1.21901199e+02,  1.49362693e+01,  2.80544577e+01, -3.97039191e-04],
                           '2b':[-8.72268415e-01, -1.26062571e+01,  1.50350544e+01,  2.81392476e+01, -4.16844530e-04],
                           '3a':[-1.39928713e-03, -5.13176986e+00,  1.50692772e+01,  2.80485850e+01, -1.72439247e-04],
                           '3b':[-5.54050193e+00,  1.03332937e+00,  1.50516414e+01,  2.80459069e+01, -2.70985724e-04],
                           '4a':[-1.68919556e+00,  1.50266899e+01,  1.50257861e+01,  2.80349982e+01, -1.18307824e-04],
                           '4b':[-1.27300158e+01,  2.84521481e+01,  1.50054374e+01,  2.80265964e+01, -2.48477998e-04]}
        elif mode=='image':
            params_dict = {'1a':[ -0.679936091154,  0.523226446148, -1118.4604275112,  608.4828378726, 0.000558578932507],
                           '1b':[ -1.111210284948,  3.129704853141,   928.2162884476, 1283.2376923699, 0.000423325155808],
                           '2a':[ -1.084928574764,  1.745666945788,   750.9150277642, 1696.6610390372, 0.000428338992284],
                           '2b':[ -0.937111596682,  3.339523725388,   918.3636725379, 1240.0037692349, 0.000457844502586],
                           '3a':[  0.057536995593, -0.442108779774,   819.2640176016,  835.9240175748, 0.000182900260607],
                           '3b':[ -0.292402520885,  0.607881819359,  2168.5716962952, 2178.7603750517, 0.000290237599955],
                           '4a':[ -0.142066614727, -0.457647901154,   140.7245984545, 1392.8826646429, 0.000109932123407],
                           '4b':[  0.001534739317,  0.659164131917,  2793.8799930654,  881.4838458998, 0.000250044843380]}
    elif field=='J1148':
        # /net/galaxy-data/export/galaxydata/kashinod/EIGER/J1148/reduction_imaging/checkWCS_v2
        if mode=='sky':
            params_dict = {'1a':[-5.22052172e+00, -3.77166435e+01,  1.77116287e+02,  5.31836441e+01, -8.61888669e-05],
                           '1b':[-1.00233481e+00, -2.48378554e+01,  1.77032728e+02,  5.29710966e+01, -1.76271525e-04],
                           '2a':[-5.52135046e+00,  1.54133311e-01,  1.77075403e+02,  5.30631157e+01, -1.31854387e-04],
                           '2b':[-8.91541757e+01, -8.00305305e+00,  1.77086762e+02,  5.28954786e+01, -6.77776314e-06],
                           '3a':[ 2.96699564e+00, -1.04888479e+01,  1.77059016e+02,  5.28348837e+01, -8.91361290e-05],
                           '3b':[ 6.19708825e-01, -3.68786148e-01,  1.77045300e+02,  5.28330907e+01, -1.20132357e-04],
                           '4a':[ 7.82339319e+00, -2.14809686e-01,  1.77017638e+02,  5.28738512e+01, -1.44369631e-04],
                           '4b':[ 2.61422477e+00, -2.06289081e+00,  1.77030841e+02,  5.28486944e+01, -1.42136092e-04]}
            
        elif mode=='image':
            params_dict = {'1a': [0.837805413135, 2.487559346508, 3108.2124012270, 1267.1472409091, -0.000006582623919],
                           '1b': [0.787059632147, 2.245594900227, 1705.2217689061, -487.8656071018, -0.000029985894132],
                           '2a': [0.512913430975, 2.461776542105, 1064.2416491705,  829.4416137441, -0.000008740332443],
                           '2b': [0.467233841778, 2.321461559777, 1407.5529270521,   18.8176487703, -0.000061814449389],
                           '3a': [0.231866383239,-1.050091020704,-1691.4364154766,  596.7843987598,  0.000135110674380],
                           '3b': [0.505132559356,-0.306546753254,  966.1245768030,-1079.5363898621,  0.000193306148084],
                           '4a': [0.255913643120,-0.636130396339, 1272.7314351786,  392.8503631849,  0.000185581786706],
                           '4b': [0.776445325796,-0.339055684816,  145.3952998526,-2566.6119666198,  0.000186876742781]}
    elif field=='J1120':
        # /net/galaxy-data/export/galaxydata/kashinod/EIGER/J1120/reduction_imaging/checkWCS_v2
        if mode=='sky':
            params_dict = {'1a':[],
                           '1b':[],
                           '2a':[],
                           '2b':[],
                           '3a':[],
                           '3b':[],
                           '4a':[],
                           '4b':[]}
            
        elif mode=='image':
            params_dict = {'1a': [  0.546840728915, 2.085508651419,   755.7441976306, 1903.0543294730, -0.000116966255265],
                           '1b': [  0.278208978687, 1.575616170670,  1618.7300826916,   38.2347480120, -0.000132133772380],
                           '2a': [  0.106012209935, 1.955256534241,  1281.8575530614, -321.0904517591, -0.000095834807390],
                           '2b': [  0.265442933350, 1.695606870906,  -301.3231313075, 1809.5191721428, -0.000106731228007],
                           '3a': [  0.124723659284, 2.546896349490, -1278.5682060299,  875.1604874917, -0.000107604547391],
                           '3b': [ -0.060101127529, 1.827324548451,  1203.0067973049,  140.1525599840, -0.000106604806396],
                           '4a': [],
                           '4b': []}
    elif field=='J0148':
        # /net/galaxy-data/export/galaxydata/kashinod/EIGER/J0148/reduction_imaging/checkWCS_v2
        if mode=='sky':
            params_dict = {'1a':[207.391216171414,-21.003666696626, 27.151822573146, 5.930165941027,  0.000070159605639],
                           '1b':[152.537887986363, -7.629500354660, 27.157305853915, 5.930759735036, -0.000097545591472],
                           '2a':[ 24.121312945620,-28.988398222381, 27.129097729427, 5.936188098798, -0.000006925410051],
                           '2b':[ -1.213112092297,-12.258138682229, 27.119083606057, 5.944811531643, -0.000109424606538],
                           '3a':[153.554439743794,-16.075405206854, 27.105529660076, 5.816029325306, -0.000030558396804],
                           '3b':[163.158792445578,-24.962983704924, 27.109603537234, 5.958668473446, -0.000017158386011],
                           '4a':[165.233242217694,-11.768772125565, 27.256719473959, 5.849131353619, -0.000023250099495],
                           '4b':[ -1.340515644704, 15.385218055429, 27.122145752164, 5.677427956358, -0.000140591341944]}

        elif mode=='image':
            params_dict = {'1a': [ -1.544828286611, 2.734060484019,   953.6005885861, -1090.5244729844, -0.000069649779291],
                           '1b': [ -1.104126070655, 2.870779963096,  2211.4687114834, -1215.8840060722,  0.000072128514995],
                           '2a': [ -0.617929388718, 0.168118151111,   416.9034036184,  1999.4210346918,  0.000029696221771],
                           '2b': [ -0.345978621152, 0.440055846913,  1602.2610155278,  -201.0166753042,  0.000111243393144],
                           '3a': [ -1.294802104057, 2.298034755265, -1442.8260173397, -1352.4870468687,  0.000029836032946],
                           '3b': [ -1.459288191226, 2.350796315496,  3608.6993246073,  3075.7115174845,  0.000024962266898],
                           '4a': [ -1.165578561294, 2.545221110188,   -79.8694290227,  1678.2717340060,  0.000015429592314],
                           '4b': [ -1.110228965745, 2.428816532077,   485.2342311744,   910.8941360321,  0.000094626985098]}
    else:
        raise ValueError('field='+field+' cannot be found.')

    return params_dict


def radec_in_this_vismod(ra0, dec0, visit=1, module='a', field=None):
    params_dict = get_wcs_offset_params_dict(field, mode='sky')
    params=params_dict[str(visit)+module.lower()]
    ra, dec = offset_wcs_in_sky(ra0, dec0, params)
    return ra, dec

def radec_in_this_vismod_with_wcs(ra0, dec0, wcs, visit=1, module='a', field=None):
    params_dict = get_wcs_offset_params_dict(field, mode='image')
    
    params=params_dict[str(visit)+module.lower()]
     
    x_radec0, y_radec0 = wcs.all_world2pix(ra0, dec0, 0)
    x_radec, y_radec = offset_wcs(x_radec0, y_radec0, params)
    #print(x_radec0, y_radec0, '==>', x_radec, y_radec)
    ra_vm, dec_vm = wcs.all_pix2world(x_radec, y_radec, 0)
    return ra_vm, dec_vm

