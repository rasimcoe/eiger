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
from astropy.coordinates import SkyCoord


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



def stack_with_reject_outliers(data,err,nobs, m = 5.):
    median=np.nanmedian(data,axis=0)
    std=np.nanstd(data,axis=0)
    mask_outliers=np.abs((data-median))/std > m
    data[mask_outliers]=np.nan
    var=err**2 
    var[mask_outliers]=np.nan
    #return np.nanmean(data,axis=0), np.nanmean(np.nansum(var,axis=0))**0.5

    return np.nansum(data,axis=0)/nobs,(np.nansum(var,axis=0)/nobs)**0.5

    
    
def stack_with_reject_outliers_median(data, m = 5.):
    median=np.nanmedian(data,axis=0)
    std=np.nanstd(data,axis=0)
    mask_outliers=np.abs((data-median))/std > m
    data[mask_outliers]=np.nan

    return np.nanmedian(data,axis=0)



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
    
    hdul = fits.open(emfile) ##emfile for COLA1 
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

def radec_in_this_vismod(ra0, dec0, wcs, visit=1, module='a', field=None):
    
#     fil = '/scratch/kashinod/EIGER/J0100/reduction_imaging/checkWCS_v2/wcs_offset_params_visit'+str(visit)+module.lower()+'.txt'
#     params = np.loadtxt(fil)
    if field=='J0100':
        params_dict = {'1a':[ -0.679936091154,  0.523226446148, -1118.4604275112,  608.4828378726, 0.000558578932507],
                       '1b':[ -1.111210284948,  3.129704853141,   928.2162884476, 1283.2376923699, 0.000423325155808],
                       '2a':[ -1.084928574764,  1.745666945788,   750.9150277642, 1696.6610390372, 0.000428338992284],
                       '2b':[ -0.937111596682,  3.339523725388,   918.3636725379, 1240.0037692349, 0.000457844502586],
                       '3a':[  0.057536995593, -0.442108779774,   819.2640176016,  835.9240175748, 0.000182900260607],
                       '3b':[ -0.292402520885,  0.607881819359,  2168.5716962952, 2178.7603750517, 0.000290237599955],
                       '4a':[ -0.142066614727, -0.457647901154,   140.7245984545, 1392.8826646429, 0.000109932123407],
                       '4b':[  0.001534739317,  0.659164131917,  2793.8799930654,  881.4838458998, 0.000250044843380]}
    else:
        raise ValueError('field='+field+' cannot be found.')
    
    params=params_dict[str(visit)+module.lower()]
     
    x_radec0, y_radec0 = wcs.all_world2pix(ra0, dec0, 0)
    x_radec, y_radec = offset_wcs(x_radec0, y_radec0, params)
    #print(x_radec0, y_radec0, '==>', x_radec, y_radec)
    ra_vm, dec_vm = wcs.all_pix2world(x_radec, y_radec, 0)
    return ra_vm, dec_vm
