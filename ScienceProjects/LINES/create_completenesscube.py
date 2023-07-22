import numpy
from astropy.io import fits as pyfits #no need to call it pyfits, sorry about that
from matplotlib import pyplot
import numpy
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
from scipy.interpolate import interp1d
from astropy.convolution import Gaussian2DKernel
import scipy.ndimage as snd
from astropy.convolution import convolve
import copy

import matplotlib
import numpy as np

from astropy.io import fits


import lmfit
from lmfit.models import Gaussian2dModel
from lmfit import Model
from scipy.integrate import simps




CATALOG='/scratch/EIGER/completeness/j0100_all_lines.fits'

#j0100_all_lines_det0p9.fits

cat=fits.open(CATALOG)


with fits.open(CATALOG) as hdul:
    orig_table = hdul[1].data
    orig_cols = orig_table.columns
    

data=cat[1].data
X_IMAGE=data.field('X_IMAGE')
Y_IMAGE=data.field('Y_IMAGE')
LAMB_MIN=data.field('LAMB_MIN')
LAMB_MAX=data.field('LAMB_MAX')
LAMBDA=data.field('LAMBDA')


FLUX=data.field('FLUX') * 9.75 #to 1E-18 erg/s/cm2
SN=data.field('SN')
ERR=FLUX/SN

#Create grid:
binwidth=150.
Xrange=np.arange(0,14000,binwidth)
Yrange=np.arange(0,7300,binwidth)

lbinwidth=0.5
Lrange=np.arange(3.15,4.05,lbinwidth)
print(Xrange,Yrange,Lrange)
print(len(Xrange))


# lbinwidth=0.4
# Lrange=[3.6]


completeness_map=np.zeros((len(Lrange),len(Yrange),len(Xrange)))+25
number_map=np.zeros((len(Lrange),len(Yrange),len(Xrange)))+25

FLUX_CUBE=0.5 #5E-18


for q in range(len(Xrange)):
    for p in range(len(Yrange)):
        thisX=Xrange[q]
        thisY=Yrange[p]  
        # binwidth_X=200 + 1.5*(thisX-7000)**2/1E5
        # binwidth_Y=250 + 1.5*(thisY-3150)**2/1E5       
        binwidth_X=150 + 1.25*(thisX-7000)**2/1E5
        binwidth_Y=200 + 1.25*(thisY-3150)**2/1E5   

        print('BINWI',binwidth_X,binwidth_Y)          
        sel_patch=(X_IMAGE>thisX-binwidth_X)*(X_IMAGE<thisX+binwidth_X)*(Y_IMAGE>thisY-binwidth_Y)*(Y_IMAGE<thisY+binwidth_Y)
        minimum_minwav=np.nanmedian(LAMB_MIN[sel_patch])/1E4
        for l in range(len(Lrange)):
            thisL=Lrange[l]
            if minimum_minwav>thisL:
                completeness_patch=np.nan
                total_num=0
            else:
                sel_data=(X_IMAGE>thisX-binwidth_X)*(X_IMAGE<thisX+binwidth_X)*(Y_IMAGE>thisY-binwidth_Y)*(Y_IMAGE<thisY+binwidth_Y)*(LAMBDA>thisL-2*thisL**-2)*(LAMBDA<thisL+2*thisL**-2) #(LAMBDA>thisL-lbinwidth/1.3)*(LAMBDA<thisL+lbinwidth/1.3)
                these_fluxes=FLUX[sel_data]
                these_SN=SN[sel_data]

                sel_flux=these_fluxes>FLUX_CUBE
                sel_sn3=(these_fluxes>FLUX_CUBE ) *(these_SN>3.)
                len_tot=len(these_fluxes[sel_flux])
                len_SN3=len(these_fluxes[sel_sn3])
                total_num=len_tot

                if len_SN3>5:
                    argsort=np.argsort(these_fluxes[sel_sn3])
                    #completeness_patch=np.nanmin(these_fluxes[sel_sn3])
                    completeness_patch=np.nanmedian(these_fluxes[sel_sn3][argsort][:5])


                else:
                    completeness_patch=np.nan
                print('Patch',q,p,completeness_patch,len_tot,total_num)
            completeness_map[l,p,q]=completeness_patch
            number_map[l,p,q]=total_num


hdu = fits.PrimaryHDU()
hdu.header['NAXIS']=3
hdu.header['CTYPE3']='WAVELENGTH'
hdu.header['CUNIT3']='Micron'
hdu.header['CRPIX1']=1.0
hdu.header['CRPIX2']=1.0
hdu.header['CRPIX3']=1.0

hdu.header['CRVAL1']=Xrange[0]
hdu.header['CDELT1']=binwidth
hdu.header['CRVAL2']=Yrange[0]
hdu.header['CDELT2']=binwidth

hdu.header['CRVAL3']=Lrange[0]
hdu.header['CDELT3']=lbinwidth

fits.writeto('test_limitingflux_evenlessoverlap.fits',completeness_map,hdu.header,overwrite=True)
fits.writeto('test_cnumber_evenlessoverlap.fits',number_map,overwrite=True)


