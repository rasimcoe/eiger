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


def get_obs_mag(filename,lamb,fluxdensity):
    d110=numpy.loadtxt(filename)
    c=2.998E18 #Angstrom/s

    wav_110=d110[:,0] #angstrom
    trans_110=d110[:,1]
    trans_110=snd.gaussian_filter(trans_110,sigma=3)


    dummy=np.arange(1.0,4.1,0.001) * 1E4
    trans_110=trans_110/numpy.nanmax(trans_110)
    interp_transmission=interp1d(wav_110,trans_110,kind='linear',fill_value='extrapolate')
    interp_flux  = interp1d(lamb,fluxdensity,kind='cubic')           
    S=interp_flux(dummy)
    T=interp_transmission(dummy)
    #sel_zero=(lamb<3E4) + (lamb>4E4)
    lam=dummy

    I1        = simps(S*T*lam,lam)                     #Denominator
    I2        = simps(  T/lam,lam)                     #Numerator
    fnu       = I1/I2 / c                         #Average flux density

    mAB       = -2.5*np.log10(fnu) - 48.6              #AB magnitude

    return mAB


def continuum_and_lines(dummyx,norm,beta,EW_Hb,EW_O3,redshift): #norm is in erg/s/cm2/AA at lambda_0=5000 A
    norm=1E-20*norm

    lamb=numpy.arange(1000,7000,0.05)*(1+redshift) 
    cont=norm*(lamb/(5000.*(1+redshift)))**beta
    lines=numpy.zeros(numpy.shape(lamb))

    sel_Hb_wav=(lamb<((4862.69)*(1+redshift)+5.))*(lamb>((4862.69)*(1+redshift)-5.))
    sel_cont_Hb_wav=(lamb<((4862.69)*(1+redshift)+10.))*(lamb>((4862.69)*(1+redshift)-10.))

    sel_O3_4960_wav=(lamb<((4960)*(1+redshift)+5.))*(lamb>((4960)*(1+redshift)-5.))
    sel_cont_O3_4960_wav=(lamb<((4960)*(1+redshift)+10.))*(lamb>((4960)*(1+redshift)-10.))

    sel_O3_5008_wav=(lamb<((5008)*(1+redshift)+5.))*(lamb>((5008)*(1+redshift)-5.))
    sel_cont_O3_5008_wav=(lamb<((5008)*(1+redshift)+10.))*(lamb>((5008)*(1+redshift)-10.))

    sel_Hg_wav=(lamb<((4340)*(1+redshift)+5.))*(lamb>((4340)*(1+redshift)-5.))
    sel_O3_4363_wav=(lamb<((4363)*(1+redshift)+5.))*(lamb>((4363)*(1+redshift)-5.))


    sel_5000=(lamb<((5000)*(1+redshift)+15.))*(lamb>((5000)*(1+redshift)-15.))

    mean_cont_Hb=numpy.nanmean(cont[sel_cont_Hb_wav])
    mean_cont_O3_4960=numpy.nanmean(cont[sel_cont_O3_4960_wav])
    mean_cont_O3_5008=numpy.nanmean(cont[sel_cont_O3_5008_wav])
    mean_5000=numpy.nanmean(cont[sel_5000])

    fnu_5000=mean_5000*(3.34E4*(5000.*(1+redshift))**2) #erg/s/cm2/A

    Hb_flux=EW_Hb*(1+redshift)*mean_cont_Hb
    O3_flux=0.25*EW_O3*(1+redshift)*mean_cont_O3_4960 + 0.75*EW_O3*(1+redshift)*mean_cont_O3_5008


    npixels=len(lines[sel_Hb_wav])
    lamb_per_pixel=lamb[1]-lamb[0]
    lines[sel_Hb_wav]+=Hb_flux/(npixels*lamb_per_pixel)

    npixels=len(lines[sel_O3_4960_wav])
    lines[sel_O3_4960_wav]+=0.25*O3_flux/(npixels*lamb_per_pixel)

    npixels=len(lines[sel_O3_5008_wav])
    lines[sel_O3_5008_wav]+=0.75*O3_flux/(npixels*lamb_per_pixel)

    #add Hgamma and O3 4363
    npixels=len(lines[sel_Hg_wav])

    lines[sel_Hg_wav]+=0.4*Hb_flux/(npixels*lamb_per_pixel)
    lines[sel_O3_4363_wav]+=0.03*0.75*O3_flux/(npixels*lamb_per_pixel)


    sedmodel=lines+cont
    fn_F356W=10**(-0.4*(get_obs_mag('/scratch/EIGER/extraction/JWST_NIRCam.F356W.dat',lamb,sedmodel)-8.9)) * 1E9 #nJy

    fn_F200W=10**(-0.4*(get_obs_mag('/scratch/EIGER/extraction/JWST_NIRCam.F200W.dat',lamb,sedmodel)-8.9)) * 1E9 #nJy
    fn_F115W=10**(-0.4*(get_obs_mag('/scratch/EIGER/extraction/JWST_NIRCam.F115W.dat',lamb,sedmodel)-8.9)) * 1E9 #nJy

    #return numpy.array([fn_F115W,fn_F200W,fn_F356W,Hb_flux*1E18,O3_flux*1E18])

    return numpy.array([fn_F115W,fn_F200W,Hb_flux*1E18,O3_flux*1E18])






CATALOG='/scratch/EIGER/extraction/testcat.fits'
SAVE_CATALOG='/scratch/EIGER/extraction/testcat_withEW.fits'


CATALOG='/scratch/EIGER/extraction/J0100_photcat_v2_O3emitters_Systems_06102022_FULL_with_EW.fits'
SAVE_CATALOG='/scratch/EIGER/extraction/J0100_photcat_v2_O3emitters_Systems_06102022_FULL_with_EW_withBetaLQ.fits'

cat=fits.open(CATALOG)


with fits.open(CATALOG) as hdul:
    orig_table = hdul[1].data
    orig_cols = orig_table.columns
    

data=cat[1].data
IDlist=data.field('NUMBER')

z_guesslist=data.field('z_O3doublet_combined')

#erg/s/cm2
f_Hb=data.field('f_Hb')#  * 1E-18
f_Hb_err=data.field('f_Hb_err')# * 1E-18
f_O3_4960=data.field('f_O3_4960')# * 1E-18
f_O3_4960_err=data.field('f_O3_4960_err')# * 1E-18
f_O3_5008=data.field('f_O3_5008')# * 1E-18
f_O3_5008_err=data.field('f_O3_5008_err') #* 1E-18


#nJy
fnu_F115W=data.field('fnu_F115W_AUTO_apcor')
fnu_F115W_err=data.field('enu_F115W_aper_model')#enu_F356W_aper_model


fnu_F200W=data.field('fnu_F200W_AUTO_apcor')
fnu_F200W_err=data.field('enu_F200W_aper_model')


fnu_F356W=data.field('fnu_F356W_AUTO_apcor')
fnu_F356W_err=data.field('enu_F356W_aper_model')



LIST_EW_O3=np.zeros(len(IDlist))
LIST_EW_Hb=np.zeros(len(IDlist))
LIST_F356W_nolines=np.zeros(len(IDlist))



LIST_EW_O3_err=np.zeros(len(IDlist))
LIST_EW_Hb_err=np.zeros(len(IDlist))
LIST_F356W_nolines_err=np.zeros(len(IDlist))

BETA=np.zeros(len(IDlist))
BETA_ERR=np.zeros(len(IDlist))

for q in range(len(IDlist)):
    try:
        thisID=IDlist[q]
        thisz=z_guesslist[q]
        thisfnu_F356W=fnu_F356W[q]
        thisfnu_F356W_err=fnu_F356W_err[q] + 0.05*fnu_F356W[q]
        thisfnu_F200W=fnu_F200W[q]
        thisfnu_F200W_err=fnu_F200W_err[q] + 0.1*fnu_F200W[q]

        thisfnu_F115W=fnu_F115W[q]
        thisfnu_F115W_err=fnu_F115W_err[q] + 0.1*fnu_F115W[q]

        thisHb_flux=f_Hb[q]
        thisHb_flux_err=f_Hb_err[q]
        thisO3_flux=f_O3_4960[q] + f_O3_5008[q]
        thisO3_flux_err=f_O3_4960_err[q] + f_O3_5008_err[q]



        this_flux_array=np.array([thisfnu_F115W,thisfnu_F200W,thisfnu_F356W,thisHb_flux,thisO3_flux])
        this_flux_err_array=np.array([thisfnu_F115W_err,thisfnu_F200W_err,thisfnu_F356W_err,thisHb_flux_err,thisO3_flux_err]) #+ 0.05*this_flux_array #Adding 5% of the flux as noise, e.g. flux calibration uncertainties


        this_flux_array=np.array([thisfnu_F115W,thisfnu_F200W,thisHb_flux,thisO3_flux])
        this_flux_err_array=np.array([thisfnu_F115W_err,thisfnu_F200W_err,thisHb_flux_err,thisO3_flux_err]) #+ 0.05*this_flux_array #Adding 5% of the flux as noise, e.g. flux calibration uncertainties

        print('ID',thisID)
        print('Observed',this_flux_array,this_flux_err_array)

        dx=2. #this is just a dummy


        model = Model(continuum_and_lines)  #defdef skewed_gaussian(x,a,x0,asym,d):
        model.set_param_hint('norm',min=0.0001,max=50.) #a_O3,c1,redshift,sigma
    #    model.set_param_hint('beta',min=-3.3,max=0.) #this  -3.3 is motivated by BPASS models
        model.set_param_hint('beta',min=-4.5,max=+1) #this  -3.3 is motivated by BPASS models

        model.set_param_hint('EW_Hb',min=2.,max=200000.)
        model.set_param_hint('EW_O3',min=2.,max=200000.)
        model.set_param_hint('redshift',min=thisz-0.005,max=thisz+0.005)


        params = model.make_params(norm=1.,beta=-2.1,EW_Hb=100.,EW_O3=300.,redshift=thisz) #a=12., x0=12, asym=0.25, d=80.,a1=12.,peaksep=200.,asym1=0.25,d1=80
        params['redshift'].vary=False
       # params['beta'].vary=False

        result = model.fit(this_flux_array, params,dummyx=dx,nan_policy='raise',weights=1./this_flux_err_array,scale_covar=False,method='least_squares')

        print(result.fit_report())
        print('Model',model.eval(result.params,dummyx=dx))

        cleaned_mag=[]
        lamb=numpy.arange(1000,7000,0.05)*(1+thisz) 

        bestn=result.params['norm'].value
        stdn=result.params['norm'].stderr
        print(bestn,stdn)

        for i in range(30):
        	thisn=np.random.normal(bestn,stdn) *1E-20
        	cont=thisn*(lamb/(5000.*(1+thisz)))**-2.
        	fn_F356W_nolines=10**(-0.4*(get_obs_mag('/scratch/EIGER/extraction/JWST_NIRCam.F356W.dat',lamb,cont)-8.9)) * 1E9 #nJy
        	cleaned_mag.append(fn_F356W_nolines)
        print(np.nanmedian(cleaned_mag),np.nanstd(cleaned_mag))



        LIST_EW_O3[q]=result.params['EW_O3'].value
        LIST_EW_Hb[q]=result.params['EW_Hb'].value
        LIST_F356W_nolines[q]=np.nanmedian(cleaned_mag)


        LIST_EW_O3_err[q]=result.params['EW_O3'].stderr
        LIST_EW_Hb_err[q]=result.params['EW_Hb'].stderr
        LIST_F356W_nolines_err[q]=np.nanstd(cleaned_mag)

        BETA[q]=result.params['beta'].value
        BETA_ERR[q]=result.params['beta'].stderr
    except:
        continue

col1 = fits.Column(name='fnu_F356W_AUTO_apcor_nolines', format='D', array=LIST_F356W_nolines)
col2 = fits.Column(name='enu_F356W_AUTO_apcor_nolines', format='D', array=LIST_F356W_nolines_err)

col3 = fits.Column(name='EW0_O3doublet', format='D', array=LIST_EW_O3)
col4 = fits.Column(name='EW0_O3doublet_err', format='D', array=LIST_EW_O3_err)

col5 = fits.Column(name='EW0_Hb', format='D', array=LIST_EW_Hb)
col6 = fits.Column(name='EW0_Hb_err', format='D', array=LIST_EW_Hb_err)

col7 = fits.Column(name='beta', format='D', array=BETA)
col8 = fits.Column(name='beta_err', format='D', array=BETA_ERR)


new_cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8])
hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)

hdu.writeto(SAVE_CATALOG, overwrite=True)
