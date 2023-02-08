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
from lmfit import Model
from scipy.integrate import simps



def gaussian(x,totflux,c,x0,sigma):
    return totflux*((sigma)**-1 * (2*np.pi)**-0.5 *np.exp(-(x-x0)**2/(2*sigma**2)))+c   


def gaussian_O3(x,a,c,redshift,sigma):
    x0_O31=(1+redshift)*4960.295
    x0_O32=(1+redshift)*5008.24
    return c+ a*(numpy.exp(-(x-x0_O31)**2/(2*sigma**2))+c  +2.98*numpy.exp(-(x-x0_O32)**2/(2*sigma**2)))

def gaussian_O3_withfudge(x,a,c,redshift,sigma,fudge):
    x0_O31=(1+redshift)*4960.295
    x0_O32=(1+redshift)*5008.24
    return c+ a*(numpy.exp(-(x-x0_O31)**2/(2*sigma**2))+c  +2.98*fudge*numpy.exp(-(x-x0_O32)**2/(2*sigma**2)))








FOLDER='/scratch/EIGER/identification/SPECTRA_O3CANDS_15SEPT/'
SPECTRA_FOLDER='/scratch/EIGER/extraction/OPTIMAL_PROFILES/'


#SPECTRA_FOLDER='/scratch/EIGER/extraction/OPTIMAL_PROFILES_SCI/'

SAVE_FOLDER='/scratch/EIGER/extraction/REDSHIFTS/'

CATALOG='/scratch/EIGER/identification/J0100_photcat_v2_CONCAT_O3candidates_HYBRID_BACKWARD_MANUALRADEC_withREDSHIFT.fits'

CATALOG='/scratch/EIGER/extraction/J0100_photcat_v2_O3emitters_Systems_28092022_FULL.fits'


SPECTRA_FOLDER='/scratch/EIGER/extraction/ALL_ONED/'

CATALOG='/scratch/EIGER/extraction/J0100_photcat_v2_O3emitters_Systems_13102022_FULL.fits'

CATALOG='/scratch/EIGER/extraction/J0100_photcat_v2_O3emitters_Systems_01112022_FULL_with_Prosp.fits'

CATALOG='/scratch/EIGER/extraction/J0100_photcat_v2_O3emitters_Systems_18112022_Prospector.fits'

with fits.open(CATALOG) as hdul:
    orig_table = hdul[1].data
    orig_cols = orig_table.columns
    
rescale_noise=True #if true, rescales the mean(err_1d) to be equal to the std(data_1d) with some outlier removal


cat=fits.open(CATALOG)

data=cat[1].data
IDlist=data.field('NUMBER')

z_fitlist=data.field('z_O3doublet_combined')
#Nclumps_Y=data.field('Nclumps_Y')
GOOD=data.field('ALL_SYSTEMS')
#logLO3=np.log10(data.field('L_O3'))

f356w=data.field('fnu_F356W_AUTO_apcor')
f200w=data.field('fnu_F200W_AUTO_apcor')
f115w=data.field('fnu_F115W_AUTO_apcor')
PHOTCONT=data.field('PHOTOMETRY_CONTAMINATED')

CONFID=data.field('CONFID_AVG')
MUV=data.field('MUV_115W')
Mstar=np.log10(data.field('mstar_50pct'))

SNHb=data.field('f_Hb')/data.field('f_Hb_err')
O3Hb=data.field('f_O3_5008')/data.field('f_Hb')

#O3Hb=data.field('f_O3_5008')/(3*data.field('f_Hb_err'))
sel=(SNHb>5)*(O3Hb<4.)


wav_stack=np.arange(4000,5500,1.)
stack=[]
Ngal=[]
F356Wstack=[]
F115Wstack=[]
F200Wstack=[]


for q in range(len(IDlist)):
	thisID=IDlist[q]
	thisz=z_fitlist[q]
	#thisL=logLO3[q]
	thisF356W=f356w[q]
	thisF200W=f200w[q]
	thisF115W=f115w[q]
	thisCRIT=PHOTCONT[q]
	thisMUV=MUV[q]
	thismass=Mstar[q]
	thisselhere=sel[q]

	print('now doing id',thisID,thisz)
	if GOOD[q]==False:
		continue
	if CONFID[q]<4.5:
		continue
	#thisNclumps_Y=Nclumps_Y[q]

	#if thisNclumps_Y>1:
	#	continue

	# if thisz<5.2 or thisz>7.89:# or thisL<42.4:
	# 	continue
	# name='full'

	# if thisz<5.5 or thisz>6.89:# or thisL<42.4:
	# 	continue
	# name='full_Hbcov_new'

	# if thisz<6.26 or thisz>6.89:# or thisL<42.4:
	# 	continue
	# name='full_Hgcov_new'

	# if thisz<6.61 or thisz>6.89:# or thisL<42.4:
	# 	continue
	# name='full_z6p8'


	# if thisz<6.3 or thisz>6.35:# or thisL<42.4:
	# 	continue
	# name='QSO'	

	# if ((thisz>6.3) * (thisz<6.35)) or thisz<6.26:
	# 	continue
	# name='nonQSO'	

	# if thisz<5.5 or thisz>6.02 or thisL<42.4:
	# 	continue
	# name='belowz6_Hbcov'		

	# if thisz<5.5 or thisL>43.1 or thisL<42.6:
	# 	continue
	# name='between42p6_43p1_Hbcov'

	# if thisz<5.5 or thisL<43.1:
	# 	continue
	# name='above43p1_Hbcov'

	# if thisz<5.5 or thisMUV>-19.5 or thisMUV<-20.2:
	# 	print('ignore',thisMUV)
	# 	continue
	# name='betweenMUV19p5_20p2_Hbcov'

	# if thisz<5.5  or thisMUV<-18.7:
	# 	print('ignore',thisMUV)
	# 	continue
	# name='belowMUV18p7_Hbcov'

	# if thisz<5.5 or thismass>9. or thismass<8.:# or thisL<42.4:
	# 	continue
	# name='midmass'

	if thisz<5.5 or thismass<7.5 or thismass>8.5:# or thisL<42.4:
		continue
	#name='mass_below7p7_new'
	name='mass_7p45_to_8p5_new'	
	#name='mass_above_9p3_new'
	# if thisz<5.5 or thisselhere==False:# or thismass>8.8:# or thisL<42.4:
	# 	continue
	# name='lowO3Hb'
	dat=fits.open(SPECTRA_FOLDER+'spectrum_1D_%s.fits'%thisID)
	data_1d=dat[1].data

	obs_wav=data_1d.field('wavelength') #in Angstroms
	sel_use=(obs_wav>31500)*(obs_wav<39500)

	thisflux=[]
	for module in ['A','B']: #Need to first coadd indiv galaxies before stacking!
		try:
			flux_tot=data_1d.field('flux_tot_%s'%module)
			flux_tot_err=data_1d.field('flux_tot_%s_err'%module)
		except:
			continue
		thisflux.append(flux_tot)

	print(np.shape(thisflux))

	flux_galaxy=np.nanmean(thisflux,axis=0)
	wav_rest=obs_wav/(1+thisz)
	flux_galaxy[~sel_use]=np.nan

	numbers=np.zeros(len(wav_rest))
	numbers[sel_use]=1.


	sel_cont=((wav_rest>4500)*(wav_rest<4820))# + ((wav_rest>5040)*(wav_rest<5400))
	med_cont=np.nanmean(flux_galaxy[sel_cont])

	#normalise
	flux_galaxy=flux_galaxy-med_cont

	DL=cosmo.luminosity_distance(thisz).value
	flux_galaxy=flux_galaxy*DL*DL * (1+thisz) 

	f=interp1d(wav_rest,flux_galaxy,kind='linear',fill_value=np.nan,bounds_error=False)
	fn=interp1d(wav_rest,numbers,kind='linear',fill_value=np.nan,bounds_error=False)

	stack.append(f(wav_stack))
	Ngal.append(fn(wav_stack))

	if thisCRIT==False:
		F356Wstack.append(thisF356W*DL*DL*(1+thisz)*1E-9*1E-23*3.0856E24 * 3.0856E24*4*numpy.pi)
		F200Wstack.append(thisF200W*DL*DL*(1+thisz)*1E-9*1E-23*3.0856E24 * 3.0856E24*4*numpy.pi)
		F115Wstack.append(thisF115W*DL*DL*(1+thisz)*1E-9*1E-23*3.0856E24 * 3.0856E24*4*numpy.pi)


F356Wstack=np.array(F356Wstack)
F200Wstack=np.array(F200Wstack)
F115Wstack=np.array(F115Wstack)

median_stack=np.nanmedian(stack,axis=0)
Nstack=np.nansum(Ngal,axis=0)

print(np.shape(stack))
indici=np.arange(0,np.shape(stack)[0],1)
print(indici)

stacks=[]


for q in range(300):
	sel=np.random.choice(indici,size=indici.shape,replace=True)
	full=[]
	for jj in range(len(sel)):
		this_index=sel[jj]
		full.append(stack[this_index])
	stacks.append(np.nanmean(full,axis=0)) ## MEAN

print(np.shape(stacks))

boot_median=np.nanmedian(stacks,axis=0)
boot_std=np.nanstd(stacks,axis=0)
print(boot_std)
#rel_trans=numpy.random.choice(trans_use,size=trans_use.shape,replace=True)


#NOW STACK PHOTOMETRY
list_f356w=[]
list_f200w=[]
list_f115w=[]
for q in range(100):
	this_selection_356w=np.random.choice(F356Wstack,size=F356Wstack.shape,replace=True)
	list_f356w.append(np.nanmedian(this_selection_356w))

	this_selection_200w=np.random.choice(F200Wstack,size=F200Wstack.shape,replace=True)
	list_f200w.append(np.nanmedian(this_selection_200w))
	
	this_selection_115w=np.random.choice(F115Wstack,size=F115Wstack.shape,replace=True)
	list_f115w.append(np.nanmedian(this_selection_115w))

boot_f356w_median=np.nanmedian(list_f356w)
boot_f200w_median=np.nanmedian(list_f200w)
boot_f115w_median=np.nanmedian(list_f115w)


boot_f356w_std=np.nanstd(list_f356w)
boot_f200w_std=np.nanstd(list_f200w)
boot_f115w_std=np.nanstd(list_f115w)

#to improve: 
#1 luminosity scaling
#2 error propagation /bootstrap
#3 sampling


col1 = fits.Column(name='restframe_wavelength', format='D', array=wav_stack)
col2 = fits.Column(name='flux_median', format='D', array=median_stack* 3.0856E24 * 3.0856E24*4*numpy.pi*1E-18)
col3 = fits.Column(name='flux_median_boot', format='D', array=boot_median* 3.0856E24 * 3.0856E24*4*numpy.pi*1E-18)
col4 = fits.Column(name='flux_median_boot_err', format='D', array=boot_std* 3.0856E24 * 3.0856E24*4*numpy.pi*1E-18)
col5 = fits.Column(name='Ngal', format='D', array=Nstack)

col6 = fits.Column(name='median_F356W_ergsHz', format='D', array=np.zeros(len(wav_stack))+boot_f356w_median)
col7 = fits.Column(name='median_F356W_ergsHz_err', format='D', array=np.zeros(len(wav_stack))+boot_f356w_std)
col8 = fits.Column(name='median_F200W_ergsHz', format='D', array=np.zeros(len(wav_stack))+boot_f200w_median)
col9 = fits.Column(name='median_F200W_ergsHz_err', format='D', array=np.zeros(len(wav_stack))+boot_f200w_std)
col10 = fits.Column(name='median_F115W_ergsHz', format='D', array=np.zeros(len(wav_stack))+boot_f115w_median)
col11= fits.Column(name='median_F115W_ergsHz_err', format='D', array=np.zeros(len(wav_stack))+boot_f115w_std)

new_cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11])
hdu = fits.BinTableHDU.from_columns(new_cols)

hdu.writeto('stacked_1D_EMLINE_median_%s.fits'%name, overwrite=True)


