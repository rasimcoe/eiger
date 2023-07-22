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

SPECTRA_FOLDER='/scratch/EIGER/identification/SPECTRA_O3CANDS_15SEPT/'

CATALOG='/scratch/EIGER/extraction/J0100_photcat_v2_O3emitters_Systems_13102022_FULL.fits'

CATALOG='/scratch/EIGER/extraction/J0100_photcat_v2_O3emitters_Systems_20102022_FULL.fits'


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
logLO3=np.log10(data.field('L_O3'))

CONTAM=data.field('GRISM_CONTAMINATED')

f356w=data.field('fnu_F356W_AUTO_apcor')
f200w=data.field('fnu_F200W_AUTO_apcor')
f115w=data.field('fnu_F115W_AUTO_apcor')
PHOTCONT=data.field('PHOTOMETRY_CONTAMINATED')

MUV=data.field('MUV_115W')
MUV=data.field('MUV_prosp')
CONFID=data.field('CONFID_AVG')
wav_stack=np.arange(4300,5500,1.)
stack=[]
Ngal=[]
F356Wstack=[]
F115Wstack=[]
F200Wstack=[]


for q in range(len(IDlist)):
	thisID=IDlist[q]
	thisz=z_fitlist[q]
	thisL=logLO3[q]
	thisF356W=f356w[q]
	thisF200W=f200w[q]
	thisF115W=f115w[q]
	thisCRIT=PHOTCONT[q]
	thisMUV=MUV[q]


	print('now doing id',thisID,thisz)
	if GOOD[q]==False:
		continue
	if CONTAM[q]==True:
		continue
	if CONFID[q]<4.5:
		continue
	#thisNclumps_Y=Nclumps_Y[q]

	#if thisNclumps_Y>1:
	#	continue

	# if thisz<5.4 or thisz>6.99:# or thisL<42.4:
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

	if thisz<5.5 or thisz>6.95:# or thisMUV>-20.5:#thisF356W<200:#thisMUV>-20.:# or thisL<42.4:
		continue
	name='bright_n'
	name='all'

	hdu= fits.open(SPECTRA_FOLDER+'v4corstacked_2D_%s.fits'%thisID)

	hd=hdu['EMLINE'].header
	data=hdu['SCI'].data

	ly,lx=np.shape(data)

	x_array=np.arange(0,lx,1)
	wav_array=x_array*hd['CDELT1'] + hd['CRVAL1']

	obs_wav=wav_array
	sel_use=(obs_wav>31000)*(obs_wav<39900)


	registered=np.zeros((ly,len(wav_stack)))
	wav_rest=obs_wav/(1+thisz)
	data[:,~sel_use]=np.nan



	#IDENTIFY STUFF
	std=numpy.nanstd(data[:,sel_use])
	sel_use2=(wav_rest>4700)*(wav_rest<5100) *(data**2 < 4*std**2)
	print(std)
	std=numpy.nanstd(data[sel_use2])
	sel_use2=(wav_rest>4700)*(wav_rest<5100)*(data**2 < 4*std**2)
	print(std)
	std=numpy.nanstd(data[sel_use2])
	sel_use2=(wav_rest>4700)*(wav_rest<5100) *(data**2 < 4*std**2)
	print(std)
	std=numpy.nanstd(data[sel_use2])
	sel_use2=(wav_rest>4700)*(wav_rest<5100) *(data**2 < 3*std**2)
	print(std)
	std=numpy.nanstd(data[sel_use2])
	sel_use2=(wav_rest>4700)*(wav_rest<5100)*(data**2 < 3*std**2)
	print(std)
	std=numpy.nanstd(data[sel_use2])
	sel_use2=(wav_rest>4700)*(wav_rest<5100) *(data**2 < 3*std**2)
	print(std)
	std=numpy.nanstd(data[sel_use2])
	med=numpy.nanmedian(data[sel_use2])
	sel_use2=(wav_rest>4700)*(wav_rest<5100) *(data**2 < 3*std**2)
	print(std)					
	
	sel_mask=data>5*std
	before_data=copy.deepcopy(data)
	data[sel_mask]=numpy.nan
	data[20:30,:]=before_data[20:30,:]

	#numbers=np.zeros(len(wav_rest))
	#numbers[sel_use]=1.

	#sel_cont=((wav_rest>4500)*(wav_rest<4820))# + ((wav_rest>5040)*(wav_rest<5400))
	#med_cont=np.nanmean(flux_galaxy[sel_cont])

	#normalise
	#flux_galaxy=flux_galaxy-med_cont

	DL=cosmo.luminosity_distance(thisz).value
	flux_galaxy=data*DL*DL * (1+thisz) *thisF200W/numpy.nanmedian(f200w)

	ff=copy.deepcopy(flux_galaxy)
	ff[19:33,:]=numpy.nan
	ff[:12,:]=numpy.nan
	ff[40:,:]=numpy.nan
	medianlevel=numpy.nanmedian(ff,axis=0)


	for q in range(ly):
		f=interp1d(wav_rest,flux_galaxy[q,:] - 0.0*medianlevel,kind='linear',fill_value=np.nan,bounds_error=False)
		registered[q,:]=f(wav_stack)
	#fn=interp1d(wav_rest,numbers,kind='linear',fill_value=np.nan,bounds_error=False)

	stack.append(registered)
	#Ngal.append(fn(wav_stack))

median_stack=np.nanmedian(stack,axis=0)
fits.writeto('2Dstacktest.fits',median_stack,overwrite=True)

#Nstack=np.nansum(Ngal,axis=0)

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
	stacks.append(np.nanmedian(full,axis=0))

print(np.shape(stacks))

boot_median=np.nanmedian(stacks,axis=0)
boot_std=np.nanstd(stacks,axis=0)
print(boot_std,np.shape(boot_std))
#rel_trans=numpy.random.choice(trans_use,size=trans_use.shape,replace=True)

fits.writeto('2Dstacktest_bootmed.fits',boot_median,overwrite=True)
fits.writeto('2Dstacktest_bootstd.fits',boot_std,overwrite=True)
fits.writeto('2Dstacktest_bootSN.fits',boot_median/boot_std,overwrite=True)

#to improve: 
#1 luminosity scaling
#2 error propagation /bootstrap
#3 sampling
print('N gals',len(indici))

#Normalisation
ff=copy.deepcopy(boot_median)
ff[20:30,:]=numpy.nan
#ff[:3,:]=numpy.nan
#ff[47:,:]=numpy.nan

ff[20:30,:]=numpy.nan
ff[:8,:]=numpy.nan
ff[42:,:]=numpy.nan
medianlevel=numpy.nanmedian(ff,axis=0)
medianlevel_std=numpy.nanstd(ff,axis=0)/len(indici)**0.5
ff=ff-medianlevel

print('-->',np.shape(medianlevel),np.shape(ff))
print('MEDIANLEVEL',numpy.nanmedian(medianlevel)/1E8)

boot_median=boot_median-medianlevel
boot_std=(boot_std**2 + medianlevel_std**2)**0.5
subset_data=(boot_median[10:40,656:664])/ 1E9
collapse_y=np.nansum(subset_data,axis=1)


x_fit=np.arange(0,len(collapse_y),1)
x_fit_oversample=np.arange(0,len(collapse_y),0.01)

model=Model(gaussian,independent_vars=('x'),prefix='m1_')
model.set_param_hint('m1_totflux',min=0.,max=40)
model.set_param_hint('m1_sigma',min=0.8,max=4.)
model.set_param_hint('m1_c',min=-2,max=2)
model.set_param_hint('m1_x0',min=21.-10,max=29.-10)

params=model.make_params(totflux=0.5,c=0,x0=26-10,sigma=2)


result=model.fit(collapse_y,x=x_fit,params=params)

print(result.fit_report())


fig, (ax1, ax2) = pyplot.subplots(1, 2)
ax1.imshow(subset_data,origin='lower')
ax2.plot(x_fit,collapse_y)
ax2.plot(x_fit_oversample,model.eval(result.params,x=x_fit_oversample))
pyplot.savefig('opt_profile_2dstack.png')
pyplot.clf()


result.params['m1_c'].value=0.

model_y=model.eval(result.params,x=x_fit)
model_y=model_y/np.nansum(model_y)
opt_weight=np.zeros(np.shape(boot_median))

for j in range(len(opt_weight[0,:])):
	opt_weight[10:40,j]=model_y





t1 = np.nansum(boot_median*opt_weight,axis=0)
t2 = np.nansum(opt_weight**2,axis=0)
sci_extracted = t1/t2	

ivar=boot_std**-2

t1 = np.nansum(ivar*opt_weight**2,axis=0)
t2 = np.nansum(opt_weight,axis=0)
sci_extracted_err = (t1/t2)**-0.5



col1 = fits.Column(name='restframe_wavelength', format='D', array=wav_stack)
col2 = fits.Column(name='flux_median_boot', format='D', array=sci_extracted* 3.0856E24 * 3.0856E24*4*numpy.pi*1E-18)
col3 = fits.Column(name='flux_median_boot_err', format='D', array=sci_extracted_err* 3.0856E24 * 3.0856E24*4*numpy.pi*1E-18)

new_cols = fits.ColDefs([col1,col2,col3])
hdu = fits.BinTableHDU.from_columns(new_cols)

hdu.writeto('stacked_1D_from2D_SCI_median_%s_O3ap_masking.fits'%name, overwrite=True)


