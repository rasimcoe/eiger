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


def rotated_2dgauss(x, y, A, x0, y0, sigma_x, sigma_y, theta):
    theta = np.radians(theta)
    sigx2 = sigma_x**2; sigy2 = sigma_y**2
    a = np.cos(theta)**2/(2*sigx2) + np.sin(theta)**2/(2*sigy2)
    b = np.sin(theta)**2/(2*sigx2) + np.cos(theta)**2/(2*sigy2)
    c = np.sin(2*theta)/(4*sigx2) - np.sin(2*theta)/(4*sigy2)
    
    expo = -a*(x-x0)**2 - b*(y-y0)**2 - 2*c*(x-x0)*(y-y0)
    return A*np.exp(expo)    

def gaussian(x,totflux,c,x0,sigma):
    return totflux*((sigma)**-1 * (2*np.pi)**-0.5 *np.exp(-(x-x0)**2/(2*sigma**2)))+c   



FOLDER='/scratch/EIGER/identification/SPECTRA_O3CANDS_15SEPT/'
SAVE_FOLDER='/scratch/EIGER/extraction/OPTIMAL_PROFILES_SCI/'
CATALOG='/scratch/EIGER/identification/J0100_photcat_v2_CONCAT_O3candidates_HYBRID_BACKWARD_MANUALRADEC.fits'


rescale_noise=True #if true, rescales the mean(err_1d) to be equal to the std(data_1d) with some outlier removal


cat=fits.open(CATALOG)

data=cat[1].data
IDlist=data.field('NUMBER')

z_guesslist=data.field('z_O3_estimate')
Nclumps_Y=data.field('Nclumps_Y')
GOOD=data.field('GOOD_CANDIDATE')

#Nclumps_Y=np.zeros(len(IDlist))+1 ##The idea is to have this in the catalog to determine how many clumps *are separated in the Y direction*, this will help the fitting


#still under heavy development, I will break up the code in functions asap


for q in range(len(IDlist)):
	thisID=IDlist[q]
	thisz=z_guesslist[q]
	print('now doing id',thisID)
	if GOOD[q]==False:
		continue
	thisNclumps_Y=Nclumps_Y[q]

	#thisNclumps_Y=3.
	#thisz=6.32724223 #this is for ID 4784
	hdu= fits.open(FOLDER+'v5corstacked_2D_%s.fits'%thisID)
	hd=hdu['EMLINE'].header
	data=hdu['EMLINE'].data

	ly,lx=np.shape(data)

	x_array=np.arange(0,lx,1)
	wav_array=x_array*hd['CDELT1'] + hd['CRVAL1']


	COLS=[]
	COLS.append(fits.Column(name='wavelength',unit='angstrom',format='E',array=wav_array))
	opt_weight=np.zeros(np.shape(data))


	for module in ['A','B']:
		data=hdu['EMLINE%s'%module].data
		errdata=hdu['ERR'].data**0.5



		#select region in 2D that contains the 5008 line ##Could also generalise to any other line 
		if thisNclumps_Y==1.:
			sel_wav=(wav_array>5008.24*(1+thisz)- 5008.24*(1+thisz)*300/3E5)*(wav_array<5008.24*(1+thisz)+ 5008.24*(1+thisz)*300/3E5) #this
		if thisNclumps_Y>1.:
			sel_wav=(wav_array>5008.24*(1+thisz)- 5008.24*(1+thisz)*1000/3E5)*(wav_array<5008.24*(1+thisz)+ 5008.24*(1+thisz)*1000/3E5) #this
				

		subset_data=data[:,sel_wav]
		sel_real=np.isfinite(subset_data)
		print(len(subset_data[sel_real]))
		if len(subset_data[sel_real])==0:
			opt_weight=np.zeros(np.shape(data))
			continue


		collapse_y=np.nansum(subset_data,axis=1)


		x_fit=np.arange(0,len(collapse_y),1)
		x_fit_oversample=np.arange(0,len(collapse_y),0.01)


		if thisNclumps_Y==1.:
			model=Model(gaussian,independent_vars=('x'),prefix='m1_')
			model.set_param_hint('m1_totflux',min=0.,max=20)
			model.set_param_hint('m1_sigma',min=0.8,max=4.)
			model.set_param_hint('m1_c',min=-0.5,max=0.5)
			model.set_param_hint('m1_x0',min=21.,max=29.)

			params=model.make_params(totflux=0.5,c=0,x0=26,sigma=2)


		if thisNclumps_Y==2.:
			model=Model(gaussian,independent_vars=('x'),prefix='m1_')+Model(gaussian,independent_vars=('x'),prefix='m2_')
			#print(model.param_names)
			model.set_param_hint('m1_sigma',min=0.8,max=3.5)
			model.set_param_hint('m2_sigma',min=0.8,max=3.5)

			model.set_param_hint('m1_totflux',min=0.,max=20)
			model.set_param_hint('m2_totflux',min=0.,max=20)

			model.set_param_hint('m1_c',min=-0.5,max=0.5)

			model.set_param_hint('m1_x0',min=5.,max=29)
			model.set_param_hint('m2_x0',min=21.,max=36)

			params=model.make_params(m1_totflux=1,m1_c=0,m1_x0=22,m1_sigma=2,m2_totflux=2.,m2_c=0,m2_x0=26,m2_sigma=2)
			params['m2_c'].vary=False


		if thisNclumps_Y==3.:
			model=Model(gaussian,independent_vars=('x'),prefix='m1_')+Model(gaussian,independent_vars=('x'),prefix='m2_')+Model(gaussian,independent_vars=('x'),prefix='m3_')
			#print(model.param_names)
			model.set_param_hint('m1_sigma',min=0.8,max=3.5)
			model.set_param_hint('m2_sigma',min=0.8,max=3.5)
			model.set_param_hint('m3_sigma',min=0.8,max=3.5)

			model.set_param_hint('m1_totflux',min=0.,max=20)
			model.set_param_hint('m2_totflux',min=0.,max=20)
			model.set_param_hint('m3_totflux',min=0.,max=20)

			model.set_param_hint('m1_c',min=-0.5,max=0.5)

			model.set_param_hint('m1_x0',min=22.,max=28)
			model.set_param_hint('m2_x0',min=5.,max=21)
			model.set_param_hint('m3_x0',min=26.,max=38)


			params=model.make_params(m1_totflux=2,m1_c=0,m1_x0=26,m1_sigma=2,m2_totflux=0.2,m2_c=0,m2_x0=24,m2_sigma=2,m3_totflux=0.2,m3_c=0,m3_x0=28,m3_sigma=2)

			params['m2_c'].vary=False
			params['m3_c'].vary=False


		if thisNclumps_Y==4.: #honestly this is only used for one object, id 19021
			model=Model(gaussian,independent_vars=('x'),prefix='m1_')+Model(gaussian,independent_vars=('x'),prefix='m2_')+Model(gaussian,independent_vars=('x'),prefix='m3_')+Model(gaussian,independent_vars=('x'),prefix='m4_')
			#print(model.param_names)
			model.set_param_hint('m1_sigma',min=0.8,max=3.5)
			model.set_param_hint('m2_sigma',min=0.8,max=3.5)
			model.set_param_hint('m3_sigma',min=0.8,max=3.5)
			model.set_param_hint('m4_sigma',min=0.8,max=3.5)

			model.set_param_hint('m1_totflux',min=0.,max=20)
			model.set_param_hint('m2_totflux',min=0.,max=20)
			model.set_param_hint('m3_totflux',min=0.,max=20)
			model.set_param_hint('m4_totflux',min=0.,max=20)

			model.set_param_hint('m1_c',min=-0.5,max=0.5)

			model.set_param_hint('m1_x0',min=22.,max=28)
			model.set_param_hint('m2_x0',min=5.,max=26)
			model.set_param_hint('m3_x0',min=26.,max=38)
			model.set_param_hint('m4_x0',min=5.,max=26)

			params=model.make_params(m1_totflux=2,m1_c=0,m1_x0=26,m1_sigma=2,m2_totflux=0.2,m2_c=0,m2_x0=10,m2_sigma=2,m3_totflux=0.2,m3_c=0,m3_x0=30,m3_sigma=2,m4_totflux=0.2,m4_c=0,m4_x0=18,m4_sigma=2)

			params['m2_c'].vary=False
			params['m3_c'].vary=False
			params['m4_c'].vary=False


		result=model.fit(collapse_y,x=x_fit,params=params)

		print(result.fit_report())


		fig, (ax1, ax2) = pyplot.subplots(1, 2)
		ax1.imshow(subset_data,origin='lower')
		ax2.plot(x_fit,collapse_y)
		ax2.plot(x_fit_oversample,model.eval(result.params,x=x_fit_oversample))
		pyplot.savefig(SAVE_FOLDER+'opt_profile_fit_%s_mod%s.png'%(thisID,module))
		pyplot.clf()


		result.params['m1_c'].value=0.

		model_y=model.eval(result.params,x=x_fit)
		model_y=model_y/np.nansum(model_y)

		for j in range(len(opt_weight[0,:])):
			opt_weight[:,j]=model_y
		#pyfits.writeto('opt_weight_mod%s.fits'%module,opt_weight,overwrite=True)

		hdu.append(fits.ImageHDU(data=opt_weight,header=hd,name='OPT_WEIGHT%s'%module))


		#NOW EXTRACT THE 1D
	# Hornes optimal extraction of model
		data=hdu['CONT%s'%module].data
		sel_wav=(wav_array[200:980]>5008.24*(1+thisz)- 5008.24*(1+thisz)*10000/3E5)*(wav_array[200:980]<5008.24*(1+thisz)+ 5008.24*(1+thisz)*10000/3E5) #this

		#BACKGROUND SUBTRACT:
		standard=np.nanstd(data[:,200:980])
		med=np.nanmedian(data[:,200:980])
		print('standard:',standard,'med',med)
		use_dat=copy.deepcopy(data[:,200:980])
		use_dat[23:29,:]=np.nan

		sel_ignore=numpy.abs(use_dat)>5*standard
		use_dat[sel_ignore]=numpy.nan
		standard=numpy.nanstd(use_dat)
		med=numpy.nanmedian(use_dat[:,sel_wav])
		print('standard:',standard,'med',med)
		sel_ignore=numpy.abs(use_dat)>5*standard
		use_dat[sel_ignore]=numpy.nan
		standard=numpy.nanstd(use_dat)	
		med=numpy.nanmedian(use_dat[:,sel_wav])
		print('standard:',standard,'med',med)
		sel_ignore=numpy.abs(use_dat)>5*standard
		use_dat[sel_ignore]=numpy.nan
		standard=numpy.nanstd(use_dat)	
		med=numpy.nanmedian(use_dat[:,sel_wav])
		print('standard:',standard,'med',med)

		#data=data-med



		t1 = np.nansum(data*opt_weight,axis=0)
		t2 = np.nansum(opt_weight**2,axis=0)
		emline_extracted = t1/t2	

		ivar=errdata**-2
		
		t1 = np.nansum(ivar*opt_weight**2,axis=0)
		t2 = np.nansum(opt_weight,axis=0)
		emline_extracted_err = (t1/t2)**-0.5


		if rescale_noise==True:
			sel_wav=(wav_array>31500)*(wav_array<39500)
			use_err=copy.deepcopy(emline_extracted_err[sel_wav])
			sel_ignore=use_err==0.
			use_err[sel_ignore]=numpy.nan
			median_err=numpy.nanmedian(use_err)

			use_dat=copy.deepcopy(emline_extracted[sel_wav])
			standard=numpy.nanstd(use_dat)
			sel_ignore=numpy.abs(use_dat)>5*standard
			use_dat[sel_ignore]=numpy.nan
			standard=numpy.nanstd(use_dat)
			sel_ignore=numpy.abs(use_dat)>3*standard
			use_dat[sel_ignore]=numpy.nan			
			standard=numpy.nanstd(use_dat)		

			emline_extracted_err=emline_extracted_err * standard/median_err

			print(standard)

		COLS.append(fits.Column(name='flux_tot_%s'%module,unit='1E18 erg/s/cm2/A',format='E',array=emline_extracted))
		COLS.append(fits.Column(name='flux_tot_%s_err'%module,unit='1E18 erg/s/cm2/A',format='E',array=emline_extracted_err))
		
	cols=fits.ColDefs(COLS)#,col13,col14])
	hdu_1D = fits.BinTableHDU.from_columns(cols)
	#hdu=fits.PrimaryHDU(numpy.arange(100.))
	hdu_1D.writeto(SAVE_FOLDER+'spectrum_1D_%s.fits'%(thisID),overwrite=True)




	#hdu.writeto('test.fits',overwrite=True)


