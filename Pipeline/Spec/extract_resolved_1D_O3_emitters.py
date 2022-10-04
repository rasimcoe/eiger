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



FOLDER='/scratch/EIGER/identification/SPECTRA_12Sept/'

CATALOG='/scratch/EIGER/identification/J0100_photcat_v2_CONCAT_O3candidates_HYBRID_BACKWARD.fits'

cat=fits.open(CATALOG)

data=cat[1].data
IDlist=data.field('NUMBER')

z_guesslist=data.field('z_O3_estimate')

Nclumps_Y=np.zeros(len(IDlist))+1 ##The idea is to have this in the catalog to determine how many clumps *are separated in the Y direction*, this will help the fitting


#still under heavy development, I will break up the code in functions asap
IDlist=[18026]
z_guesslist=[6.76541752431055]
Nclumps_Y=[3]


for q in range(len(IDlist)):
	thisID=IDlist[q]
	thisz=z_guesslist[q]
	thisNclumps_Y=Nclumps_Y[q]

	hdu= fits.open(FOLDER+'v3corstacked_2D_%s.fits'%thisID)
	hd=hdu['EMLINE'].header
	data=hdu['EMLINE'].data

	data_A=hdu['EMLINEA'].data
	data_B=hdu['EMLINEB'].data

	data_A_orig=copy.deepcopy(data_A)

	pyfits.writeto('input_EMA.fits',data_A,overwrite=True)

	print(hd['CRVAL1'],hd['CDELT1'])

	ly,lx=np.shape(data)

	x_array=np.arange(0,lx,1)
	wav_array=x_array*hd['CDELT1'] + hd['CRVAL1']

	print(wav_array)
	#select 5008
	sel_wav=(wav_array>5008.24*(1+thisz)-200)*(wav_array<5008.24*(1+thisz)+200)

	use_image=data_A[:,sel_wav]
	ly2,lx2=np.shape(use_image)
	Yg, Xg = numpy.mgrid[:ly2, :lx2]



	if thisNclumps_Y==1:
		model=Model(rotated_2dgauss,independent_vars=('x','y'))
		params=model.make_params(A=0.1,x0=lx2/2,y0=26,sigma_x=2,sigma_y=2,theta=0)


	if thisNclumps_Y==3:
		model=Model(rotated_2dgauss,independent_vars=('x','y'),prefix='m1_') + Model(rotated_2dgauss,independent_vars=('x','y'),prefix='m2_') +  Model(rotated_2dgauss,independent_vars=('x','y'),prefix='m3_') 
		#print(model.param_names)
		model.set_param_hint('m1_sigma_x',min=0.5,max=3.5)
		model.set_param_hint('m2_sigma_x',min=0.5,max=3.5)
		model.set_param_hint('m3_sigma_x',min=0.5,max=3.5)

		model.set_param_hint('m1_sigma_y',min=0.5,max=4)
		model.set_param_hint('m2_sigma_y',min=0.5,max=4)
		model.set_param_hint('m3_sigma_y',min=0.5,max=4)

		params=model.make_params(m1_A=0.15,m1_x0=lx2/2,m1_y0=26,m1_sigma_x=2,m1_sigma_y=2,m1_theta=0, m2_A=0.1,m2_x0=lx2/2,m2_y0=29,m2_sigma_x=2,m2_sigma_y=2,m2_theta=0, m3_A=0.1,m3_x0=lx2/2,m3_y0=23,m3_sigma_x=2,m3_sigma_y=2,m3_theta=0)




	result=model.fit(use_image,x=Xg,y=Yg,params=params)

	model_image=model.eval(result.params,x=Xg,y=Yg)

	result_1=copy.deepcopy(result)
	result_1.params['m1_A'].value=0
	model23_image=model.eval(result_1.params,x=Xg,y=Yg)
	#print('Total flux',np.sum(model23_image))

	print(result.fit_report())


	fig, (ax1, ax2,ax3) = pyplot.subplots(1, 3)
	ax1.imshow(use_image)
	imgs = ax1.get_images()
	if len(imgs) > 0:
   		vmin, vmax = imgs[0].get_clim()

	ax1.imshow(use_image,vmin=0.5*vmin,vmax=0.5*vmax)

	ax2.imshow(model_image,vmin=0.5*vmin,vmax=0.5*vmax)
	ax3.imshow(use_image-model23_image,vmin=0.5*vmin,vmax=0.5*vmax)
	#pyplot.show()
	#pyplot.clf()



	data_A[:,sel_wav]=data_A[:,sel_wav]-model23_image



	#NOW APPLY THIS MODEL TO 4960
	sel_wav=(wav_array>4960.295*(1+thisz)-200)*(wav_array<4960.295*(1+thisz)+200)

	use_image=data_A[:,sel_wav]
	ly2,lx2=np.shape(use_image)
	Yg, Xg = numpy.mgrid[:ly2, :lx2]



	#model=Model(rotated_2dgauss,independent_vars=('x','y'))
	#params=model.make_params(A=0.1,x0=lx2/2,y0=26,sigma_x=2,sigma_y=2,theta=0)

	model=Model(rotated_2dgauss,independent_vars=('x','y'),prefix='m1_') + Model(rotated_2dgauss,independent_vars=('x','y'),prefix='m2_') +  Model(rotated_2dgauss,independent_vars=('x','y'),prefix='m3_') 

	params=result.params#model.make_params(m1_A=0.1,m1_x0=lx2/2,m1_y0=26,m1_sigma_x=2,m1_sigma_y=2,m1_theta=0, m2_A=0.1,m2_x0=lx2/2,m2_y0=29,m2_sigma_x=2,m2_sigma_y=2,m2_theta=0, m3_A=0.1,m3_x0=lx2/2,m3_y0=23,m3_sigma_x=2,m3_sigma_y=2,m3_theta=0)

	# model.set_param_hint('m1_A',0.3*params['m1_A'].value)
	# model.set_param_hint('m2_A',0.3*params['m2_A'].value)
	# model.set_param_hint('m3_A',0.3*params['m3_A'].value)

	params['m1_sigma_y'].vary=False
	params['m1_y0'].vary=False
	params['m1_theta'].vary=False

	params['m2_sigma_y'].vary=False
	params['m2_y0'].vary=False
	params['m2_theta'].vary=False

	params['m3_sigma_y'].vary=False
	params['m3_y0'].vary=False
	params['m3_theta'].vary=False

	model.set_param_hint('m1_sigma_x',min=0.98*params['m1_sigma_x'].value,max=1.02*params['m1_sigma_x'].value)
	model.set_param_hint('m2_sigma_x',min=0.98*params['m2_sigma_x'].value,max=1.02*params['m1_sigma_x'].value)
	model.set_param_hint('m3_sigma_x',min=0.98*params['m3_sigma_x'].value,max=1.02*params['m1_sigma_x'].value)

	model.set_param_hint('m1_sigma_y',min=0.98*params['m1_sigma_y'].value,max=1.02*params['m1_sigma_y'].value)
	model.set_param_hint('m2_sigma_y',min=0.98*params['m2_sigma_y'].value,max=1.02*params['m1_sigma_y'].value)
	model.set_param_hint('m3_sigma_y',min=0.98*params['m3_sigma_y'].value,max=1.02*params['m1_sigma_y'].value)

	model.set_param_hint('m1_x0',min=0.98*params['m1_x0'].value,max=1.02*params['m1_x0'].value)
	model.set_param_hint('m2_x0',min=0.98*params['m2_x0'].value,max=1.02*params['m1_x0'].value)
	model.set_param_hint('m3_x0',min=0.98*params['m3_x0'].value,max=1.02*params['m1_x0'].value)

	model.set_param_hint('m1_y0',min=0.98*params['m1_y0'].value,max=1.02*params['m1_y0'].value)
	model.set_param_hint('m2_y0',min=0.98*params['m2_y0'].value,max=1.02*params['m1_y0'].value)
	model.set_param_hint('m3_y0',min=0.98*params['m3_y0'].value,max=1.02*params['m1_y0'].value)

	result=model.fit(use_image,x=Xg,y=Yg,params=params)

	model_image=model.eval(result.params,x=Xg,y=Yg)

	result.params['m1_A'].value=0
	model23_image=model.eval(result.params,x=Xg,y=Yg)

	print(result.fit_report())


	fig, (ax1, ax2,ax3) = pyplot.subplots(1, 3)
	ax1.imshow(use_image)
	imgs = ax1.get_images()
	if len(imgs) > 0:
   		vmin, vmax = imgs[0].get_clim()

	ax1.imshow(use_image,vmin=0.5*vmin,vmax=0.5*vmax)

	ax2.imshow(model23_image,vmin=0.5*vmin,vmax=0.5*vmax)
	ax3.imshow(use_image-model_image,vmin=0.5*vmin,vmax=0.5*vmax)
	pyplot.show()


	data_A[:,sel_wav]=data_A[:,sel_wav]-model23_image

	pyfits.writeto('modified_EMA.fits',data_A,overwrite=True)

	pyfits.writeto('residual_EMA.fits',data_A_orig-data_A,overwrite=True)

	# pyplot.imshow(use_image)
	# pyplot.show()




