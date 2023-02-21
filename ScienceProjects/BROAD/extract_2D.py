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



FOLDER='/scratch/EIGER/BROAD/SPECTRA_COLSEL/' #FOLDER WITH SPECTRA
SAVEFOLDER='/scratch/EIGER/BROAD/STAMP/'


CATALOG='/scratch/EIGER/BROAD/BROADsel_allfields_17022023_zguess.fits'

FIELDNAMES=['1120','0100','1148','0148']



with fits.open(CATALOG) as hdul:
    orig_table = hdul[1].data
    orig_cols = orig_table.columns
    

IDlist=orig_table.field('NUMBER')
RAlist=orig_table.field('ALPHA_J2000_det')
DEClist=orig_table.field('DELTA_J2000_det')
FIELDlist=orig_table.field('FIELD')
Zlist=orig_table.field('zguess')


for q in [4]:#range(len(IDlist)):

	thisField=FIELDlist[q]
	thisz=Zlist[q]

	thisID=IDlist[q]
	thisNclumps=3
	print('Now doing',thisField,thisID)
	#try:
	hdu= fits.open(FOLDER+'stacked_2D_%s_%s.fits'%(thisField,thisID))
	#except:
		#continue
	hd=hdu['EMLINE'].header
	data=hdu['EMLINE'].data

	data_A=hdu['EMLINEA'].data
	data_B=hdu['EMLINEB'].data
	scidata_A=hdu['SCIA'].data
	scidata_B=hdu['SCIB'].data

	data_A_orig=copy.deepcopy(data_A)

	x_peak=int((6564.633*(1+thisz) -3E4)/9.75)

	pyfits.writeto('input_EMA.fits',data_A,overwrite=True)

	print(hd['CRVAL1'],hd['CDELT1'])

	ly,lx=np.shape(data)

	x_array=np.arange(0,lx,1)
	wav_array=x_array*hd['CDELT1'] + hd['CRVAL1']

	print(wav_array)
	#select 5008
	sel_wav=(wav_array>6564.633*(1+thisz) - 300)*(wav_array< 6564.633*(1+thisz)   +300)


	use_image=data_A[:,sel_wav]
	ly2,lx2=np.shape(use_image)
	Yg, Xg = numpy.mgrid[:ly2, :lx2]



	if thisNclumps==1:
		model=Model(rotated_2dgauss,independent_vars=('x','y'))
		params=model.make_params(A=0.1,x0=lx2/2,y0=26,sigma_x=2,sigma_y=2,theta=0)


	if thisNclumps==3:
		model=Model(rotated_2dgauss,independent_vars=('x','y'),prefix='m1_') + Model(rotated_2dgauss,independent_vars=('x','y'),prefix='m2_') +  Model(rotated_2dgauss,independent_vars=('x','y'),prefix='m3_') 
		#print(model.param_names)
		model.set_param_hint('m1_sigma_x',min=0.3,max=4.5)
		model.set_param_hint('m2_sigma_x',min=0.3,max=7.5)
		model.set_param_hint('m3_sigma_x',min=0.3,max=10.5)

		model.set_param_hint('m1_sigma_y',min=0.2,max=4)
		model.set_param_hint('m2_sigma_y',min=0.2,max=2.5)
		model.set_param_hint('m3_sigma_y',min=0.2,max=1.5)

		model.set_param_hint('m3_theta',min=-10,max=10)		

		params=model.make_params(m1_A=0.0,m1_x0=lx2/2,m1_y0=20,m1_sigma_x=0.5,m1_sigma_y=0.5,m1_theta=0, m2_A=0.1,m2_x0=lx2/2,m2_y0=25,m2_sigma_x=2,m2_sigma_y=2,m2_theta=0, m3_A=0.1,m3_x0=lx2/2,m3_y0=25,m3_sigma_x=5,m3_sigma_y=2,m3_theta=0)


		params['m3_theta'].vary=False
		params['m1_theta'].vary=False
	result=model.fit(use_image,x=Xg,y=Yg,params=params)

	model_image=model.eval(result.params,x=Xg,y=Yg)

	##In case multiple objects
	#result_1=copy.deepcopy(result)
	#result_1.params['m1_A'].value=0
	#model23_image=model.eval(result_1.params,x=Xg,y=Yg)
	#print('Total flux',np.sum(model23_image))

	print(result.fit_report())


	fig, (ax1, ax2,ax3) = pyplot.subplots(1, 3)
	ax1.imshow(use_image)
	imgs = ax1.get_images()
	if len(imgs) > 0:
   		vmin, vmax = imgs[0].get_clim()

	ax1.imshow(use_image,vmin=0.3*vmin,vmax=0.3*vmax)

	ax2.imshow(model_image,vmin=0.3*vmin,vmax=0.3*vmax)
	ax3.imshow(use_image-model_image,vmin=0.3*vmin,vmax=0.3*vmax)
	pyplot.show()
	pyplot.clf()



	data_A[:,sel_wav]=data_A[:,sel_wav]-model_image

	pyfits.writeto('modified_EMA.fits',data_A,overwrite=True)

	pyfits.writeto('residual_EMA.fits',data_A_orig-data_A,overwrite=True)

	# pyplot.imshow(use_image)
	# pyplot.show()




