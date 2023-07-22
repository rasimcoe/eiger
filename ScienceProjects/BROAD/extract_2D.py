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
from matplotlib import colors
import matplotlib
import numpy as np

from astropy.io import fits

import matplotlib

matplotlib.rcParams.update({ 'font.family': 'serif','mathtext.fontset':'cm'}) # change for diff style

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

def rotated_2dgauss_totflux(x, y, A, x0, y0, sigma_x, sigma_y, theta):
    theta = np.radians(theta)
    sigx2 = sigma_x**2; sigy2 = sigma_y**2
    a = np.cos(theta)**2/(2*sigx2) + np.sin(theta)**2/(2*sigy2)
    b = np.sin(theta)**2/(2*sigx2) + np.cos(theta)**2/(2*sigy2)
    c = np.sin(2*theta)/(4*sigx2) - np.sin(2*theta)/(4*sigy2)
    
    expo = -a*(x-x0)**2 - b*(y-y0)**2 - 2*c*(x-x0)*(y-y0)
    return A*(sigma_x)**-0.5 * (sigma_y)**-0.5 * (2*np.pi)**-0.5 *np.exp(expo)    


def read_fitsdata(hdu,module='A',extension='EMLINE'):
	thisdata=hdu[extension+module].data
	thisdata_err=hdu['ERR%s'%module].data**0.5  #THIS IS BECAUSE THE ERR extension is actually VARIANCE -- needs to be fixed

	#RENORMALISE NOISE
	use_err=copy.deepcopy(thisdata_err[:,200:980])
	sel_ignore=use_err==0.
	use_err[sel_ignore]=numpy.nan

	use_dat=copy.deepcopy(thisdata[:,200:980])
	standard=numpy.nanstd(use_dat)
	sel_ignore=numpy.abs(use_dat)>5*standard
	use_dat[sel_ignore]=numpy.nan
	standard=numpy.nanstd(use_dat)
	sel_ignore=numpy.abs(use_dat)>3*standard
	use_dat[sel_ignore]=numpy.nan
	standard=numpy.nanstd(use_dat-numpy.nanmedian(use_dat))		
	sel_ignore=numpy.abs(use_dat)>3*standard
	use_dat[sel_ignore]=numpy.nan
	standard=numpy.nanstd(use_dat-numpy.nanmedian(use_dat))	

	sel_ignore=numpy.abs(use_dat)>3*standard
	use_dat[sel_ignore]=numpy.nan
	standard=numpy.nanstd(use_dat-numpy.nanmedian(use_dat))	
	
	errdata=thisdata_err * standard/numpy.nanmedian(use_err) #Renormalising
	
	unmasked_data=copy.deepcopy(thisdata)
	unmasked_data[:,200:980]=use_dat

	
	return thisdata,errdata,unmasked_data
	

def read_fitsdata_SCI(hdu,module='A',extension='EMLINE',thisz=5.0):
	thisdata=hdu[extension+module].data
	hd=hdu[extension+module].header
	print(hd['CRVAL1'],hd['CDELT1'])


	
	if extension=='SCI': #kernel filter!
		ly,lx=np.shape(thisdata)
		n_y=ly
		n_x=lx
		continua_img = np.zeros((n_y, n_x))

		x_array=np.arange(0,lx,1)
		wav_array=x_array*hd['CDELT1'] + hd['CRVAL1']

		print(wav_array)
		sel_wav=(wav_array>6564.633*(1+thisz) - 260)*(wav_array< 6564.633*(1+thisz)   +260)
		orig_data=copy.deepcopy(thisdata)
		
		#mask line!
		thisdata[:,sel_wav]=numpy.nan
	
		#Mid kernel
		kx=141
		kx_gap=31
		
		#Widekernel
		kx=251
		kx_gap=71		
	#
		#kx=55
		#kx_gap=13
		ky=1
		kernel=str(kx)+'x'+str(ky)
		print('Kernel size: ', kernel, flush=True)

		idx_x = np.fromfunction(lambda  i, j: i+j, (n_x, kx), dtype=np.int64) - kx // 2
		idx_y = np.fromfunction(lambda  i, j: i+j, (n_y, ky), dtype=np.int64) - ky // 2

	# Exclude the center gap
		if kx_gap>0:
	    		print('Gap size: ', kx_gap, flush=True)
	    		cols_remain=np.append(np.arange((kx - kx_gap)/2., dtype=np.int64),
	                           np.flip(kx-1-np.arange((kx - kx_gap)/2., dtype=np.int64)))
	    		idx_x = idx_x[:, cols_remain]

		idx_x[idx_x < 0]=0
		idx_x[idx_x > n_x-1]=n_x-1
		idx_y[idx_y < 0]=0
		idx_y[idx_y > n_y-1]=n_y-1

		print( ': Loop start (',flush=True)
		for iy in np.arange(n_y):
	    	#print('Loop: ', iy, ' /', n_y)
			med_tmp = np.nanmedian(thisdata[idx_y[iy,0]:idx_y[iy,-1]+1,idx_x], axis=[0,2])
			continua_img[iy,:] = med_tmp[:]
		print( ': Loop end (', flush=True)
		fits.writeto('cont.fits',continua_img,overwrite=True)
		
		
		thisdata=orig_data-continua_img
	
	
	thisdata_err=hdu['ERR%s'%module].data**0.5  #THIS IS BECAUSE THE ERR extension is actually VARIANCE -- needs to be fixed

	#RENORMALISE NOISE
	use_err=copy.deepcopy(thisdata_err[:,200:980])
	sel_ignore=use_err==0.
	use_err[sel_ignore]=numpy.nan

	use_dat=copy.deepcopy(thisdata[:,200:980])
	standard=numpy.nanstd(use_dat)
	sel_ignore=numpy.abs(use_dat)>5*standard
	use_dat[sel_ignore]=numpy.nan
	standard=numpy.nanstd(use_dat)
	sel_ignore=numpy.abs(use_dat)>3*standard
	use_dat[sel_ignore]=numpy.nan
	standard=numpy.nanstd(use_dat-numpy.nanmedian(use_dat))		
	sel_ignore=numpy.abs(use_dat)>3*standard
	use_dat[sel_ignore]=numpy.nan
	standard=numpy.nanstd(use_dat-numpy.nanmedian(use_dat))	

	sel_ignore=numpy.abs(use_dat)>3*standard
	use_dat[sel_ignore]=numpy.nan
	standard=numpy.nanstd(use_dat-numpy.nanmedian(use_dat))	
	
	errdata=thisdata_err * standard/numpy.nanmedian(use_err) #Renormalising
	
	unmasked_data=copy.deepcopy(thisdata)
	unmasked_data[:,200:980]=use_dat

	
	return thisdata,errdata,unmasked_data
	



FOLDER='/scratch/EIGER/BROAD/SPECTRA_COLSEL/' #FOLDER WITH SPECTRA
SAVE_FOLDER='/scratch/EIGER/BROAD/EXTRACTIONS/'


CATALOG='/scratch/EIGER/BROAD/BROADsel_allfields_17022023_zguess.fits'
CATALOG='/scratch/EIGER/BROAD/BROADsel_allfields_23022023_zguess_classed.fits'
FIELDNAMES=['1120','0100','1148','0148']



with fits.open(CATALOG) as hdul:
    orig_table = hdul[1].data
    orig_cols = orig_table.columns
    

IDlist=orig_table.field('NUMBER')
RAlist=orig_table.field('ALPHA_J2000_det')
DEClist=orig_table.field('DELTA_J2000_det')
FIELDlist=orig_table.field('FIELD')
Zlist=orig_table.field('zguess')
Nlist=orig_table.field('Ncomponents')
MODlist=orig_table.field('MODULE')


for q in [15]:#range(len(IDlist)):

	thisField=FIELDlist[q]
	thisz=Zlist[q]

	thisID=IDlist[q]
	thisNclumps=3.#Nlist[q]
	thisMod=MODlist[q]
	if thisMod=='BOTH':
		thisMod=''
	print('Now doing',thisField,thisID)
	#try:
	hdu= fits.open(FOLDER+'stacked_2D_%s_%s.fits'%(thisField,thisID))
	thisID=str(thisID)+'_widekernel'
	#except:
		#continue
	hd=hdu['EMLINE'].header
	data=hdu['EMLINE'].data

	#data_A,errdata_A,unmasked_data_A=read_fitsdata(hdu,module=thisMod,extension='EMLINE')
	data_A,errdata_A,unmasked_data_A=read_fitsdata_SCI(hdu,module=thisMod,extension='SCI',thisz=thisz)
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
	err_image=errdata_A[:,sel_wav]
	unmasked_image=unmasked_data_A[:,sel_wav]
	use_image=use_image-numpy.nanmedian(unmasked_image)
	
	ly2,lx2=np.shape(use_image)
	Yg, Xg = numpy.mgrid[:ly2, :lx2]

	print(ly2,lx2)

	if thisNclumps==1:
		model=Model(rotated_2dgauss,independent_vars=('x','y'))
		model.set_param_hint('sigma_x',min=1.0,max=10.5)
		model.set_param_hint('sigma_y',min=0.5,max=2.5)
		model.set_param_hint('x0',min=20.,max=40.)
		model.set_param_hint('y0',min=23.,max=27.)						
		params=model.make_params(A=0.1,x0=lx2/2,y0=26,sigma_x=2,sigma_y=2,theta=0)

	if thisNclumps==2:
		model=Model(rotated_2dgauss,independent_vars=('x','y'),prefix='m1_') + Model(rotated_2dgauss,independent_vars=('x','y'),prefix='m2_') 
		#print(model.param_names)
		model.set_param_hint('m1_sigma_x',min=1.0,max=2.5)
		model.set_param_hint('m2_sigma_x',min=1.1,max=13.5)
		model.set_param_hint('m1_A',min=0.0,max=2.)
		model.set_param_hint('m2_A',min=0.0,max=2.)

		model.set_param_hint('m1_x0',min=20.,max=40.)
		model.set_param_hint('m2_x0',min=30.,max=40.)
	
		model.set_param_hint('m1_y0',min=15.,max=40.)
		model.set_param_hint('m2_y0',min=23.,max=27.)
							
		model.set_param_hint('m1_sigma_y',min=0.5,max=2.5)
		model.set_param_hint('m2_sigma_y',min=0.5,max=2.5)

		model.set_param_hint('m1_theta',min=0,max=360)		
		model.set_param_hint('m2_theta',min=-10,max=10)		

		params=model.make_params(m1_A=0.05,m1_x0=lx2/2,m1_y0=25,m1_sigma_x=2.,m1_sigma_y=1.0,m1_theta=20, m2_A=0.1,m2_x0=lx2/2,m2_y0=25,m2_sigma_x=8,m2_sigma_y=1.0,m2_theta=0)
	

		params['m2_theta'].vary=False


	if thisNclumps==3:
		model=Model(rotated_2dgauss,independent_vars=('x','y'),prefix='m1_') + Model(rotated_2dgauss,independent_vars=('x','y'),prefix='m2_') +  Model(rotated_2dgauss,independent_vars=('x','y'),prefix='m3_') 
		#print(model.param_names)
		model.set_param_hint('m1_sigma_x',min=1.1,max=2.5)
		model.set_param_hint('m2_sigma_x',min=2.5,max=13.5)
		model.set_param_hint('m3_sigma_x',min=1.1,max=5.5)

		model.set_param_hint('m1_sigma_y',min=0.5,max=4.)
		model.set_param_hint('m2_sigma_y',min=0.5,max=2.5)
		model.set_param_hint('m3_sigma_y',min=0.5,max=4.5)

		model.set_param_hint('m1_A',min=0.0,max=2.)
		model.set_param_hint('m2_A',min=0.0,max=2.)
		model.set_param_hint('m3_A',min=0.0,max=2.)
		
		model.set_param_hint('m1_x0',min=20.,max=40.)
		model.set_param_hint('m2_x0',min=20.,max=40.)
		model.set_param_hint('m3_x0',min=20.,max=40.)	
		
		model.set_param_hint('m1_y0',min=20.,max=35.)
		model.set_param_hint('m2_y0',min=23.,max=27.)
		model.set_param_hint('m3_y0',min=15.,max=40.)					
					
		model.set_param_hint('m1_theta',min=0,max=360)		
		model.set_param_hint('m2_theta',min=-10,max=10)		
		model.set_param_hint('m3_theta',min=0,max=360)		

		params=model.make_params(m1_A=0.0,m1_x0=lx2/2,m1_y0=25,m1_sigma_x=0.5,m1_sigma_y=2.,m1_theta=0, m2_A=0.1,m2_x0=lx2/2,m2_y0=25,m2_sigma_x=7.,m2_sigma_y=1.,m2_theta=0, m3_A=0.1,m3_x0=24.7,m3_y0=20,m3_sigma_x=1,m3_sigma_y=2,m3_theta=0)


		params['m3_theta'].vary=False
		params['m2_theta'].vary=False
		#params['m1_theta'].vary=False
	result=model.fit(use_image,x=Xg,y=Yg,params=params,weights=err_image**-1,method='least_squares')

	model_image=model.eval(result.params,x=Xg,y=Yg)

	fwhm_x=2.355*result.params['m2_sigma_x'].value
	fwhm_vel=3E5 * fwhm_x * 9.75 / (6564.633*(1+thisz))

	fwhm_x_narrow=2.355*result.params['m1_sigma_x'].value
	fwhm_vel_narrow=3E5 * fwhm_x_narrow * 9.75 / (6564.633*(1+thisz))

	##In case multiple objects
	#result_1=copy.deepcopy(result)
	#result_1.params['m1_A'].value=0
	#model23_image=model.eval(result_1.params,x=Xg,y=Yg)
	#print('Total flux',np.sum(model23_image))

	if thisNclumps==2:
		result_1=copy.deepcopy(result)
		result_2=copy.deepcopy(result)		
		result_1.params['m1_A'].value=0
		result_2.params['m2_A'].value=0

		model_broad=model.eval(result_1.params,x=Xg,y=Yg)	
		model_narrow=model.eval(result_2.params,x=Xg,y=Yg)	
		model_third=model_narrow-model_narrow

	if thisNclumps==3:
		result_1=copy.deepcopy(result)
		result_2=copy.deepcopy(result)
		result_3=copy.deepcopy(result)				
		result_1.params['m1_A'].value=0
		result_1.params['m3_A'].value=0
				
		result_2.params['m2_A'].value=0
		result_2.params['m3_A'].value=0
		
		result_3.params['m1_A'].value=0
		result_3.params['m2_A'].value=0
				
		model_broad=model.eval(result_1.params,x=Xg,y=Yg)	
		model_narrow=model.eval(result_2.params,x=Xg,y=Yg)
		model_third=model.eval(result_3.params,x=Xg,y=Yg)

	print(result.fit_report())
	print('Total narrow,broad,third',9.75*numpy.nansum(model_narrow),9.75*numpy.nansum(model_broad),9.75*numpy.nansum(model_third))				
	#STOP



	#SAVE DIAGNOSTIC FIGURE
	fig, (ax1, ax2,ax3,ax4,ax5,ax6) = pyplot.subplots(6,1,figsize=(1.3,6.))
	ax1.imshow(use_image)
	imgs = ax1.get_images()
	fig.suptitle(thisField+'-'+str(IDlist[q]))
	if len(imgs) > 0:
   		vmin, vmax = imgs[0].get_clim()
   		MIN,MAX=vmin,vmax
	vlmn=0.05
	vlmx=0.15
	ax1.imshow(5*use_image,cmap='cubehelix',norm=colors.PowerNorm(gamma=0.33,vmin=-0.005,vmax=0.6),origin='lower') #DATA

	ax2.imshow(5*model_narrow,cmap='cubehelix',norm=colors.PowerNorm(gamma=0.33,vmin=-0.005,vmax=0.6),origin='lower') #NARROW
	ax3.imshow(5*model_broad,cmap='cubehelix',norm=colors.PowerNorm(gamma=0.33,vmin=-0.005,vmax=0.6),origin='lower') #BROAD
	ax4.imshow(5*model_third,cmap='cubehelix',norm=colors.PowerNorm(gamma=0.33,vmin=-0.005,vmax=0.6),origin='lower') #NARROW
	ax5.imshow(5*model_image,cmap='cubehelix',norm=colors.PowerNorm(gamma=0.33,vmin=-0.005,vmax=0.6),origin='lower') #MODE	
	ax6.imshow(5*(use_image-model_image),cmap='cubehelix',norm=colors.PowerNorm(gamma=0.33,vmin=-0.005,vmax=0.6),origin='lower') #RESIDUAL



	print(np.shape(use_image))
	ax2.plot([1,60],[25.5,25.5],color='white',ls='--',lw=1,alpha=0.6)
	ax3.plot([1,60],[25.5,25.5],color='white',ls='--',lw=1,alpha=0.6)
	ax4.plot([1,60],[25.5,25.5],color='white',ls='--',lw=1,alpha=0.6)
	ax5.plot([1,60],[25.5,25.5],color='white',ls='--',lw=1,alpha=0.6)	
	
	
	
	ax1.text(5,40,'DATA',color='white')
	ax2.text(5,40,'NARROW',color='white')
	ax3.text(5,40,'BROAD',color='white')
	ax4.text(5,40,'THIRD',color='white')
	ax5.text(5,40,'MODEL',color='white')
	ax6.text(5,40,'RESIDUAL',color='white')

	#ax1.text(5,10,'z ~ '+str(thisz),color='white')
	#ax2.text(5,10,str(int(fwhm_vel_narrow))+' km/s',color='white')	
	#ax3.text(5,10,str(int(fwhm_vel))+' km/s',color='white')

	ax6.minorticks_on()
	ax1.tick_params(left=False,right=False,bottom=False,top=False,labelleft=False,labelbottom=False) # Get ticks to look nice
	ax2.tick_params(left=False,right=False,bottom=False,top=False,labelleft=False,labelbottom=False) # Get ticks to look nice
	ax3.tick_params(left=False,right=False,bottom=False,top=False,labelleft=False,labelbottom=False) # Get ticks to look nice
	ax4.tick_params(left=False,right=False,bottom=False,top=False,labelleft=False,labelbottom=False) # Get ticks to look nice
	ax5.tick_params(left=False,right=False,bottom=False,top=False,labelleft=False,labelbottom=False) # Get ticks to look nice
	ax6.tick_params(which='both',left=False,right=False,top=False,labelleft=False) # Get ticks to look nice
	
	ax6.set_xticks((30-2*13,30.5,30.5+2*13),(-2000,0,2000))
	ax6.set_xlabel(r'$\Delta v$ [km s$^{-1}$]')

	pyplot.tight_layout()
	#pyplot.show()
	#stop
	pyplot.savefig(SAVE_FOLDER+'fits_%s_%s.pdf'%(thisField,thisID),dpi=150)
	pyplot.clf()
	STOP



	#NOW SAVE "CLEANED" 2D
	 
	data_A[:,sel_wav]=data_A[:,sel_wav]-model_third

	fits.writeto(SAVE_FOLDER+'cleaned_2d_for_extraction_%s_%s.fits'%(thisField,thisID),data_A,header=hd,overwrite=True)
	fits.writeto(SAVE_FOLDER+'cleaned_2d_for_extraction_%s_%s_err.fits'%(thisField,thisID),errdata_A,header=hd,overwrite=True)



	#NOW EXTRACT THE 1D SPECTRUM:
	#CREATE OPTIMAL WEIGHT
	model_y=model_broad+model_narrow

	optmodel=numpy.nansum(model_y,axis=1)
	optmodel=optmodel/np.nansum(optmodel)
	COLS=[]

	opt_weight=numpy.zeros(numpy.shape(data_A))
	for j in range(len(data_A[0,:])):
		opt_weight[:,j]=optmodel
		#pyfits.writeto('opt_weight_mod%s.fits'%module,opt_weight,overwrite=True)


	# Hornes optimal extraction of model
	t1 = np.nansum(data_A*opt_weight,axis=0)
	t2 = np.nansum(opt_weight**2,axis=0)
	emline_extracted = t1/t2	

	ivar=errdata_A**-2
		
	t1 = np.nansum(ivar*opt_weight**2,axis=0)
	t2 = np.nansum(opt_weight,axis=0)
	emline_extracted_err = (t1/t2)**-0.5
	
	COLS.append(fits.Column(name='wavelength',unit='angstrom',format='E',array=wav_array))
	COLS.append(fits.Column(name='flux',unit='1E18 erg/s/cm2/A',format='E',array=emline_extracted))
	COLS.append(fits.Column(name='flux_err',unit='1E18 erg/s/cm2/A',format='E',array=emline_extracted_err))
		
	cols=fits.ColDefs(COLS)#,col13,col14])
	hdu_1D = fits.BinTableHDU.from_columns(cols)
	hdu=fits.PrimaryHDU(numpy.arange(100.))
	hdu_1D.writeto(SAVE_FOLDER+'cleaned_spectrum_1D_%s_%s.fits'%(thisField,thisID),overwrite=True)






