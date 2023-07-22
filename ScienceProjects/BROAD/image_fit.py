from astropy.io import fits as pyfits #no need to call it pyfits, sorry about that
from matplotlib import pyplot
import numpy
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
from scipy.interpolate import interp1d
from astropy.convolution import Gaussian2DKernel
import scipy.ndimage as snd
from astropy.convolution import convolve
import os
from astropy import units as u
from photutils.psf.matching import resize_psf

import matplotlib
pyplot.rcParams['xtick.labelsize']=16
pyplot.rcParams['ytick.labelsize']=16
pyplot.rcParams['axes.labelsize']=16
pyplot.rcParams['image.origin']='lower'
from matplotlib.colors import LogNorm
matplotlib.rcParams['axes.linewidth'] = 1.5
matplotlib.rcParams['xtick.major.size'] = 6
matplotlib.rcParams['xtick.minor.size'] = 4
matplotlib.rcParams['xtick.major.width'] = 2.
matplotlib.rcParams['xtick.minor.width'] = 1.5
matplotlib.rcParams['ytick.major.size'] = 6
matplotlib.rcParams['ytick.minor.size'] = 4
matplotlib.rcParams['ytick.major.width'] = 2.
matplotlib.rcParams['ytick.minor.width'] = 1.5
matplotlib.rcParams['image.origin']='lower'
matplotlib.rcParams.update({'font.size': 16, 'font.family': 'serif','mathtext.fontset':'cm'}) # change for diff style

import numpy as np
from scipy.ndimage import rotate

from astropy.io import fits
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

import copy
from PIL import Image

import pyimfit
import webbpsf #export WEBBPSF_PATH=/scratch/EIGER/webbpsf/webbpsf-data
def get_model_JWST(FILTER, oversamp_fact=3, fov_arcsec=2, rot_angle=0.89):
    #get F356W PSF with extra fact 2 in oversamp
    instrument = webbpsf.NIRCam()
    instrument.filter=FILTER
    instrument.options['parity'] = 'odd'
    instrument.pixelscale = 0.03

    psf_model = instrument.calc_psf(oversample=oversamp_fact, fov_arcsec=fov_arcsec)

    psf_data_rot = rotate(psf_model['OVERSAMP'].data, rot_angle*u.deg, reshape=False)
    
    pix_clip = np.ceil(psf_model['OVERSAMP'].shape[0] * np.tan(rot_angle*u.deg)).astype(int)
    psf_data_rot_crop = psf_data_rot[pix_clip:-pix_clip,pix_clip:-pix_clip]

    if oversamp_fact == 1:
        return psf_data_rot_crop
    else:
        #trim array so final array is odd after resizing
        factor = np.floor(psf_data_rot_crop.shape[0]/oversamp_fact)
        trim = int((psf_data_rot_crop.shape[0]-(factor-(1 - (factor % 2)))*oversamp_fact)/2)
        psf_data_rot_crop_match = psf_data_rot_crop[trim:psf_data_rot_crop.shape[0]-trim,trim:psf_data_rot_crop.shape[1]-trim]
        det_samp = resize_psf(psf_data_rot_crop_match, 0.03/oversamp_fact, 0.03)
        return det_samp


def qso_galaxy_model(x0_qso, y0_qso, x0_host,y0_host, I_tot, PA_disk, ell_disk, I_0, h):
    model = pyimfit.ModelDescription()
    # define the limits on X0 and Y0 as +/-10 pixels relative to initial values

    qso = pyimfit.make_imfit_function('PointSource', label='qso')
    qso.I_tot.setValue(I_tot, [0, 10])
    qso.X0.setValue(x0_qso, [23,27])
    qso.Y0.setValue(y0_qso, [23,27])


    disk = pyimfit.make_imfit_function('Exponential', label='disk')
    disk.PA.setValue(PA_disk, [0, 180])
    disk.ell.setValue(ell_disk, [0, 1])
    disk.I_0.setValue(I_0, [0, 10*I_0])
    disk.h.setValue(h, [0, 10*h])
    disk.X0.setValue(x0_host, [20,30])
    disk.Y0.setValue(y0_host, [20,30])

    model.addFunction(qso)
    model.addFunction(disk)
    
    print('---->',model.params)

    return model



FOLDER='/scratch/EIGER/BROAD/SPECTRA_COLSEL/' #FOLDER WITH SPECTRA
SAVEFOLDER='/scratch/EIGER/BROAD/STAMP/'


CATALOG='/scratch/EIGER/BROAD/BROADsel_allfields_23022023_zguess_classed_withfits.fits'

FIELDNAMES=['1120','0100','1148','0148']

configfile='point_exp.config'

with fits.open(CATALOG) as hdul:
    orig_table = hdul[1].data
    orig_cols = orig_table.columns
    

IDlist=orig_table.field('NUMBER')
RAlist=orig_table.field('ALPHA_J2000_det')
DEClist=orig_table.field('DELTA_J2000_det')
FIELDlist=orig_table.field('FIELD')
Zlist=orig_table.field('zguess')
CONF=orig_table.field('HIGH_CONFID')

#



for q in range(len(IDlist)):
	thisField=FIELDlist[q]
	thisz=Zlist[q]
	if CONF[q] == False:
		continue
	RGBFILE='/scratch/EIGER/identification/stiff_bin1_j%s.tif'%thisField[1:]
	directimage='/scratch/EIGER/identification/j%s_F115W.fits'%thisField[1:]
	directimage_wht='/scratch/EIGER/identification/j%s_F115W_wht.fits'%thisField[1:]


	hdu = fits.open(directimage)
	wcs = WCS(hdu['SCI'].header)
	data_F356W=hdu['SCI'].data
	
	weight=fits.getdata(directimage_wht)

	thisID=IDlist[q]
	print('Now doing',thisField,thisID)


	###NOW LOAD
	#THEN DO PYIMFIT
	
	###ALTERNATIVE: USE RGB image
	thisx,thisy=wcs.all_world2pix(RAlist[q], DEClist[q], 0)#SkyCoord(RAlist[q],DEClist[q],unit="deg")
	#thisz=Zlist[q]
	print(thisx,thisy)
	thumbsize=25
	image2=data_F356W[int(thisy)-thumbsize:int(thisy)+thumbsize,int(thisx)-thumbsize:int(thisx)+thumbsize]
	whtcut=weight[int(thisy)-thumbsize:int(thisy)+thumbsize,int(thisx)-thumbsize:int(thisx)+thumbsize]
	#image2=copy.deepcopy(image)	
	fits.writeto('/scratch/EIGER/BROAD/STAMP/F115W_image_%s_%s.fits'%(thisField,thisID),image2,overwrite=True)
	fits.writeto('/scratch/EIGER/BROAD/STAMP/F115W_image_%s_%s_wht.fits'%(thisField,thisID),1./np.sqrt(whtcut),overwrite=True)

STOP
for q in [0]:
	#/net/theia/scratch/mattheej/MUSE_VR7/imfit-1.6.1/imfit image.fits -c point_exp.config --psf PSF_F356W.fits --noise wht.fits  --save-model bestfit.fits --save-residual residual.fits

	STOP
	psf_model = get_model_JWST('F356W', oversamp_fact=1, fov_arcsec=1.5, rot_angle=0.89)
	fits.writeto('PSF_F356W.fits',psf_model,overwrite=True)
	#psf_model=fits.getdata('PSF_F200W.fits')
	
	
	model_desc = pyimfit.ModelDescription.load(configfile)
	
	
	fitter_psf = pyimfit.Imfit(model_desc, psf_model)
	fitter_psf.loadData(image2, error=1./np.sqrt(whtcut))#, mask=total_mask)
	result_psf = fitter_psf.doFit()
	

	#if imfit_fitter.fitConverged is True:
	#	print("Fit converged: chi^2 = {0}, reduced chi^2 = {1}".format(imfit_fitter.fitStatistic,imfit_fitter.reducedFitStatistic))
	#	print("Best-fit parameter values:")
	#	print(imfit_fitter.getRawParameters())
	#print('Results',fitter_psf.getRawParameters())
	print('Results',result_psf)
	print(fitter_psf.getModelFluxes())
	bestmodel= fitter_psf.getModelImage(newParameters=result_psf.params)
    	
    	
    	
	fig, (ax1, ax2,ax3,ax4,ax5,ax6) = pyplot.subplots(1, 6,figsize=(10,2.5))
	ax1.imshow(image2)
	imgs = ax1.get_images()
	fig.suptitle(thisField+'-'+str(thisID))
	if len(imgs) > 0:
   		vmin, vmax = imgs[0].get_clim()

	ax1.imshow(image2,vmin=0.01*vmin,vmax=0.05*vmax,origin='lower') #DATA

	ax2.imshow(bestmodel,vmin=0.01*vmin,vmax=0.05*vmax,origin='lower') #NARROW
	ax3.imshow(bestmodel,vmin=0.01*vmin,vmax=0.05*vmax,origin='lower') #BROAD
	ax4.imshow(bestmodel,vmin=0.01*vmin,vmax=0.05*vmax,origin='lower') #NARROW
	ax5.imshow(bestmodel,vmin=0.01*vmin,vmax=0.05*vmax,origin='lower') #MODE	
	ax6.imshow(image2-bestmodel,vmin=0.01*vmin,vmax=0.05*vmax,origin='lower') #RESIDUAL
	pyplot.show()


	
	STOP
	
	
	
	
	fig=pyplot.figure(figsize=(5, 5))


	axstamp = pyplot.axes([0.14,0.13,0.84,0.84])

	#F356W stamp:
	#im = axstamp.imshow(stampcut,cmap='viridis',vmin=-0.007,vmax=0.046,origin='lower',aspect='equal',interpolation='none')


	#RGB stamp:
	im = axstamp.imshow(image2,origin='lower',aspect='equal')
	axstamp.plot([thumbsize+1,thumbsize+1],[0,thumbsize-12],lw=2,color='white',ls='--')
	axstamp.plot([0,thumbsize-12],[thumbsize+1,thumbsize+1],lw=2,color='white',ls='--')
	axstamp.set_xticks((50-33,50,50+33),(-1,0,1))
	axstamp.set_yticks((50-33,50,50+33),(-1,0,1))
	axstamp.set_xlabel(r'$\Delta$x [arcsec]')
	axstamp.set_ylabel(r'$\Delta$y [arcsec]')
	
	axstamp.text(8,86,thisField+'-'+str(thisID),color='white',fontsize=20)	
	



	pyplot.tight_layout()
	#pyplot.savefig('test.png',dpi=120)
	pyplot.savefig(SAVEFOLDER+'stamp_%s_%s.png'%(thisField,thisID),dpi=120)
	pyplot.clf()



