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

from astropy.io import fits
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

import copy
from PIL import Image

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

#
#OPEN THE RGB IMAGE -- THIS MUST HAVE THE SAME DIMENSIONS AS THE DIRECT IMAGE


for q in range(len(IDlist)):
	thisField=FIELDlist[q]
	thisz=Zlist[q]
	RGBFILE='/scratch/EIGER/identification/stiff_bin1_j%s.tif'%thisField[1:]
	directimage='/scratch/EIGER/identification/j%s_F356W.fits'%thisField[1:]
	img=Image.open(RGBFILE)
	img = np.asarray(img)
	img=img[::-1,:,:]


	hdu = fits.open(directimage)
	wcs = WCS(hdu['SCI'].header)

	thisID=IDlist[q]

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

	stamp=hdu['STAMP'].data ##F356W image only


	x_peak=int((6564.633*(1+thisz) -3E4)/9.75)


	start=100#int(thisx)-150
	end=-180#int(thisx)+60
	X_5008=150-0.5
	print(start,end)
#	#2D stacks
	image=data[10:-10,start:end]
	imageA=data_A[10:-10,start:end] 
	imageB=data_B[10:-10,start:end]
	scimageA=scidata_A[10:-10,start:end] 
	scimageB=scidata_B[10:-10,start:end]
	ly,lx=np.shape(image)


	#STAMP: F356W image
	ly2,lx2=np.shape(stamp)
	thumbsize=43
	stampcut=stamp[int(ly2/2)-thumbsize:int(ly2/2)+thumbsize,int(lx2/2)-thumbsize:int(lx2/2)+thumbsize]
	
	
	###ALTERNATIVE: USE RGB image
	thisx,thisy=wcs.all_world2pix(RAlist[q], DEClist[q], 0)#SkyCoord(RAlist[q],DEClist[q],unit="deg")
	#thisz=Zlist[q]
	print(thisx,thisy)
	thumbsize=50
	image2=img[int(thisy)-thumbsize:int(thisy)+thumbsize,int(thisx)-thumbsize:int(thisx)+thumbsize,:]
	#image2=copy.deepcopy(image)	
	
	
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



