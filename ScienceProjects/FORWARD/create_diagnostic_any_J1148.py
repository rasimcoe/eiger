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


CATALOG='/scratch/EIGER/identification/J1148+5251_photcat_v1_noisemodel_short_withCandidateDoublets.fits' 

CATALOG='/scratch/EIGER/identification/J1148_REDcrit_mag26.fits'

FOLDER='/scratch/EIGER/identification/SPECTRA_J1148/' #FOLDER WITH SPECTRA

field='J1148'
with fits.open(CATALOG) as hdul:
    orig_table = hdul[1].data
    orig_cols = orig_table.columns
    

IDlist=orig_table.field('NUMBER')
RAlist=orig_table.field('ALPHA_J2000_det')
DEClist=orig_table.field('DELTA_J2000_det')


#
#OPEN THE RGB IMAGE -- THIS MUST HAVE THE SAME DIMENSIONS AS THE DIRECT IMAGE
img=Image.open('/scratch/EIGER/identification/stiff_bin1_j1148.tif')
img = np.asarray(img)
img=img[::-1,:,:]

directimage='/scratch/EIGER/identification/j1148_r.fits' #this is F356W purely SCI extension 
hdu = fits.open(directimage)
wcs = WCS(hdu['SCI'].header)

for q in range(len(IDlist)):

	thisID=IDlist[q]

	try:
		hdu= fits.open(FOLDER+'stacked_2D_J1148_%s.fits'%thisID)
	except:
		continue
	hd=hdu['EMLINE'].header
	data=hdu['EMLINE'].data

	data_A=hdu['EMLINEA'].data
	data_B=hdu['EMLINEB'].data
	scidata_A=hdu['SCIA'].data
	scidata_B=hdu['SCIB'].data

	stamp=hdu['STAMP'].data ##F356W image only


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
	
	
	
	
	
	
	
	fig=pyplot.figure(figsize=(28.9, 5.7))

	ax = pyplot.axes([0.01,0.04,0.7,0.18])
	ax2 = pyplot.axes([0.01,0.24,0.7,0.18])
	ax3 = pyplot.axes([0.01,0.44,0.7,0.18])
	ax4 = pyplot.axes([0.01,0.64,0.7,0.18])
	ax5 = pyplot.axes([0.01,0.84,0.7,0.18])

	axstamp = pyplot.axes([0.51,0.15,0.73,0.73])

	image=snd.gaussian_filter(image,sigma=0.5)
	im = ax5.imshow(image,cmap='viridis',vmin=-0.003,vmax=0.009,origin='lower',aspect='equal',interpolation='none')


	image=snd.gaussian_filter(imageA,sigma=0.5)
	im = ax4.imshow(image,cmap='viridis',vmin=-0.003,vmax=0.009,origin='lower',aspect='equal',interpolation='none')

	image=snd.gaussian_filter(imageB,sigma=0.5)
	im = ax2.imshow(image,cmap='viridis',vmin=-0.003,vmax=0.009,origin='lower',aspect='equal',interpolation='none')

	image=snd.gaussian_filter(scimageA,sigma=0.5)
	im = ax3.imshow(image,cmap='viridis',vmin=-0.003,vmax=0.009,origin='lower',aspect='equal',interpolation='none')
	image=snd.gaussian_filter(scimageB,sigma=0.5)
	im = ax.imshow(image,cmap='viridis',vmin=-0.003,vmax=0.009,origin='lower',aspect='equal',interpolation='none')

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
	


	ax.plot([1,lx-1],[ly/2.,ly/2.],color='white',ls=':',lw=2,alpha=0.6)

	ax2.plot([1,lx-1],[ly/2.,ly/2.],color='white',ls=':',lw=2,alpha=0.6)

	ax3.plot([1,lx-1],[ly/2.,ly/2.],color='white',ls=':',lw=2,alpha=0.6)
	ax4.plot([1,lx-1],[ly/2.,ly/2.],color='white',ls=':',lw=2,alpha=0.6)
	ax5.plot([1,lx-1],[ly/2.,ly/2.],color='white',ls=':',lw=2,alpha=0.6)
	

	#ax.contour(SN,levels=numpy.array([2,3,4,5]),colors='k',linewidths=(1),interpolation='nearest',extent=extent)
	
	#F356W stamp
	#axstamp.plot([thumbsize,thumbsize],[0,thumbsize-8],lw=2,color='white',ls='--')
	#axstamp.plot([0,thumbsize-8],[thumbsize,thumbsize],lw=2,color='white',ls='--')

	ax.tick_params(left=False,bottom=False,top=False,right=False,labelbottom=False,labelleft=False) # Get ticks to look nice
	ax2.tick_params(left=False,bottom=False,top=False,right=False,labelbottom=False,labelleft=False) # Get ticks to look nice
	ax3.tick_params(left=False,bottom=False,top=False,right=False,labelbottom=False,labelleft=False) # Get ticks to look nice
	ax4.tick_params(left=False,bottom=False,top=False,right=False,labelbottom=False,labelleft=False) # Get ticks to look nice
	ax5.tick_params(left=False,bottom=False,top=False,right=False,labelbottom=False,labelleft=False) # Get ticks to look nice



	pyplot.tight_layout()
	#pyplot.savefig('test.png',dpi=120)
	pyplot.savefig('/scratch/EIGER/identification/VISCHECK_NOW/spectrum_J1148_%s.png'%thisID,dpi=120)
	pyplot.clf()


