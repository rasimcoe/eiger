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

CATALOG='/scratch/EIGER/identification/J1148_list1.fits'
CATALOG='/scratch/EIGER/identification/J1148_O3DK_missed_by_JM.fits'
FOLDER='/scratch/EIGER/identification/SPECTRA_J1148/' #FOLDER WITH SPECTRA

field='J1148'
with fits.open(CATALOG) as hdul:
    orig_table = hdul[1].data
    orig_cols = orig_table.columns
    

IDlist=orig_table.field('NUMBER_2')
RAlist=orig_table.field('ALPHA_J2000_det_2')
DEClist=orig_table.field('DELTA_J2000_det_2')
redshift=orig_table.field('z_O3')
exp_sep=(5008.28-4960.)*(1+redshift)
lamblist=(1+redshift)*5008.28



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
	if lamblist[q]==0.:
		continue
	thisx=(lamblist[q]-30000)/9.75
	
	wav_4960=4960*(lamblist[q]/5008.2)
	thisx_4960=(wav_4960-30000)/9.75
	
	dx=thisx-thisx_4960
	wav_4862=4862*(lamblist[q]/5008.2)
	thisx_4862=(wav_4862-30000)/9.75
	
	dx_Hb=thisx-thisx_4862

	hdu= fits.open(FOLDER+'stacked_2D_J1148_%s.fits'%thisID)
	hd=hdu['EMLINE'].header
	data=hdu['EMLINE'].data

	data_A=hdu['EMLINEA'].data
	data_B=hdu['EMLINEB'].data
	stamp=hdu['STAMP'].data ##F356W image only


	start=int(thisx)-150
	end=int(thisx)+60
	X_5008=150-0.5
	print(start,end)
#	#2D stacks
	image=data[15:-15,start:end]
	imageA=data_A[15:-15,start:end] 
	imageB=data_B[15:-15,start:end]
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
	
	
	
	
	
	
	
	fig=pyplot.figure(figsize=(17.9, 5.7))

	ax = pyplot.axes([0.01,0.06,0.7,0.32])
	ax2 = pyplot.axes([0.01,0.38,0.7,0.32])
	ax3 = pyplot.axes([0.01,0.67,0.7,0.32])

	axstamp = pyplot.axes([0.51,0.15,0.73,0.73])

	image=snd.gaussian_filter(image,sigma=0.5)
	im = ax3.imshow(image,cmap='viridis',vmin=-0.003,vmax=0.009,origin='lower',aspect='equal',interpolation='none')


	image=snd.gaussian_filter(imageA,sigma=0.5)
	im = ax2.imshow(image,cmap='viridis',vmin=-0.003,vmax=0.009,origin='lower',aspect='equal',interpolation='none')

	image=snd.gaussian_filter(imageB,sigma=0.5)
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
	


	ax.plot([X_5008,X_5008],[0,ly/2.-4],color='tab:red',lw=3,ls='-',alpha=0.9)
	ax2.plot([X_5008,X_5008],[0,ly/2.-4],color='tab:red',lw=3,ls='-',alpha=0.9)
	ax3.plot([X_5008,X_5008],[0,ly/2.-4],color='tab:red',lw=3,ls='-',alpha=0.9)

	ax.plot([X_5008-dx,X_5008-dx],[0,ly/2.-4],color='tab:red',lw=2,ls='-',alpha=0.9)
	ax2.plot([X_5008-dx,X_5008-dx],[0,ly/2.-4],color='tab:red',lw=2,ls='-',alpha=0.9)
	ax3.plot([X_5008-dx,X_5008-dx],[0,ly/2.-4],color='tab:red',lw=2,ls='-',alpha=0.9)
	
	ax.plot([X_5008-dx_Hb,X_5008-dx_Hb],[0,ly/2.-4],color='tab:orange',lw=2,ls='-',alpha=0.9)
	ax2.plot([X_5008-dx_Hb,X_5008-dx_Hb],[0,ly/2.-4],color='tab:orange',lw=2,ls='-',alpha=0.9)
	ax3.plot([X_5008-dx_Hb,X_5008-dx_Hb],[0,ly/2.-4],color='tab:orange',lw=2,ls='-',alpha=0.9)

	#ax.contour(SN,levels=numpy.array([2,3,4,5]),colors='k',linewidths=(1),interpolation='nearest',extent=extent)
	
	#F356W stamp
	#axstamp.plot([thumbsize,thumbsize],[0,thumbsize-8],lw=2,color='white',ls='--')
	#axstamp.plot([0,thumbsize-8],[thumbsize,thumbsize],lw=2,color='white',ls='--')

	#ax.tick_params(left=False,bottom=False,top=False,right=False,labelbottom=False,labelleft=False) # Get ticks to look nice



	pyplot.tight_layout
	pyplot.savefig('/scratch/EIGER/identification/VISCHECK_JMmissed/spectrum_J1148_%s.png'%thisID,dpi=120)	
	#pyplot.savefig('/scratch/EIGER/identification/VISCHECK_J1148/spectrum_O3candidate_%s.png'%thisID,dpi=120)
	pyplot.clf()


