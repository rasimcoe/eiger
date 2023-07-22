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


def gaussian(x,totflux,c,x0,sigma):
    return totflux*((sigma)**-1 * (2*np.pi)**-0.5 *np.exp(-(x-x0)**2/(2*sigma**2)))+c   




FOLDER='/scratch/EIGER/identification/SPECTRA_O3CANDS_15SEPT/'
SPECTRA_FOLDER='/scratch/EIGER/extraction/OPTIMAL_PROFILES/'
SAVE_FOLDER='/scratch/EIGER/extraction/REDSHIFTS/'

CATALOG='/scratch/EIGER/identification/J0100_photcat_v2_CONCAT_O3candidates_HYBRID_BACKWARD_MANUALRADEC.fits'
SAVE_CATALOG='/scratch/EIGER/identification/J0100_photcat_v2_CONCAT_O3candidates_HYBRID_BACKWARD_MANUALRADEC_withREDSHIFT.fits'

with fits.open(CATALOG) as hdul:
    orig_table = hdul[1].data
    orig_cols = orig_table.columns
    
rescale_noise=True #if true, rescales the mean(err_1d) to be equal to the std(data_1d) with some outlier removal


cat=fits.open(CATALOG)

data=cat[1].data
IDlist=data.field('NUMBER')

z_guesslist=data.field('z_O3_estimate')
Nclumps_Y=data.field('Nclumps_Y')
GOOD=data.field('GOOD_CANDIDATE')

redshift_O3_4960_A=np.zeros(len(IDlist))
redshift_O3_4960_A_err=np.zeros(len(IDlist))
redshift_O3_5008_A=np.zeros(len(IDlist))
redshift_O3_5008_A_err=np.zeros(len(IDlist))

redshift_O3_4960_B=np.zeros(len(IDlist))
redshift_O3_4960_B_err=np.zeros(len(IDlist))
redshift_O3_5008_B=np.zeros(len(IDlist))
redshift_O3_5008_B_err=np.zeros(len(IDlist))



#still under heavy development, I will break up the code in functions asap


for q in range(len(IDlist)):
	thisID=IDlist[q]
	thisz=z_guesslist[q]
	print('now doing id',thisID)
	if GOOD[q]==False:
		continue
	thisNclumps_Y=Nclumps_Y[q]


	dat=fits.open(SPECTRA_FOLDER+'spectrum_1D_%s.fits'%thisID)
	data_1d=dat[1].data

	obs_wav=data_1d.field('wavelength') #in Angstroms
	for module in ['A','B']:
		try:
			flux_tot=data_1d.field('flux_tot_%s'%module)
			flux_tot_err=data_1d.field('flux_tot_%s_err'%module)
		except:
			continue

		model=Model(gaussian,independent_vars=('x'))
		model.set_param_hint('totflux',min=0.,max=20)
		model.set_param_hint('sigma',min=1E-5,max=0.008)
		model.set_param_hint('c',min=-0.5,max=0.5)
		model.set_param_hint('x0',min=thisz-0.02,max=thisz+0.02)


		for line in [4960.295,5008.24]:
			z_array=obs_wav/line -1

			if line==4960.295:
				sel_include=(z_array>thisz-0.09)*(z_array<thisz+0.04)
			if line==5008.24:
				sel_include=(z_array>thisz-0.04)*(z_array<thisz+0.09)


			params=model.make_params(totflux=0.005,c=0,x0=thisz,sigma=0.001)

			result=model.fit(flux_tot[sel_include],x=z_array[sel_include],params=params,weights=1./flux_tot_err[sel_include],nan_policy='propagate')
			

			print(result.fit_report())

			print(result.params['x0'].value,result.params['x0'].stderr)

			z_plot=np.arange(thisz-0.05,thisz+0.05,0.001)
			pyplot.fill_between(z_array[sel_include],-flux_tot_err[sel_include],flux_tot_err[sel_include],lw=0,alpha=0.4,color='tab:blue')
			pyplot.plot(z_array[sel_include],flux_tot[sel_include],color='tab:blue')
			pyplot.plot(z_plot,model.eval(result.params,x=z_plot),lw=2,color='k')
			pyplot.savefig(SAVE_FOLDER+'redshift_fit_%s_%s_mod%s.png'%(thisID,int(line),module))
			pyplot.clf()

			#I know this is ugly, but it does the job
			if module=='A':
				if line==4960.295:
					redshift_O3_4960_A[q]=result.params['x0'].value
					redshift_O3_4960_A_err[q]=result.params['x0'].stderr
				if line==5008.24:
					redshift_O3_5008_A[q]=result.params['x0'].value
					redshift_O3_5008_A_err[q]=result.params['x0'].stderr				
			if module=='B':
				if line==4960.295:
					redshift_O3_4960_B[q]=result.params['x0'].value
					redshift_O3_4960_B_err[q]=result.params['x0'].stderr
				if line==5008.24:
					redshift_O3_5008_B[q]=result.params['x0'].value
					redshift_O3_5008_B_err[q]=result.params['x0'].stderr				

print(redshift_O3_5008_B)





col1 = fits.Column(name='z_O3_4960_A', format='D', array=redshift_O3_4960_A)
col2 = fits.Column(name='z_O3_4960_B', format='D', array=redshift_O3_4960_B)
col3 = fits.Column(name='z_O3_5008_A', format='D', array=redshift_O3_5008_A)
col4 = fits.Column(name='z_O3_5008_B', format='D', array=redshift_O3_5008_B)
col5 = fits.Column(name='z_O3_4960_A_err', format='D', array=redshift_O3_4960_A_err)
col6 = fits.Column(name='z_O3_4960_B_err', format='D', array=redshift_O3_4960_B_err)
col7 = fits.Column(name='z_O3_5008_A_err', format='D', array=redshift_O3_5008_A_err)
col8 = fits.Column(name='z_O3_5008_B_err', format='D', array=redshift_O3_5008_B_err)

new_cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8])
hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)

hdu.writeto(SAVE_CATALOG, overwrite=True)


