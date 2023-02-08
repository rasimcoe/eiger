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
SAVE_FOLDER='/scratch/EIGER/extraction/REDSHIFTS/'

CATALOG='/scratch/EIGER/identification/J0100_photcat_v2_O3emitters_CONVERGED.fits'
SAVE_CATALOG=CATALOG


SPECTRA_FOLDER='/scratch/EIGER/extraction/OPTIMAL_PROFILES_EXTRA/'
SAVE_FOLDER='/scratch/EIGER/extraction/REDSHIFTS_EXTRA/'

FOLDER='/scratch/EIGER/identification/EXTRA_10OCT/'
CATALOG='/scratch/EIGER/identification/EXTRA_12Oct_24O3.fits'
CATALOG='/scratch/EIGER/identification/EXTRA_12Oct_FINAL82.fits'

SAVE_CATALOG='/scratch/EIGER/identification/EXTRA_12Oct_FINAL82withflux.fits'
#SAVE_CATALOG='/scratch/EIGER/identification/J0100_photcat_v2_CONCAT_O3candidates_HYBRID_BACKWARD_MANUALRADEC_withREDSHIFT.fits'

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
Nclumps_spec=data.field('Nclumps_spec')

redshift_O3doublet_A=np.zeros(len(IDlist))
redshift_O3doublet_A_err=np.zeros(len(IDlist))
redshift_O3doublet_B=np.zeros(len(IDlist))
redshift_O3doublet_B_err=np.zeros(len(IDlist))

fudge_A=np.zeros(len(IDlist))
fudge_B=np.zeros(len(IDlist))

fudge_A_err=np.zeros(len(IDlist))
fudge_B_err=np.zeros(len(IDlist))

#still under heavy development, I will break up the code in functions asap


for q in range(len(IDlist)):
	thisID=IDlist[q]
	thisz=z_guesslist[q]
	print('now doing id',thisID)
	if GOOD[q]==False:
		continue
	thisNclumps_Y=Nclumps_Y[q]

	thisNclumps_spec=Nclumps_spec[q]

	dat=fits.open(SPECTRA_FOLDER+'spectrum_1D_%s.fits'%thisID)
	data_1d=dat[1].data

	obs_wav=data_1d.field('wavelength') #in Angstroms
	for module in ['A','B']:
		try:
			flux_tot=data_1d.field('flux_tot_%s'%module)
			flux_tot_err=data_1d.field('flux_tot_%s_err'%module)
		except:
			continue

		sel_include=(obs_wav>4920*(1+thisz))*(obs_wav<5060*(1+thisz))



		if thisNclumps_spec==1.:
			model=Model(gaussian_O3_withfudge,independent_vars=('x'))
			model.set_param_hint('a',min=0.,max=200) 
			model.set_param_hint('sigma',min=0.01,max=200)
			model.set_param_hint('c',min=-0.5,max=0.5)
			model.set_param_hint('redshift',min=thisz-0.02,max=thisz+0.02)
			model.set_param_hint('fudge',min=0.5,max=2.0)
			params=model.make_params(a=0.5,c=0.,redshift=thisz,sigma=10.,fudge=1.0)


		if thisNclumps_spec==2.: #2 or 4 a the moment
			model=Model(gaussian_O3_withfudge,independent_vars=('x'),prefix='m1_') + Model(gaussian_O3_withfudge,independent_vars=('x'),prefix='m2_')
			model.set_param_hint('m1_a',min=0.,max=200) 
			model.set_param_hint('m1_sigma',min=7.,max=16)
			model.set_param_hint('m1_c',min=-0.5,max=0.5)
			model.set_param_hint('m1_redshift',min=thisz-0.02,max=thisz+0.02)
			model.set_param_hint('m1_fudge',min=0.95,max=1.05)

			model.set_param_hint('m2_a',min=0.,max=200) 
			model.set_param_hint('m2_sigma',min=7.,max=16)
			model.set_param_hint('m2_c',min=-0.5,max=0.5)
			model.set_param_hint('m2_redshift',min=thisz-0.02,max=thisz+0.02)
			model.set_param_hint('m2_fudge',min=0.95,max=1.05)			
			params=model.make_params(m1_a=0.5,m1_c=0.,m1_redshift=thisz+0.005,m1_sigma=10.,m1_fudge=1.0,m2_a=0.5,m2_c=0.,m2_redshift=thisz-0.005,m2_sigma=10.,m2_fudge=1.0)

			params['m2_c'].vary=False
			params['m1_fudge'].vary=False
			params['m2_fudge'].vary=False




		if thisNclumps_spec==3.: #2 or 4 a the moment
			model=Model(gaussian_O3_withfudge,independent_vars=('x'),prefix='m1_') + Model(gaussian_O3_withfudge,independent_vars=('x'),prefix='m2_')+ Model(gaussian_O3_withfudge,independent_vars=('x'),prefix='m3_')
			model.set_param_hint('m1_a',min=0.,max=200) 
			model.set_param_hint('m1_sigma',min=5.,max=35)
			model.set_param_hint('m1_c',min=-0.5,max=0.5)
			model.set_param_hint('m1_redshift',min=thisz-0.04,max=thisz+0.02)
			model.set_param_hint('m1_fudge',min=0.95,max=1.05)

			model.set_param_hint('m2_a',min=0.,max=200) 
			model.set_param_hint('m2_sigma',min=5.,max=35)
			model.set_param_hint('m2_c',min=-0.5,max=0.5)
			model.set_param_hint('m2_redshift',min=thisz-0.02,max=thisz+0.02)
			model.set_param_hint('m2_fudge',min=0.95,max=1.05)		


			model.set_param_hint('m3_a',min=0.,max=200) 
			model.set_param_hint('m3_sigma',min=5.,max=35)
			model.set_param_hint('m3_c',min=-0.5,max=0.5)
			model.set_param_hint('m3_redshift',min=thisz-0.02,max=thisz+0.02)
			model.set_param_hint('m3_fudge',min=0.95,max=1.05)						
			params=model.make_params(m1_a=0.5,m1_c=0.,m1_redshift=thisz+0.005,m1_sigma=10.,m1_fudge=1.0,m2_a=0.5,m2_c=0.,m2_redshift=thisz-0.005,m2_sigma=10.,m2_fudge=1.0,m3_a=0.5,m3_c=0.,m3_redshift=thisz-0.00,m3_sigma=10.,m3_fudge=1.0)

			params['m2_c'].vary=False
			params['m3_c'].vary=False

			params['m1_fudge'].vary=False
			params['m2_fudge'].vary=False
			params['m3_fudge'].vary=False




		result=model.fit(flux_tot[sel_include],x=obs_wav[sel_include],params=params,weights=1./flux_tot_err[sel_include],nan_policy='propagate')
		

		print(result.fit_report())
		if thisNclumps_spec==1:
			redshift,redshifterr=result.params['redshift'].value,result.params['redshift'].stderr

		if thisNclumps_spec>1:
			redshift1,redshift1err=result.params['m1_redshift'].value,result.params['m1_redshift'].stderr
			redshift2,redshift2err=result.params['m2_redshift'].value,result.params['m2_redshift'].stderr
			tot1=result.params['m1_a']*result.params['m1_sigma']  
			tot2=result.params['m2_a']*result.params['m2_sigma']
			if tot1>tot2:
				redshift,redshifterr=redshift1,redshift1err
			else:
				redshift,redshifterr=redshift2,redshift2err
			print('TOT1,TOT2',tot1,tot2,redshift1,redshift2)



		pyplot.plot(obs_wav[sel_include],flux_tot[sel_include],color='tab:blue')
		pyplot.fill_between(obs_wav[sel_include],-flux_tot_err[sel_include],flux_tot_err[sel_include],lw=0,alpha=0.4,color='tab:blue')

		xx=numpy.arange(4920*(1+thisz),5060*(1+thisz),1.)
		pyplot.plot(xx,model.eval(result.params,x=xx),lw=2,color='k')
		pyplot.savefig(SAVE_FOLDER+'redshift_fit_O3doublet_%s_mod%s.png'%(thisID,module))
		pyplot.clf()
		

		if module == 'A':
			redshift_O3doublet_A[q]=redshift
			redshift_O3doublet_A_err[q]=redshifterr
			fudge_A[q]=result.params['fudge'].value
			fudge_A_err[q]=result.params['fudge'].stderr

		if module == 'B':
			redshift_O3doublet_B[q]=redshift
			redshift_O3doublet_B_err[q]=redshifterr
			fudge_B[q]=result.params['fudge'].value
			fudge_B_err[q]=result.params['fudge'].stderr


col1 = fits.Column(name='z_O3doublet_A_n', format='D', array=redshift_O3doublet_A)
col2 = fits.Column(name='z_O3doublet_B_n', format='D', array=redshift_O3doublet_B)

col3 = fits.Column(name='z_O3doublet_A_err_n', format='D', array=redshift_O3doublet_A_err)
col4 = fits.Column(name='z_O3doublet_B_err_n', format='D', array=redshift_O3doublet_B_err)



z_combined=np.zeros(len(redshift_O3doublet_A))

sel_A=redshift_O3doublet_A>0
sel_B=redshift_O3doublet_B>0
sel_AB=(redshift_O3doublet_A>0)*(redshift_O3doublet_B>0)

z_combined[sel_A]=redshift_O3doublet_A[sel_A]
z_combined[sel_B]=redshift_O3doublet_B[sel_B]
z_combined[sel_AB]=0.5*(redshift_O3doublet_A[sel_AB]+redshift_O3doublet_B[sel_AB])


col5 = fits.Column(name='z_O3doublet_combined_n', format='D', array=z_combined)

col6 = fits.Column(name='fudge_A', format='D', array=fudge_A)
col7 = fits.Column(name='fudge_A_err', format='D', array=fudge_A_err)

col8 = fits.Column(name='fudge_B', format='D', array=fudge_B)
col9 = fits.Column(name='fudge_B_err', format='D', array=fudge_B_err)



new_cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9])

hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)

hdu.writeto(SAVE_CATALOG, overwrite=True)


