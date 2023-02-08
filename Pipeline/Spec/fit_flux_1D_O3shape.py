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
SAVE_FOLDER='/scratch/EIGER/extraction/LINES/'

CATALOG='/scratch/EIGER/extraction/J0100_photcat_v2_O3emitters_CONVERGED_manualZcheck.fits'
SAVE_CATALOG=CATALOG[:-5]+'_withflux.fits'






SAVE_FOLDER='/scratch/EIGER/extraction/LINES_EXTRA/'
SPECTRA_FOLDER='/scratch/EIGER/extraction/OPTIMAL_PROFILES_EXTRA/'

FOLDER='/scratch/EIGER/identification/EXTRA_10OCT/'
CATALOG='/scratch/EIGER/identification/EXTRA_12Oct_24O3.fits'

#CATALOG='/scratch/EIGER/identification/EXTRA_12Oct_FINAL82.fits'

SAVE_CATALOG=CATALOG[:-5]+'_withflux)negallowed.fits'

#SAVE_CATALOG='/scratch/EIGER/identification/J0100_photcat_v2_CONCAT_O3candidates_HYBRID_BACKWARD_MANUALRADEC_withREDSHIFT.fits'

with fits.open(CATALOG) as hdul:
    orig_table = hdul[1].data
    orig_cols = orig_table.columns
    
rescale_noise=True #if true, rescales the mean(err_1d) to be equal to the std(data_1d) with some outlier removal


cat=fits.open(CATALOG)

data=cat[1].data
IDlist=data.field('NUMBER')

z_guesslist=data.field('z_O3doublet_combined')
Nclumps_Y=data.field('Nclumps_Y')
Nclumps_spec=data.field('Nclumps_spec')

flux_5008_A=np.zeros(len(IDlist))
flux_5008_A_err=np.zeros(len(IDlist))
flux_5008_B=np.zeros(len(IDlist))
flux_5008_B_err=np.zeros(len(IDlist))


flux_4960_A=np.zeros(len(IDlist))
flux_4960_A_err=np.zeros(len(IDlist))
flux_4960_B=np.zeros(len(IDlist))
flux_4960_B_err=np.zeros(len(IDlist))


flux_4862_A=np.zeros(len(IDlist))
flux_4862_A_err=np.zeros(len(IDlist))
flux_4862_B=np.zeros(len(IDlist))
flux_4862_B_err=np.zeros(len(IDlist))



#still under heavy development, I will break up the code in functions asap


for q in range(len(IDlist)):
	thisID=IDlist[q]
	thisz=z_guesslist[q]
	print('now doing id',thisID)
	thisNclumps_Y=Nclumps_Y[q]

	thisNclumps_spec=Nclumps_spec[q]

	dat=fits.open(SPECTRA_FOLDER+'spectrum_1D_%s.fits'%thisID)
	data_1d=dat[1].data

	obs_wav=data_1d.field('wavelength') #in Angstroms
	for module in ['A','B']:
		for i in range(3):
			try:
				flux_tot=data_1d.field('flux_tot_%s'%module)
				flux_tot_err=data_1d.field('flux_tot_%s_err'%module)

				if np.nanmax(flux_tot_err)==0:
					flux_tot_err=np.zeros(len(flux_tot_err)) + np.nanstd(flux_tot[(obs_wav> 3.15E4)*(obs_wav<3.95E4)])


				thisline=5008.24
				sel_include=(obs_wav>4975*(1+thisz))*(obs_wav<5050*(1+thisz))

				if thisNclumps_spec==1.:
					model=Model(gaussian,independent_vars=('x'))
					model.set_param_hint('totflux',min=0.,max=200) 
					model.set_param_hint('sigma',min=0.1,max=50)
					model.set_param_hint('c',min=-0.05,max=0.05)
					model.set_param_hint('x0',min=(1+thisz)*thisline - 3,max=(1+thisz)*thisline + 3)
					params=model.make_params(totflux=10,sigma=15.,c=0,x0=(1+thisz)*thisline)

				if thisNclumps_spec==2.:
					model=Model(gaussian,independent_vars=('x'),prefix='m1_') + Model(gaussian,independent_vars=('x'),prefix='m2_') 
					model.set_param_hint('m1_totflux',min=0.1,max=200) 
					model.set_param_hint('m1_sigma',min=0.1,max=26)
					model.set_param_hint('m1_c',min=-0.05,max=0.05)
					model.set_param_hint('m1_x0',min=(1+thisz)*thisline - 24,max=(1+thisz)*thisline + 24)

					model.set_param_hint('m2_totflux',min=0.1,max=200) 
					model.set_param_hint('m2_sigma',min=0.1,max=26)
					model.set_param_hint('m2_c',min=-0.05,max=0.05)
					model.set_param_hint('m2_x0',min=(1+thisz)*thisline - 70,max=(1+thisz)*thisline + 70)			
					params=model.make_params(m1_totflux=10,m1_sigma=15.,m1_c=0,m1_x0=(1+thisz)*thisline,m2_totflux=10,m2_sigma=15.,m2_c=0,m2_x0=(1+thisz)*thisline-5.)

					params['m2_c'].vary=False


				if thisNclumps_spec==3.:
					model=Model(gaussian,independent_vars=('x'),prefix='m1_') + Model(gaussian,independent_vars=('x'),prefix='m2_') + Model(gaussian,independent_vars=('x'),prefix='m3_') 
					model.set_param_hint('m1_totflux',min=0.01,max=200) 
					model.set_param_hint('m1_sigma',min=0.1,max=50)
					model.set_param_hint('m1_c',min=-0.05,max=0.05)
					model.set_param_hint('m1_x0',min=(1+thisz)*thisline - 10,max=(1+thisz)*thisline + 10)

					model.set_param_hint('m2_totflux',min=0.01,max=200) 
					model.set_param_hint('m2_sigma',min=0.1,max=50)
					model.set_param_hint('m2_c',min=-0.05,max=0.05)
					model.set_param_hint('m2_x0',min=(1+thisz)*thisline - 70,max=(1+thisz)*thisline + 70)	


					model.set_param_hint('m3_totflux',min=0.01,max=200) 
					model.set_param_hint('m3_sigma',min=0.1,max=50)
					model.set_param_hint('m3_c',min=-0.05,max=0.05)
					model.set_param_hint('m3_x0',min=(1+thisz)*thisline - 90,max=(1+thisz)*thisline + 90)					
					params=model.make_params(m1_totflux=10,m1_sigma=15.,m1_c=0,m1_x0=(1+thisz)*thisline,m2_totflux=10,m2_sigma=15.,m2_c=0,m2_x0=(1+thisz)*thisline-5.,m3_totflux=10,m3_sigma=15.,m3_c=0,m3_x0=(1+thisz)*thisline+5.)

					params['m2_c'].vary=False
					params['m3_c'].vary=False


				result=model.fit(flux_tot[sel_include],x=obs_wav[sel_include],params=params,weights=1./flux_tot_err[sel_include],nan_policy='propagate',method='differential_evolution')
				

				print(result.fit_report())



				pyplot.plot(obs_wav[sel_include],flux_tot[sel_include],color='tab:blue',drawstyle='steps-mid')
				pyplot.fill_between(obs_wav[sel_include],-flux_tot_err[sel_include],flux_tot_err[sel_include],lw=0,alpha=0.4,color='tab:blue')

				xx=numpy.arange(min(obs_wav[sel_include]),max(obs_wav[sel_include]),0.01)
				pyplot.plot(xx,model.eval(result.params,x=xx),lw=2,color='k',alpha=0.5)
				pyplot.savefig(SAVE_FOLDER+'ID_%s_line_%s_mod%s.png'%(thisID,int(thisline),module))
				pyplot.clf()

				basicresult=copy.deepcopy(result)

				if thisNclumps_spec==1:
					O3_5008_sigma=result.params['sigma'].value
					O3_5008_flux=result.params['totflux'].value
					O3_5008_flux_err=result.params['totflux'].stderr

				if thisNclumps_spec==2:
					O3_5008_sigma_1=result.params['m1_sigma'].value
					O3_5008_sigma_2=result.params['m2_sigma'].value
					O3_5008_flux=result.params['m1_totflux'].value + result.params['m2_totflux'].value
					O3_5008_flux_err=(result.params['m1_totflux'].stderr**2  + result.params['m2_totflux'].stderr**2)**0.5 

				if thisNclumps_spec==3:
					O3_5008_sigma_1=result.params['m1_sigma'].value
					O3_5008_sigma_2=result.params['m2_sigma'].value
					O3_5008_sigma_3=result.params['m3_sigma'].value
					O3_5008_flux=result.params['m1_totflux'].value + result.params['m2_totflux'].value + + result.params['m3_totflux'].value
					O3_5008_flux_err=(result.params['m1_totflux'].stderr**2  + result.params['m2_totflux'].stderr**2 + result.params['m3_totflux'].stderr**2)**0.5 



				#NOW 4960 and Hbeta
				thisline=4862.69#4960.295 #Hbeta:4862.69
				LINES=[4862.69,4960.295]
				SELECTIONS=[(obs_wav>4800*(1+thisz))*(obs_wav<4900*(1+thisz)),(obs_wav>4935*(1+thisz))*(obs_wav<4985*(1+thisz))]

				#SELECTIONS=[(obs_wav>4800*(1+thisz))*(obs_wav<4880*(1+thisz)),(obs_wav>4935*(1+thisz))*(obs_wav<4985*(1+thisz))]

				for qq in range(len(LINES)):
					thisline=LINES[qq]
					sel_include=SELECTIONS[qq]

					if thisNclumps_spec==1.:
						model=Model(gaussian,independent_vars=('x'))
						model.set_param_hint('totflux',min=-10.,max=200) 
						model.set_param_hint('sigma',min=0.1,max=50)
						model.set_param_hint('c',min=-0.05,max=0.05)
						model.set_param_hint('x0',min=(1+thisz)*thisline - 3,max=(1+thisz)*thisline + 3)
						params=model.make_params(totflux=10,sigma=O3_5008_sigma*thisline/5008.24,c=0,x0=(1+thisz)*thisline)
						params['sigma'].vary=False

					if thisNclumps_spec==2.:
						model=Model(gaussian,independent_vars=('x'),prefix='m1_') + Model(gaussian,independent_vars=('x'),prefix='m2_') 
						model.set_param_hint('m1_totflux',min=-10.1,max=200) 
						model.set_param_hint('m1_sigma',min=0.1,max=50)
						model.set_param_hint('m1_c',min=-0.05,max=0.05)
						model.set_param_hint('m1_x0',min=(1+thisz)*thisline - 12,max=(1+thisz)*thisline + 12)

						model.set_param_hint('m2_totflux',min=-10.1,max=200) 
						model.set_param_hint('m2_sigma',min=0.1,max=50)
						model.set_param_hint('m2_c',min=-0.05,max=0.05)
						model.set_param_hint('m2_x0',min=(1+thisz)*thisline - 70,max=(1+thisz)*thisline + 70)			
						params=model.make_params(m1_totflux=10,m1_sigma=O3_5008_sigma_1*thisline/5008.24,m1_c=0,m1_x0=(1+thisz)*thisline,m2_totflux=10,m2_sigma=O3_5008_sigma_2*thisline/5008.24,m2_c=0,m2_x0=(1+thisz)*thisline-5.)

						params['m2_c'].vary=False
						params['m1_sigma'].vary=False
						params['m2_sigma'].vary=False


					if thisNclumps_spec==3.:
						model=Model(gaussian,independent_vars=('x'),prefix='m1_') + Model(gaussian,independent_vars=('x'),prefix='m2_') + Model(gaussian,independent_vars=('x'),prefix='m3_') 
						model.set_param_hint('m1_totflux',min=-10.001,max=200) 
						model.set_param_hint('m1_sigma',min=0.1,max=50)
						model.set_param_hint('m1_c',min=-0.05,max=0.05)
						model.set_param_hint('m1_x0',min=(1+thisz)*thisline - 30,max=(1+thisz)*thisline + 30)

						model.set_param_hint('m2_totflux',min=-10.001,max=200) 
						model.set_param_hint('m2_sigma',min=0.1,max=50)
						model.set_param_hint('m2_c',min=-0.05,max=0.05)
						model.set_param_hint('m2_x0',min=(1+thisz)*thisline - 70,max=(1+thisz)*thisline + 70)	


						model.set_param_hint('m3_totflux',min=-10.001,max=200) 
						model.set_param_hint('m3_sigma',min=0.1,max=50)
						model.set_param_hint('m3_c',min=-0.05,max=0.05)
						model.set_param_hint('m3_x0',min=(1+thisz)*thisline - 90,max=(1+thisz)*thisline + 90)					
						params=model.make_params(m1_totflux=10,m1_sigma=O3_5008_sigma_1*thisline/5008.24,m1_c=0,m1_x0=(1+thisz)*thisline,m2_totflux=10,m2_sigma=O3_5008_sigma_2*thisline/5008.24,m2_c=0,m2_x0=(1+thisz)*thisline-5.,m3_totflux=10,m3_sigma=O3_5008_sigma_3*thisline/5008.24,m3_c=0,m3_x0=(1+thisz)*thisline+5.)

						params['m2_c'].vary=False
						params['m3_c'].vary=False
						params['m1_sigma'].vary=False
						params['m2_sigma'].vary=False
						params['m3_sigma'].vary=False


					result=model.fit(flux_tot[sel_include],x=obs_wav[sel_include],params=params,weights=1./flux_tot_err[sel_include],nan_policy='propagate',method='differential_evolution')

					print(result.fit_report())

					pyplot.plot(obs_wav[sel_include],flux_tot[sel_include],color='tab:blue',drawstyle='steps-mid')
					pyplot.fill_between(obs_wav[sel_include],-flux_tot_err[sel_include],flux_tot_err[sel_include],lw=0,alpha=0.4,color='tab:blue')

					xx=numpy.arange(min(obs_wav[sel_include]),max(obs_wav[sel_include]),0.01)
					pyplot.plot(xx,model.eval(result.params,x=xx),lw=2,color='k',alpha=0.5)
					pyplot.savefig(SAVE_FOLDER+'ID_%s_line_%s_mod%s.png'%(thisID,int(thisline),module))
					pyplot.clf()

					if thisline==4862.69:
						if thisNclumps_spec==1:
							Hb_flux=result.params['totflux'].value
							Hb_flux_err=result.params['totflux'].stderr

						if thisNclumps_spec==2:
							Hb_flux=result.params['m1_totflux'].value + result.params['m2_totflux'].value
							Hb_flux_err=(result.params['m1_totflux'].stderr**2  + result.params['m2_totflux'].stderr**2)**0.5 

						if thisNclumps_spec==3:
							Hb_flux=result.params['m1_totflux'].value + result.params['m2_totflux'].value + + result.params['m3_totflux'].value
							Hb_flux_err=(result.params['m1_totflux'].stderr**2  + result.params['m2_totflux'].stderr**2 + result.params['m3_totflux'].stderr**2)**0.5 

					if thisline==4960.295:
						if thisNclumps_spec==1:
							O3_4960_flux=result.params['totflux'].value
							O3_4960_flux_err=result.params['totflux'].stderr

						if thisNclumps_spec==2:
							O3_4960_flux=result.params['m1_totflux'].value + result.params['m2_totflux'].value
							O3_4960_flux_err=(result.params['m1_totflux'].stderr**2  + result.params['m2_totflux'].stderr**2)**0.5 

						if thisNclumps_spec==3:
							O3_4960_flux=result.params['m1_totflux'].value + result.params['m2_totflux'].value + + result.params['m3_totflux'].value
							O3_4960_flux_err=(result.params['m1_totflux'].stderr**2  + result.params['m2_totflux'].stderr**2 + result.params['m3_totflux'].stderr**2)**0.5 

			except:
				continue
			break

			
		if module == 'A':
			flux_5008_A[q]=O3_5008_flux
			flux_5008_A_err[q]=O3_5008_flux_err

			flux_4960_A[q]=O3_4960_flux
			flux_4960_A_err[q]=O3_4960_flux_err

			flux_4862_A[q]=Hb_flux
			flux_4862_A_err[q]=Hb_flux_err



		if module == 'B':
			flux_5008_B[q]=O3_5008_flux
			flux_5008_B_err[q]=O3_5008_flux_err

			flux_4960_B[q]=O3_4960_flux
			flux_4960_B_err[q]=O3_4960_flux_err

			flux_4862_B[q]=Hb_flux
			flux_4862_B_err[q]=Hb_flux_err

		#need to reset!
		print(thisID,'modA: ')
		print(flux_4862_A[q],flux_4960_A[q],flux_5008_A[q],flux_4862_A_err[q],flux_4960_A_err[q],flux_5008_A_err[q])

		print(thisID,'modB: ')
		print(flux_4862_B[q],flux_4960_B[q],flux_5008_B[q],flux_4862_B_err[q],flux_4960_B_err[q],flux_5008_B_err[q])		
		O3_5008_flux,O3_5008_flux_err,O3_4960_flux,O3_4960_flux_err,Hb_flux,Hb_flux_err=np.nan,np.nan,np.nan,np.nan,np.nan,np.nan


		print('COMB')
		print(0.5*(flux_4862_A[q]+flux_4862_B[q]),0.5*(flux_4960_A[q]+flux_4960_B[q]),0.5*(flux_5008_A[q]+flux_5008_B[q]))


f_4862_COMB=np.zeros(len(flux_5008_A))
f_4862_COMB_err=np.zeros(len(flux_5008_A))

f_4960_COMB=np.zeros(len(flux_5008_A))
f_4960_COMB_err=np.zeros(len(flux_5008_A))

f_5008_COMB=np.zeros(len(flux_5008_A))
f_5008_COMB_err=np.zeros(len(flux_5008_A))

sel_A=flux_5008_A>-1E9
sel_B=flux_5008_B>-1E9
sel_AB=(flux_5008_A>-1E9)*(flux_5008_B>-1E9)

f_4862_COMB[sel_A]=flux_4862_A[sel_A]
f_4862_COMB[sel_B]=flux_4862_B[sel_B]
f_4862_COMB[sel_AB]=0.5*(flux_4862_A[sel_AB]+flux_4862_B[sel_AB])


f_4960_COMB[sel_A]=flux_4960_A[sel_A]
f_4960_COMB[sel_B]=flux_4960_B[sel_B]
f_4960_COMB[sel_AB]=0.5*(flux_4960_A[sel_AB]+flux_4960_B[sel_AB])


f_5008_COMB[sel_A]=flux_5008_A[sel_A]
f_5008_COMB[sel_B]=flux_5008_B[sel_B]
f_5008_COMB[sel_AB]=0.5*(flux_5008_A[sel_AB]+flux_5008_B[sel_AB])



f_4862_COMB_err[sel_A]=flux_4862_A_err[sel_A]
f_4862_COMB_err[sel_B]=flux_4862_B_err[sel_B]
f_4862_COMB_err[sel_AB]=0.5*(flux_4862_A_err[sel_AB]**2+flux_4862_B_err[sel_AB]**2)**0.5


f_4960_COMB_err[sel_A]=flux_4960_A_err[sel_A]
f_4960_COMB_err[sel_B]=flux_4960_B_err[sel_B]
f_4960_COMB_err[sel_AB]=0.5*(flux_4960_A_err[sel_AB]**2+flux_4960_B_err[sel_AB]**2)**0.5


f_5008_COMB_err[sel_A]=flux_5008_A_err[sel_A]
f_5008_COMB_err[sel_B]=flux_5008_B_err[sel_B]
f_5008_COMB_err[sel_AB]=0.5*(flux_5008_A_err[sel_AB]**2+flux_5008_B_err[sel_AB]**2)**0.5


col1 = fits.Column(name='f_Hb', format='D', array=f_4862_COMB)
col2 = fits.Column(name='f_O3_4960', format='D', array=f_4960_COMB)
col3 = fits.Column(name='f_O3_5008', format='D', array=f_5008_COMB)
col4 = fits.Column(name='f_Hb_err', format='D', array=f_4862_COMB_err)
col5 = fits.Column(name='f_O3_4960_err', format='D', array=f_4960_COMB_err)
col6 = fits.Column(name='f_O3_5008_err', format='D', array=f_5008_COMB_err)


#add A and B
col7 = fits.Column(name='f_Hb_A', format='D', array=flux_4862_A)
col8 = fits.Column(name='f_O3_4960_A', format='D', array=flux_4960_A)
col9 = fits.Column(name='f_O3_5008_A', format='D', array=flux_5008_A)
col10 = fits.Column(name='f_Hb_A_err', format='D', array=flux_4862_A_err)
col11 = fits.Column(name='f_O3_4960_A_err', format='D', array=flux_4960_A_err)
col12 = fits.Column(name='f_O3_5008_A_err', format='D', array=flux_5008_A_err)


col13 = fits.Column(name='f_Hb_B', format='D', array=flux_4862_B)
col14 = fits.Column(name='f_O3_4960_B', format='D', array=flux_4960_B)
col15 = fits.Column(name='f_O3_5008_B', format='D', array=flux_5008_B)
col16 = fits.Column(name='f_Hb_B_err', format='D', array=flux_4862_B_err)
col17 = fits.Column(name='f_O3_4960_B_err', format='D', array=flux_4960_B_err)
col18 = fits.Column(name='f_O3_5008_B_err', format='D', array=flux_5008_B_err)




new_cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18])
hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)

hdu.writeto(SAVE_CATALOG, overwrite=True)


