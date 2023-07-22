import numpy
from astropy.io import fits as pyfits #no need to call it pyfits, sorry about that
from matplotlib import pyplot
import numpy
from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import Planck18
#cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
cosmo=Planck18
from scipy.interpolate import interp1d
from astropy.convolution import Gaussian2DKernel
import scipy.ndimage as snd
from astropy.convolution import convolve
import copy
from scipy.integrate import simps

import matplotlib
import numpy as np

from astropy.io import fits


import lmfit
from lmfit import Model


##JM PLOT STYLE:
from matplotlib import pyplot

pyplot.rcParams['xtick.labelsize']=17
pyplot.rcParams['ytick.labelsize']=17
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



def gaussian_HaN2_onecomp(x,flux_Ha_narrow,N2Ha_narrow,c,redshift_narrow,sigma_narrow):
    x0_Ha_n=(1+redshift_narrow)*6564.633
    x0_N2r_n=(1+redshift_narrow)*6585.42
    x0_N2b_n=(1+redshift_narrow)*6549.91

    
    Halpha_N=flux_Ha_narrow*(sigma_narrow)**-1 * (2*np.pi)**-0.5*numpy.exp(-(x-x0_Ha_n)**2/(2*sigma_narrow**2))

    
    N2b_N=2.94**-1 * N2Ha_narrow*flux_Ha_narrow*(sigma_narrow)**-1 * (2*np.pi)**-0.5*numpy.exp(-(x-x0_N2b_n)**2/(2*sigma_narrow**2))

  
    N2r_N=N2Ha_narrow*flux_Ha_narrow*(sigma_narrow)**-1 * (2*np.pi)**-0.5*numpy.exp(-(x-x0_N2r_n)**2/(2*sigma_narrow**2))

    

    return c+Halpha_N+N2b_N+N2r_N


def gaussian_HaN2_twocomp(x,flux_Ha_narrow,flux_Ha_broad,N2Ha_narrow,N2Ha_broad,c,redshift_narrow,redshift_broad,sigma_narrow,sigma_broad):
    redshift_broad=redshift_narrow
    x0_Ha_n=(1+redshift_narrow)*6564.633
    x0_Ha_b=(1+redshift_broad)*6564.633    

    x0_N2r_n=(1+redshift_narrow)*6585.42
    x0_N2r_b=(1+redshift_broad)*6585.42 

    x0_N2b_n=(1+redshift_narrow)*6549.91
    x0_N2b_b=(1+redshift_broad)*6549.91
    
    Halpha_N=flux_Ha_narrow*(sigma_narrow)**-1 * (2*np.pi)**-0.5*numpy.exp(-(x-x0_Ha_n)**2/(2*sigma_narrow**2))
    Halpha_B=flux_Ha_broad*(sigma_broad)**-1 * (2*np.pi)**-0.5*numpy.exp(-(x-x0_Ha_b)**2/(2*sigma_broad**2))    
    
    N2b_N=2.94**-1 * N2Ha_narrow*flux_Ha_narrow*(sigma_narrow)**-1 * (2*np.pi)**-0.5*numpy.exp(-(x-x0_N2b_n)**2/(2*sigma_narrow**2))
    N2b_B=2.94**-1 * N2Ha_broad*flux_Ha_broad*(sigma_broad)**-1 * (2*np.pi)**-0.5*numpy.exp(-(x-x0_N2b_b)**2/(2*sigma_broad**2))    
  
    N2r_N=N2Ha_narrow*flux_Ha_narrow*(sigma_narrow)**-1 * (2*np.pi)**-0.5*numpy.exp(-(x-x0_N2r_n)**2/(2*sigma_narrow**2))
    N2r_B=N2Ha_broad*flux_Ha_broad*(sigma_broad)**-1 * (2*np.pi)**-0.5*numpy.exp(-(x-x0_N2r_b)**2/(2*sigma_broad**2))         
    

    return c+Halpha_N+Halpha_B+N2b_N+N2b_B+N2r_N+N2r_B


#totflux*((sigma)**-1 * (2*np.pi)**-0.5

#Halpha 6564.633
#NII 6585.42, 6549.91  #Brighter line is a factor 2.94 higher


FOLDER='/scratch/EIGER/BROAD/SPECTRA_COLSEL/' #FOLDER WITH SPECTRA
SAVE_FOLDER='/scratch/EIGER/BROAD/EXTRACTIONS/'


CATALOG='/scratch/EIGER/BROAD/BROADsel_allfields_23022023_zguess_classed.fits'

SAVE_CATALOG='/scratch/EIGER/BROAD/BROADsel_allfields_23022023_zguess_classed_withfits.fits'
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
HIGHC=orig_table.field('HIGH_CONFID')
MULTI=orig_table.field('Multiply_totalnarrow')
KERNELlist=orig_table.field('KERNEL')


FLUX_NARROW=numpy.zeros(len(IDlist))
FWHM_NARROW=numpy.zeros(len(IDlist))
REDSHIFT_NARROW=numpy.zeros(len(IDlist))
N2Ha_NARROW=numpy.zeros(len(IDlist))

FLUX_BROAD=numpy.zeros(len(IDlist))
FWHM_BROAD=numpy.zeros(len(IDlist))
REDSHIFT_BROAD=numpy.zeros(len(IDlist))
N2Ha_BROAD=numpy.zeros(len(IDlist))

FLUX_NARROW_ERR=numpy.zeros(len(IDlist))
FWHM_NARROW_ERR=numpy.zeros(len(IDlist))
REDSHIFT_NARROW_ERR=numpy.zeros(len(IDlist))
N2Ha_NARROW_ERR=numpy.zeros(len(IDlist))

FLUX_BROAD_ERR=numpy.zeros(len(IDlist))
FWHM_BROAD_ERR=numpy.zeros(len(IDlist))
REDSHIFT_BROAD_ERR=numpy.zeros(len(IDlist))
N2Ha_BROAD_ERR=numpy.zeros(len(IDlist))


L_HA_NARROW=numpy.zeros(len(IDlist))
L_HA_BROAD=numpy.zeros(len(IDlist))
SMBH_MASS=numpy.zeros(len(IDlist))
SFR_NARROW=numpy.zeros(len(IDlist))





for q in [4]:#range(len(IDlist)):
	thisField=FIELDlist[q]
	thisz=Zlist[q]
	thisMulti=MULTI[q]
	
	if HIGHC[q]==False:
		continue

	thisID=IDlist[q]
	if KERNELlist[q]=='wide':
		thisID=str(thisID)+'_widekernel'
	if KERNELlist[q]=='mid':
		thisID=str(thisID)+'_midkernel'
	
	thisNclumps=Nlist[q]
	thisMod=MODlist[q]
	print('Now doing',thisField,thisID)
	if thisMod=='BOTH':
		thisMod=''
	dat=fits.open(SAVE_FOLDER+'cleaned_spectrum_1D_%s_%s.fits'%(thisField,thisID))
	data_1d=dat[1].data

	obs_wav=data_1d.field('wavelength') #in Angstroms
	flux_tot=data_1d.field('flux')
	flux_tot_err=data_1d.field('flux_err')



	###Load 2D 
	data_2D=fits.getdata(SAVE_FOLDER+'cleaned_2d_for_extraction_%s_%s.fits'%(thisField,thisID))
	hd=fits.getheader(SAVE_FOLDER+'cleaned_2d_for_extraction_%s_%s.fits'%(thisField,thisID))
	ly,lx=np.shape(data_2D)
	x_array=np.arange(0,lx,1)
	wav_array=x_array*hd['CDELT1'] + hd['CRVAL1']









	sel_include=(obs_wav>6440*(1+thisz))*(obs_wav<6690*(1+thisz))
	sel_plot=(obs_wav>6480*(1+thisz))*(obs_wav<6650*(1+thisz))
	sel_plot=(obs_wav>6450*(1+thisz))*(obs_wav<6680*(1+thisz))		
		
		
	#FIRST DO A SINGLE COMPONENT FIT	

	onemodel=Model(gaussian_HaN2_onecomp,independent_vars=('x')) #prefix='m1_'
	onemodel.set_param_hint('flux_Ha_narrow',min=0.,max=200.)
	onemodel.set_param_hint('N2Ha_narrow',min=0.,max=0.5)
	onemodel.set_param_hint('c',min=-0.5,max=0.5)
	onemodel.set_param_hint('redshift_narrow',min=thisz-0.02,max=thisz+0.02)
	onemodel.set_param_hint('sigma_narrow',min=7.,max=150.)

	params=onemodel.make_params(flux_Ha_narrow=10.,N2Ha_narrow=0.02,c=0,redshift_narrow=thisz,sigma_narrow=12.)
	result_one=onemodel.fit(flux_tot[sel_include],x=obs_wav[sel_include],params=params,weights=1./flux_tot_err[sel_include],nan_policy='propagate',method='ampgo')
	
			
	chisq_one=result_one.redchi
		
		
	#TWO COMPONENT	

	
	#model=Model(gaussian_HaN2_twocomp,independent_vars=('x')) #prefix='m1_'
	#model.set_param_hint('flux_Ha_narrow',min=0.,max=20.)
	#model.set_param_hint('flux_Ha_broad',min=0.,max=100.)	
	#model.set_param_hint('N2Ha_narrow',min=0.,max=0.5)
	#model.set_param_hint('N2Ha_broad',min=0.,max=0.5)
	#model.set_param_hint('c',min=-0.5,max=0.5)
	#model.set_param_hint('redshift_narrow',min=thisz-0.012,max=thisz+0.012)
	#model.set_param_hint('redshift_broad',min=thisz-0.02,max=thisz+0.02)
	#model.set_param_hint('sigma_narrow',min=7.,max=35.)
	#model.set_param_hint('sigma_broad',min=30.,max=300.)
	#params=model.make_params(flux_Ha_narrow=0.2,flux_Ha_broad=2.,N2Ha_narrow=0.02,N2Ha_broad=0.02,c=0, \
#redshift_narrow=thisz,redshift_broad=thisz,sigma_narrow=12.,sigma_broad=70.)
#(x,flux_Ha_narrow,flux_Ha_broad,N2Ha_narrow,N2Ha_broad,c,redshift_narrow,redshift_broad,sigma_narrow,sigma_broad
	#params['redshift_broad'].vary=False
	
	
	
		
	model=Model(gaussian_HaN2_twocomp,independent_vars=('x')) #prefix='m1_'  	redshift_ABS,flux_Ha_ABS,sigma_ABS
	model.set_param_hint('flux_Ha_narrow',min=0.,max=20.)
	model.set_param_hint('flux_Ha_broad',min=0.,max=100.)	
	model.set_param_hint('flux_Ha_ABS',min=-10.,max=0.)		
	model.set_param_hint('N2Ha_narrow',min=0.,max=0.5)
	model.set_param_hint('N2Ha_broad',min=0.,max=0.5)
	model.set_param_hint('c',min=-0.5,max=0.5)
	model.set_param_hint('redshift_narrow',min=thisz-0.012,max=thisz+0.012)
	model.set_param_hint('redshift_broad',min=thisz-0.02,max=thisz+0.02)
	model.set_param_hint('redshift_ABS',min=thisz-0.04,max=thisz+0.04)	

	
	params=model.make_params(flux_Ha_narrow=0.2,flux_Ha_broad=2.,N2Ha_narrow=0.0,N2Ha_broad=0.0,c=0, \
redshift_narrow=thisz,redshift_broad=5.01267,sigma_narrow=12.,sigma_broad=70.)
#(x,flux_Ha_narrow,flux_Ha_broad,N2Ha_narrow,N2Ha_broad,c,redshift_narrow,redshift_broad,sigma_narrow,sigma_broad

	params['redshift_broad'].vary=False
	#params['redshift_narrow'].vary=False
	#params['N2Ha_narrow'].vary=False
	#params['N2Ha_broad'].vary=False




	result=model.fit(flux_tot[sel_include],x=obs_wav[sel_include],params=params,weights=1./flux_tot_err[sel_include],nan_policy='propagate',method='ampgo')
	
	
	chisq_two=result.redchi	
	best_z=result.params['redshift_narrow'].value
	
	print(result.fit_report())
	print(result_one.fit_report())	

	fig=pyplot.figure(figsize=(6.7,5.5))
	ax = pyplot.axes([0.145,0.28,0.80,0.58])

	axinset = pyplot.axes([0.145,0.11,0.80,0.16])
	axinset2 = pyplot.axes([0.145,0.87,0.80,0.12])


	imageA=data_2D[17:-17,sel_plot]
	print(numpy.shape(imageA))

	ly,lx=numpy.shape(imageA)
	imageA=snd.gaussian_filter(imageA,sigma=0.5)

	im = axinset2.imshow(imageA,cmap='viridis',vmin=-0.006,vmax=0.016,origin='lower',aspect='auto',interpolation='nearest')
	axinset2.plot([0,9],[ly/2.,ly/2.],color='white',ls=':',lw=2,alpha=0.9)
	axinset2.plot([lx-10,lx-1],[ly/2.,ly/2.],color='white',ls=':',lw=2,alpha=0.9)

	flux_tot=flux_tot*(1+best_z)
	flux_tot_err=flux_tot_err*(1+best_z)


	dv=3E5*(obs_wav[sel_plot]/(1+best_z)-6564.633)/6564.633
	
	ax.plot(dv,flux_tot[sel_plot],drawstyle='steps-mid',color='tab:blue',lw=2,label='%s-%s'%(thisField,thisID))
	ax.fill_between(dv,-flux_tot_err[sel_plot],flux_tot_err[sel_plot],lw=0,alpha=0.4,color='tab:blue')

	xx=numpy.arange(min(obs_wav[sel_plot]),max(obs_wav[sel_plot]),1.)
	xxplot=xx/(1+best_z)
	xxplot=3E5*(xxplot-6564.633)/6564.633
	
	
	#Sub-components:
	result_narrow=copy.deepcopy(result)
	result_narrow_Ha=copy.deepcopy(result)
	result_broad=copy.deepcopy(result)
	result_broad_Ha=copy.deepcopy(result)

	result_narrow.params['flux_Ha_broad'].value=0.

	
	result_broad.params['flux_Ha_narrow'].value=0.

	
	result_narrow_Ha.params['flux_Ha_broad'].value=0.
	result_narrow_Ha.params['N2Ha_narrow'].value=0.

	
	result_broad_Ha.params['flux_Ha_narrow'].value=0.
	result_broad_Ha.params['N2Ha_broad'].value=0.	

	
			
	result_ABS.params['flux_Ha_broad'].value=0.	
	result_ABS.params['flux_Ha_narrow'].value=0.	
	result_ABS.params['N2Ha_narrow'].value=0.
	result_ABS.params['N2Ha_broad'].value=0.
			
	mymodel=model.eval(result.params,x=xx)
	mymodel_narrow_Ha=model.eval(result_narrow_Ha.params,x=xx)	
	mymodel_narrow=model.eval(result_narrow.params,x=xx)
	mymodel_broad_Ha=model.eval(result_broad_Ha.params,x=xx)	
	mymodel_broad=model.eval(result_broad.params,x=xx)
	mymodel_ABS=model.eval(result_ABS.params,x=xx)
		
	mymodel_single=onemodel.eval(result_one.params,x=xx)


	ax.plot(xxplot,(1+best_z)*mymodel_narrow_Ha,lw=1,color='tab:red',label=r'H$\alpha$')
	ax.plot(xxplot,(1+best_z)*mymodel_broad_Ha,lw=1,color='tab:red',ls='--')	
	
	ax.plot(xxplot,(1+best_z)*(mymodel_narrow-mymodel_narrow_Ha),lw=1,color='tab:green',label='[NII]')
	ax.plot(xxplot,(1+best_z)*(mymodel_broad-mymodel_broad_Ha),lw=1,color='tab:green',ls='--')	
		
	ax.plot(xxplot,(1+best_z)*model.eval(result.params,x=xx),lw=2,color='k')	
	

		#DUMMYS
	ax.plot(obs_wav[sel_plot]/(1+best_z) +1000,flux_tot[sel_plot],drawstyle='steps-mid',color='grey',lw=1,label='Narrow')
	ax.plot(obs_wav[sel_plot]/(1+best_z) +1000,flux_tot[sel_plot],drawstyle='steps-mid',color='grey',lw=1,ls='--',label='Broad')


	
	
	
	mean_noise=numpy.nanmean(flux_tot_err[sel_plot])
	
	
	axinset.fill_between(dv,-flux_tot_err[sel_plot]/mean_noise,flux_tot_err[sel_plot]/mean_noise,lw=0,alpha=0.4,color='tab:blue')	
	axinset.plot(dv,numpy.zeros(len(obs_wav[sel_plot])),color='k',ls=':',lw=2,alpha=0.8)	
	axinset.plot(dv,(flux_tot[sel_plot]-(1+best_z)*onemodel.eval(result_one.params,x=obs_wav[sel_plot]))/mean_noise,lw=2,color='tab:purple',drawstyle='steps-mid')	
	axinset.plot(dv,(flux_tot[sel_plot]-(1+best_z)*model.eval(result.params,x=obs_wav[sel_plot]))/mean_noise,lw=2,color='tab:blue',drawstyle='steps-mid')
	
	
	
		
	ax.minorticks_on()
	axinset.minorticks_on()

	ax.tick_params(which='both',direction='in',right=True,top=True,labelbottom=False) # Get ticks to look nice


	axinset.tick_params(which='both',direction='in',right=False,top=False) # Get ticks to look nice
	axinset2.tick_params(which='both',direction='in',right=False,top=False,bottom=False,labelbottom=False,labelleft=False,left=False) # Get ticks to look nice

	axinset.text(-4500,3.9,r'One component $\chi_{\nu}^2$ =%.2f'%chisq_one,color='tab:purple',fontsize=12)
	axinset.text(800,3.9,r'Two components $\chi_{\nu}^2$ =%.2f'%chisq_two,color='tab:blue',fontsize=12)
	
	#axinset.set_xlabel(r'$\lambda_{\rm 0}$ [${\AA}$]',fontsize=20,labelpad=-5)
	axinset.set_xlabel(r'$\Delta v$ [km s$^{-1}$]',fontsize=20,labelpad=-5)	
	ax.set_ylabel(r'$f_{\lambda}$ [$10^{-18}$ erg s$^{-1}$ cm$^{-2}$ ${\AA}^{-1}$]',fontsize=17)
	axinset.set_ylabel(r'$\Delta \sigma$',fontsize=17)
	#ax.set_ylim(-0.04,0.34)
	ax.set_xlim(6480,6650)
	axinset.set_xlim(6480,6650)
	
	ax.set_xlim(-5200,5200)
	axinset.set_xlim(-5200,5200)	
	axinset.set_ylim(-7,9)
	ax.legend(loc=2,frameon=False,fontsize=13)
	pyplot.tight_layout()
	#pyplot.show()
	#STIOPPPPP
	pyplot.savefig(SAVE_FOLDER+'profile_fit_%s_%s.png'%(thisField,IDlist[q]))
	pyplot.clf()
	
	
	
	#SAVE VALUES:
	FLUX_NARROW[q]=thisMulti*result.params['flux_Ha_narrow'].value
	FWHM_NARROW[q]=3E5*2.355*result.params['sigma_narrow'].value /((1+result.params['redshift_narrow'].value) * 6564.633)
	REDSHIFT_NARROW[q]=result.params['redshift_narrow'].value
	N2Ha_NARROW[q]=result.params['N2Ha_narrow'].value

	FLUX_BROAD[q]=result.params['flux_Ha_broad'].value
	fwhm_broad=3E5*2.355*result.params['sigma_broad'].value/((1+result.params['redshift_broad'].value) * 6564.633)
	FWHM_BROAD[q]=fwhm_broad
	REDSHIFT_BROAD[q]=result.params['redshift_broad'].value
	N2Ha_BROAD[q]=result.params['N2Ha_broad'].value

	FLUX_NARROW_ERR[q]=thisMulti*result.params['flux_Ha_narrow'].stderr
	FWHM_NARROW_ERR[q]=3E5*2.355*result.params['sigma_narrow'].stderr /((1+result.params['redshift_narrow'].value) * 6564.633)
	REDSHIFT_NARROW_ERR[q]=result.params['redshift_narrow'].stderr
	N2Ha_NARROW_ERR[q]=result.params['N2Ha_narrow'].stderr

	FLUX_BROAD_ERR[q]=result.params['flux_Ha_broad'].stderr
	FWHM_BROAD_ERR[q]=3E5*2.355*result.params['sigma_broad'].stderr /((1+result.params['redshift_broad'].value) * 6564.633)
	REDSHIFT_BROAD_ERR[q]=result.params['redshift_broad'].stderr
	N2Ha_BROAD_ERR[q]=result.params['N2Ha_broad'].stderr

	##ADD Luminosity, SFR, SMBH mass
	DL=cosmo.luminosity_distance(result.params['redshift_narrow'].value).value
	print('z',result.params['redshift_narrow'].value)
	 
	L_Ha_Narrow=thisMulti*4*3.1415*(DL*3.0856E24)**2 * result.params['flux_Ha_narrow'].value * 1E-18
	L_Ha_Broad=4*3.1415*(DL*3.0856E24)**2 * result.params['flux_Ha_broad'].value * 1E-18	
		
	logMBH=6.57+numpy.log10(1.075) + 0.47*numpy.log10(L_Ha_Broad/1E42) +2.06*numpy.log10(fwhm_broad/1000.) #1.075= epsilon Reines&Volonteri 2015
	SFR=10**-41.27 * L_Ha_Narrow ##Kennicutt Evans 2012 conversion
	print('DL',DL,'L_Ha N,B',L_Ha_Narrow,L_Ha_Broad,'FWHMS',fwhm_broad)
	print('SMBH, SFR',logMBH,SFR)
	
	L_HA_NARROW[q]=L_Ha_Narrow
	L_HA_BROAD[q]=L_Ha_Broad
	SMBH_MASS[q]=logMBH
	SFR_NARROW[q]=SFR
	STOP



col1 = fits.Column(name='f_Ha_narrow', format='D', array=FLUX_NARROW)
col2 = fits.Column(name='f_Ha_narrow_err', format='D', array=FLUX_NARROW_ERR)
col3 = fits.Column(name='f_Ha_broad', format='D', array=FLUX_BROAD)
col4 = fits.Column(name='f_Ha_broad_err', format='D', array=FLUX_BROAD_ERR)
col5 = fits.Column(name='z_Ha_narrow', format='D', array=REDSHIFT_NARROW)
col6 = fits.Column(name='z_Ha_narrow_err', format='D', array=REDSHIFT_NARROW_ERR)
col7 = fits.Column(name='z_Ha_broad', format='D', array=REDSHIFT_BROAD)
col8 = fits.Column(name='z_Ha_broad_err', format='D', array=REDSHIFT_BROAD_ERR)
col9 = fits.Column(name='FWHM_narrow', format='D', array=FWHM_NARROW)
col10 = fits.Column(name='FWHM_narrow_err', format='D', array=FWHM_NARROW_ERR)
col11 = fits.Column(name='FWHM_broad', format='D', array=FWHM_BROAD)
col12 = fits.Column(name='FWHM_broad_err', format='D', array=FWHM_BROAD_ERR)
col13 = fits.Column(name='N2Ha_narrow', format='D', array=N2Ha_NARROW)
col14 = fits.Column(name='N2Ha_narrow_err', format='D', array=N2Ha_NARROW_ERR)
col15 = fits.Column(name='N2Ha_broad', format='D', array=N2Ha_BROAD)
col16 = fits.Column(name='N2Ha_broad_err', format='D', array=N2Ha_BROAD_ERR)
col17 = fits.Column(name='LHA_NARROW', format='D', array=L_HA_NARROW)
col18 = fits.Column(name='LHA_BROAD', format='D', array=L_HA_BROAD)
col19 = fits.Column(name='MBH', format='D', array=SMBH_MASS)
col20 = fits.Column(name='SFR', format='D', array=SFR_NARROW)





new_cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18,col19,col20])

hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)

hdu.writeto(SAVE_CATALOG, overwrite=True)



