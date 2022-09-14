import numpy
import astropy
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from matplotlib import pyplot
from astropy.visualization import simple_norm
import grismconf
import numpy as np
from jwst import datamodels
from jwst.pipeline import calwebb_image2
import os
from lmfit import Model
import time
from scipy import stats
import urllib
import shutil
from astropy.io import fits
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.stats import SigmaClip, sigma_clipped_stats
from scipy import ndimage
from astropy.coordinates import SkyCoord
from astropy import units as u
from stcal.dqflags import interpret_bit_flags
from jwst.resample.resample_utils import build_mask
from scipy import ndimage
from photutils.segmentation import detect_sources, SourceCatalog
from photutils.utils import circular_footprint
from jwst.pipeline import Image3Pipeline
import warnings
warnings.filterwarnings("ignore")

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return numpy.isnan(y), lambda z: z.nonzero()[0]


def xylamb_dispersed_to_radec(ratefile,x_line,y_line,obs_lamb_line,configfile='../simulate_mirage/MIRAGE_DATA/mirage_data/nircam/GRISM_NIRCAM/current/NIRCAM_F356W_modA_R.conf'):
	#delta_lamb_int =0.0144 O3_5008 to Hbeta
	#ratefile='jw01243001001_01101_00001_nrca5_flatfieldstep.fits'
	grism_with_wcs = datamodels.open(ratefile)

	pix_to_world = grism_with_wcs.meta.wcs.get_transform('detector','world')

	C = grismconf.Config(configfile)
	t= C.INVDISPL('+1',0,0,obs_lamb_line)
	dx=C.DISPX('+1',0,0,t) 
	dy=C.DISPY('+1',0,0,t) 
	x0=x_line-dx
	y0=y_line-dy
	RA,DEC,a,b=pix_to_world(x0,y0,0,0)
	return RA,DEC


def pix_dispersed_guess_doublet_to_radec(ratefile,x_line,y_line,lamb_line,delta_x_obs,delta_z_min=0.16,delta_z_max=0.16,delta_z=0.04,delta_lamb_int=0.0144,configfile='../simulate_mirage/MIRAGE_DATA/mirage_data/nircam/GRISM_NIRCAM/current/NIRCAM_F356W_modA_R.conf'):
	#delta_lamb_int =0.0144 O3_5008 to Hbeta
	#ratefile='jw01243001001_01101_00001_nrca5_flatfieldstep.fits'
	grism_with_wcs = datamodels.open(ratefile)

	pix_to_world = grism_with_wcs.meta.wcs.get_transform('detector','world')


	C = grismconf.Config(configfile)
	dLdx=C.DDISPL('+1',0,0,0)/C.DDISPX('+1',0,0,0)

	zguess=-1+  (delta_x_obs*dLdx)/delta_lamb_int

	print('Redshift guess:',zguess)
	l_obs_line=lamb_line * (1+zguess)
	t= C.INVDISPL('+1',0,0,l_obs_line)
	dx=C.DISPX('+1',0,0,t) 
	dy=C.DISPY('+1',0,0,t) 
	x0=x_line-dx
	y0=y_line-dy
	RA,DEC,a,b=pix_to_world(x0,y0,0,0)

	#now a range:
	redshift_guesses=numpy.arange(zguess-delta_z_min,zguess+delta_z_max,delta_z)
	l_obs_line=lamb_line * (1+redshift_guesses)
	t= C.INVDISPL('+1',0,0,l_obs_line)
	dx=C.DISPX('+1',0,0,t) 
	dy=C.DISPY('+1',0,0,t) 
	x0=x_line-dx
	y0=y_line-dy
	RAlist,DEClist,a,b=pix_to_world(x0,y0,0,0)

	return RA,DEC,RAlist,DEClist,redshift_guesses


def create_cutout_circlelist(directimage,RA,DEC,RAlist,DEClist,cutout_filename='STAMP.fits',stampsize=200):
	#directimage='reduced/test_10bright_mag26_z4to7/imaging_F356W/step_i2d.fits'
	hdu = fits.open(directimage)
	wcs = WCS(hdu['SCI'].header)
	positions=SkyCoord(RA,DEC,unit="deg")
	size = (stampsize,stampsize)
	cutout = Cutout2D(hdu['SCI'].data, position=positions, size=size, wcs=wcs)
	hdu['SCI'].data = cutout.data
	hdu['SCI'].header.update(cutout.wcs.to_header())
	fits.writeto(cutout_filename,cutout.data,hdu['SCI'].header,overwrite=True)


	w=WCS(hdu['SCI'].header)
	clist=SkyCoord(RAlist,DEClist,unit="deg")

	xlist,ylist=w.world_to_pixel(clist)
	return cutout.data,xlist,ylist


def plot_stamp(data,xlist,ylist,redshift_guesses,filename='stamp.pdf',CMAP='viridis',percentile=99):
	#data=2D fits file that is the stamp. 
	#list,ylist,redshift_guesses: the pixel positions where the source should lie for various redshifts
	fig=pyplot.figure(1, figsize=(7.4, 6.2))
	ax = pyplot.axes([0.08,0.1,0.9,0.9])
	ax.imshow(data,norm=simple_norm(data,percent=percentile),cmap=CMAP,origin='lower')
	ax.scatter(xlist,ylist,edgecolor='white',marker='o',s=220,facecolors='none',lw=2,alpha=0.95)

	for i, txt in enumerate(redshift_guesses):
		ax.annotate('%.2f'%txt, (xlist[i], ylist[i]),horizontalalignment='left',fontsize=16,color='white',xytext=(-9, -20), textcoords='offset points')

	ax.set_ylabel(r'DEC [px]',fontsize=22)
	ax.set_xlabel(r'RA [px]',fontsize=22)
	pyplot.tight_layout()
	pyplot.savefig(filename)


def create_globalsky(fits_file_list,masking=False,name_to_replace_with_mask='emline',visit=1,module='a'):
	SKY=[]
	if masking==False:
		savename='globalsky_visit%s_mod%s_unmasked.fits'%(visit,module)
	else:
		savename='globalsky_visit%s_mod%s_masked.fits'%(visit,module)

	for fits_file in fits_file_list:
		hdul = fits.open(fits_file)
		imheader = hdul[1].header
		sci_img = np.copy(hdul[1].data)
		err_img = np.copy(hdul[2].data)

		n_y, n_x = sci_img.shape

		sci_img_original = np.copy(sci_img)
		sci_img_bg = np.copy(sci_img)

		idx_negative_err = np.where(err_img <= 0)
		sci_img_bg[idx_negative_err]=np.nan

		if masking==True:
			maskfile=fits_file.replace(name_to_replace_with_mask,'mask')
			mask = fits.getdata(maskfile, 0)
			idx_maskedout = np.where(mask==1)
			sci_img[idx_maskedout]=np.nan
			err_img[idx_maskedout]=np.nan
			sci_img_bg[idx_maskedout]=np.nan

		err_ptiles = np.nanpercentile(err_img[:,:], [0.15, 1., 2.5, 16., 50., 84., 97.5, 99., 99.85])
		sci_ptiles = np.nanpercentile(sci_img[:,:], [0.15, 1., 2.5, 16., 50., 84., 97.5, 99., 99.85])
		print('ERR Ptiles: ', err_ptiles)
		print('SCI Ptiles: ', sci_ptiles)

		idx_bad = np.where((sci_img<sci_ptiles[1])|(sci_img>sci_ptiles[7])|
		                   (err_img<err_ptiles[1])|(err_img>err_ptiles[7]))

		sci_img_bg[idx_bad]=np.nan
		SKY.append(sci_img_bg)
	globalsky=np.nanmedian(SKY,axis=0)
	fits.writeto(savename,globalsky,overwrite=True)
	print('Global sky created: ',savename)

def filter_times_linear_modB(dummyvar,x_end,norm,slope):
	C = grismconf.Config('../simulate_mirage/MIRAGE_DATA/mirage_data/nircam/GRISM_NIRCAM/current/NIRCAM_F356W_modB_R.conf')
	n_x=2048
	dLdx=C.DDISPL('+1',1000,1000,0.5)/C.DDISPX('+1',1000,1000,0.5) #in the current V3 config this doesnt depend on the position, unit is micron. V4 it changes but this is some sort of average
	lam=np.arange(3.06,4.08,abs(dLdx))
	n_curve=len(lam)

	lam=lam[::-1]
	s = C.SENS["+1"](lam)
	s=s/np.max(s)
	x_array=np.arange(0,n_x,1)

	#Start of Variable part:
	#print(x_end,norm,slope)
	x_start=int(x_end)-n_curve

	x_dummy=np.append(np.arange(0,x_start,1),np.arange(x_start,int(x_end),1))
	x_dummy=np.append(x_dummy,np.arange(int(x_end),n_x,1))

	y_dummy=np.zeros(len(x_dummy))

	x_curve=np.arange(0.,n_curve,1.)

	y_dummy[(x_dummy>=x_start)*(x_dummy<int(x_end))]+=norm*s*(1+slope*(x_curve-n_curve/2.))
	ykeep=y_dummy[(x_dummy>=0)*(x_dummy<n_x)]
	return ykeep

def filter_times_linear_modA(dummyvar,x_end,norm,slope):
	C = grismconf.Config('../simulate_mirage/MIRAGE_DATA/mirage_data/nircam/GRISM_NIRCAM/current/NIRCAM_F356W_modA_R.conf')
	n_x=2048
	dLdx=C.DDISPL('+1',1000,1000,0.5)/C.DDISPX('+1',1000,1000,0.5) #in the current V3 config this doesnt depend on the position, unit is micron. V4 it changes but this is some sort of average
	lam=np.arange(3.06,4.08,abs(dLdx))
	n_curve=len(lam)

	s = C.SENS["+1"](lam)
	s=s/np.max(s)
	x_array=np.arange(0,n_x,1)

	#Start of Variable part:
	#print(x_end,norm,slope)
	x_start=int(x_end)-n_curve

	x_dummy=np.append(np.arange(0,x_start,1),np.arange(x_start,int(x_end),1))
	x_dummy=np.append(x_dummy,np.arange(int(x_end),n_x,1))

	y_dummy=np.zeros(len(x_dummy))

	x_curve=np.arange(0.,n_curve,1.)

	y_dummy[(x_dummy>=x_start)*(x_dummy<int(x_end))]+=norm*s*(1+slope*(x_curve-n_curve/2.))

	ykeep=y_dummy[(x_dummy>=0)*(x_dummy<n_x)]
	return ykeep

def continuum_finding(rate_file,output_dir,module='a',configfile='../simulate_mirage/MIRAGE_DATA/mirage_data/nircam/GRISM_NIRCAM/current/NIRCAM_F356W_modA_R.conf'):
	print('Running Continuum Finding for file',rate_file)
	updated_filename=rate_file.replace('rate','rate_uhead')
	cal_filename=updated_filename.replace('rate_uhead','rate_uhead_cal')
	globalMedSubt_filename=cal_filename.replace('rate_uhead_cal','globalMedSubt')	

	filter_continuum_filename=cal_filename.replace('rate_uhead_cal','filter_continuum')	
	filter_continuum_res_filename=cal_filename.replace('rate_uhead_cal','filter_continuum_residual')	

	wavelength_filename=cal_filename.replace('rate_uhead_cal','continuum_wavs')	

	#Load the sensitivity curve (first order)

	C = grismconf.Config(configfile)
	dLdx=C.DDISPL('+1',1000,1000,0.5)/C.DDISPX('+1',1000,1000,0.5) #in the current V3 config this doesnt depend on the position, unit is micron. V4 it changes but this is some sort of average

	print(globalMedSubt_filename,dLdx)
	lam=np.arange(3.06,4.08,abs(dLdx))
	n_curve=len(lam)
	print(len(lam))

	s = C.SENS["+1"](lam)
	s=s/np.max(s)



	hdul_rate = fits.open(rate_file) 
	dq_img=np.copy(hdul_rate['DQ'].data)

	hdul = fits.open(globalMedSubt_filename) 
	imheader = hdul[1].header
	sci_img = np.copy(hdul[1].data)



	sel_mask=dq_img>0.5
	#sci_img[sel_mask]=np.nan
	n_y, n_x = sci_img.shape[0], sci_img.shape[1]

	print(n_x)
	x_array=np.arange(0,n_x,1)

	newimage=np.zeros(np.shape(sci_img))
	chisqimage=np.zeros(np.shape(sci_img))
	wavimage=np.zeros(np.shape(sci_img))

	for i_y in range(n_y): 
		thisrow=sci_img[i_y,:]
		sel=~np.isfinite(thisrow)
		thisrow[sel]=0.#np.nan#np.nanmedian(thisrow[~sel])

		sel_dq=(dq_img[i_y,:]>0)*(dq_img[i_y,:]<1E7) #dont change the borders

		#remove outliers
		med,std=np.nanmedian(thisrow),np.nanstd(thisrow)
		thisrow[thisrow-med > 5*std]=np.nan

		maxrow=np.nanmax(thisrow[~sel_dq])

		sel_thresh=thisrow>0.5
		if len(thisrow[sel_thresh])<10:
			#print('ignore',maxrow)
			continue
		print(i_y,maxrow)

		thisrow[sel_dq]=np.nan
		nans, x= nan_helper(thisrow)
		thisrow[nans]= np.interp(x(nans), x(~nans), thisrow[~nans])

		if module=='a':
			model=Model(filter_times_linear_modA) #filter_times_linear(dummyvar,x_end,norm,slope)
		if module=='b':
			model=Model(filter_times_linear_modB) #filter_times_linear(dummyvar,x_end,norm,slope)

		model.set_param_hint('x_end',min=1,max=2048+1016+1)
		#model.set_param_hint('x_end',min=1770,max=1800)

		model.set_param_hint('norm',min=0.,max=10.)
		model.set_param_hint('slope',min=-0.0027,max=0.0027)

		params = model.make_params(x_end=1000.,norm=0.3,slope=0.) #a=12., x0=12, asym=0.25, d=80.,a1=12.,peaksep=200.,asym1=0.25,d1=80
		#params['norm'].vary=False
		#params['slope'].vary=False
		a=time.time()

		result = model.fit(thisrow, params,dummyvar=thisrow,nan_policy='raise',method='differential_evolution')
		b=time.time()
		chisq=result.chisqr
		#print(result.params['slope'].value)
		print(rate_file,'TIME USED',b-a)
		print(result.fit_report())
		pyplot.plot(x_array,thisrow)
		pyplot.plot(x_array,model.eval(result.params,dummyvar=x_array),color='k',lw=3,ls='--')
		pyplot.savefig('fitting_%s.png'%i_y)
		pyplot.clf()
		newimage[i_y,:]=model.eval(result.params,dummyvar=x_array)
		chisqimage[i_y,:]+=chisq

		x_end=result.params['x_end'].value
		x_start=int(x_end)-n_curve
		x_dummy=np.append(np.arange(0,x_start,1),np.arange(x_start,int(x_end),1))
		x_dummy=np.append(x_dummy,np.arange(int(x_end),n_x,1))

		y_dummy=np.zeros(len(x_dummy))
		x_curve=np.arange(0.,n_curve,1.)

		y_dummy[(x_dummy>=x_start)*(x_dummy<int(x_end))]+=lam
		ykeep=y_dummy[(x_dummy>=0)*(x_dummy<n_x)]
		wavimage[i_y,:]=ykeep
	#pyplot.show()
	
	hdu_tmp = fits.ImageHDU(data = newimage, header = imheader)
	hdu_tmp.writeto(filter_continuum_filename, overwrite=True)
	# hdu_tmp = fits.ImageHDU(data = chisqimage, header = imheader)
	# hdu_tmp.writeto('test_chisq.fits', overwrite=True)

	hdu_tmp = fits.ImageHDU(data = sci_img-newimage, header = imheader)
	hdu_tmp.writeto(filter_continuum_res_filename, overwrite=True)
	#hdu_tmp = fits.ImageHDU(data = sci_img, header = imheader)
	#hdu_tmp.writeto('test_contfit_data.fits', overwrite=True)

	hdu_tmp = fits.ImageHDU(data = wavimage, header = imheader)
	hdu_tmp.writeto(wavelength_filename, overwrite=True)



def globalMed_filtering(rate_file,output_dir,run_image2=False,masking=False):
	#FIRST DEFINE SOME FILENAMES:
	updated_filename=rate_file.replace('rate','rate_uhead')
	cal_filename=updated_filename.replace('rate_uhead','rate_uhead_cal')
	globalMed_filename=cal_filename.replace('rate_uhead_cal','globalMed')
	globalMedSubt_filename=cal_filename.replace('rate_uhead_cal','globalMedSubt')

	mask_file=cal_filename.replace('rate_uhead_cal','mask')

	#STEP 1: Modify header such that a grism image can be treated as an image
	hdul = fits.open(rate_file)
	module=hdul[0].header['MODULE']
	print('PUPIL   : ', hdul[0].header['pupil'])
	print('EXP_TYPE: ', hdul[0].header['exp_type'])

	print('Modigy header...')
	hdul[0].header['pupil'] = 'CLEAR'
	hdul[0].header['exp_type'] = 'NRC_IMAGE'
	print('PUPIL   : ', hdul[0].header['pupil'])
	print('EXP_TYPE: ', hdul[0].header['exp_type'])

	hdul.writeto(updated_filename, overwrite=True)
	hdul.close()	

	#RUN IMAGE2
	if run_image2==True:
		image2 = calwebb_image2.Image2Pipeline()
		image2.output_dir = output_dir
		image2.save_results = True
		###image2.photom.skip=True #JM new 23 Jun ## REMOVED AGAIN
		image2.run(updated_filename)

	hdul = fits.open(cal_filename) #### JM 31 march
	imheader = hdul[1].header
	sci_img = np.copy(hdul[1].data)
	err_img = np.copy(hdul[2].data)
	dq_img = np.copy(hdul['DQ'].data)

	n_y, n_x = sci_img.shape[0], sci_img.shape[1]
	sci_img_original = np.copy(sci_img)
	sci_img_bg = np.copy(sci_img)

	if masking==True:
		print('Read mask_fil: ', mask_file)
		mask = fits.getdata(mask_file,0)
		idx_maskedout = np.where(mask>0)  ## Pixels of bright emission lines
		idx_unmasked = np.where(mask==0)   
		print('# idx_maskedout: ', idx_maskedout[0].size)
		print('# idx_unmasked : ', idx_unmasked[0].size)
		sci_img_bg[idx_maskedout]=np.nan
	else:
		idx_unmasked = np.where(np.isfinite(sci_img_bg))
	idx_negative_err = np.where(err_img <= 0)
	print('Pixels negative error: ', idx_negative_err[0].size, '(', updated_filename, ')', flush=True)
	sci_img_bg[idx_negative_err]=np.nan

	err_ptiles = np.nanpercentile(err_img[idx_unmasked], [0.15, 1., 2.5, 16., 50., 84., 97.5, 99., 99.85, 99.99])
	sci_ptiles = np.nanpercentile(sci_img[idx_unmasked], [0.15, 1., 2.5, 16., 50., 84., 97.5, 99., 99.85, 99.99])
	print('ERR Ptiles: ', err_ptiles)
	print('SCI Ptiles: ', sci_ptiles)
	# Bad pixels if outside the 1-99th percentile
	idx_bad = np.where((sci_img<sci_ptiles[1])|(sci_img>sci_ptiles[7])|
	                   (err_img<err_ptiles[1])|(err_img>err_ptiles[7]))

	sci_img_bg[idx_bad]=np.nan

	sci_img_medians_at_x = np.zeros((n_y, n_x))  # Medians at each "X", thus "vertically" striped
	sci_img_medians_at_y = np.zeros((n_y, n_x))  # Medians at each "Y", thus "horizontally" striped

	print('start global xy median')
	sci_img_bg_tmp = sci_img_bg - sci_img_medians_at_y

	for i_x in range(n_x):
	    sci_img_medians_at_x[:,i_x] = np.nanmedian(sci_img_bg_tmp[:,i_x])


	sci_img_cleaned = sci_img_original - sci_img_medians_at_x # - sci_img_medians_at_y
	medians_img =  sci_img_medians_at_x # + sci_img_medians_at_y



	hdu_tmp = fits.ImageHDU(data = sci_img_cleaned, header = imheader)
	hdu_tmp.writeto(globalMedSubt_filename, overwrite=True)

	hdul[1].data = medians_img
	hdul[1].writeto(globalMed_filename, overwrite=True)



def kernel_filtering(rate_file,output_dir,masking=False,cont_filt=False,kx=51,ky=1,kx_gap=9,apply_masterbias=False,masterbias='eiger_masterbias_F356W_GRISMR_nrcalong.fits'):
	#kx ky and kx_gap must be a odd numbers
	#FIRST DEFINE SOME FILENAMES:
	updated_filename=rate_file.replace('rate','rate_uhead')
	cal_filename=updated_filename.replace('rate_uhead','rate_uhead_cal')
	globalMed_filename=cal_filename.replace('rate_uhead_cal','globalMed')
	globalMedSubt_filename=cal_filename.replace('rate_uhead_cal','globalMedSubt')

	filter_continuum_filename=cal_filename.replace('rate_uhead_cal','filter_continuum')	
	filter_continuum_res_filename=cal_filename.replace('rate_uhead_cal','filter_continuum_residual')	

	emlines_file=cal_filename.replace('rate_uhead_cal','emline')
	continua_file=cal_filename.replace('rate_uhead_cal','continua')
	flatfield_file=cal_filename.replace('rate_uhead_cal','flatfieldstep')
	emlines_oldhead_file=cal_filename.replace('rate_uhead_cal','emline_oldhead')
	continua_oldhead_file=cal_filename.replace('rate_uhead_cal','continua_oldhead')

	mask_file=cal_filename.replace('rate_uhead_cal','mask')

	#hdul = fits.open(globalMedSubt_filename) #### JM 31 march

	hdul = fits.open(cal_filename) #### JM 31 march


	if apply_masterbias==True:
		emlines_file=cal_filename.replace('rate_uhead_cal','emline_bsub')



	#ADDING THIS 2 MAY, I May need to add more
	hdul2 = fits.open(cal_filename) #### JM 2  may
	sci_img = np.copy(hdul2[1].data)

	err_img = np.copy(hdul2[2].data)
	dq_img = np.copy(hdul2['DQ'].data)


	n_y, n_x = sci_img.shape[0], sci_img.shape[1]
	sci_img_original = np.copy(sci_img)
	sci_img_bg = np.copy(sci_img)

	if masking==True:
		print('Read mask_fil: ', mask_file)
		mask = fits.getdata(mask_file,0)
		idx_maskedout = np.where(mask>0)  ## Pixels of bright emission lines
		idx_unmasked = np.where(mask==0)   
		print('# idx_maskedout: ', idx_maskedout[0].size)
		print('# idx_unmasked : ', idx_unmasked[0].size)
		sci_img_bg[idx_maskedout]=np.nan
	else:
		idx_unmasked = np.where(np.isfinite(sci_img_bg))
	idx_negative_err = np.where(err_img <= 0)
	print('Pixels negative error: ', idx_negative_err[0].size, '(', updated_filename, ')', flush=True)
	sci_img_bg[idx_negative_err]=np.nan

	err_ptiles = np.nanpercentile(err_img[idx_unmasked], [0.15, 1., 2.5, 16., 50., 84., 97.5, 99., 99.85, 99.99])
	sci_ptiles = np.nanpercentile(sci_img[idx_unmasked], [0.15, 1., 2.5, 16., 50., 84., 97.5, 99., 99.85, 99.99])
	# Bad pixels if outside the 1-99th percentile
	idx_bad = np.where((sci_img<sci_ptiles[1])|(sci_img>sci_ptiles[7])|
	                   (err_img<err_ptiles[1])|(err_img>err_ptiles[7]))

	sci_img_bg[idx_bad]=np.nan

	sci_img_medians_at_x = np.zeros((n_y, n_x))  # Medians at each "X", thus "vertically" striped
	sci_img_medians_at_y = np.zeros((n_y, n_x))  # Medians at each "Y", thus "horizontally" striped
	#END ADDING 2 May









	if cont_filt==False:
		sci_img_new_bg = fits.getdata(globalMedSubt_filename)
	if cont_filt==True:
		sci_img_new_bg = fits.getdata(filter_continuum_res_filename)
	sci_img_new_bg[idx_negative_err]=np.nan


	if apply_masterbias==True:
		masterbiasdata=fits.getdata(masterbias)
		sci_img_new_bg=sci_img_new_bg-masterbiasdata

	# Masking
	if masking==True:
		sci_img_new_bg[idx_maskedout]=np.nan
	#####sci_img_new_bg[idx_maskedout]=np.nan ##not done yet

	idx_bad_new = np.where((sci_img<sci_ptiles[1])|(sci_img>sci_ptiles[9])|
	                       (err_img<err_ptiles[1])|(err_img>err_ptiles[9]))
	sci_img_new_bg[idx_bad_new] = np.nan

	# Empty array for continua image
	continua_img = np.zeros((n_y, n_x))

	# Set the kernel
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
	    med_tmp = np.nanmedian(sci_img_new_bg[idx_y[iy,0]:idx_y[iy,-1]+1,idx_x], axis=[0,2])
	    continua_img[iy,:] = med_tmp[:]
	print( ': Loop end (', flush=True)

	# Emission-line image (continuum extracted)
	if cont_filt==True:
		sci_img_emlines = fits.getdata(filter_continuum_res_filename) - continua_img
		continua_img+=fits.getdata(filter_continuum_filename)
	else:
		sci_img_emlines = fits.getdata(globalMedSubt_filename) - continua_img ##check this now

	# Save files
	hdul[1].data = sci_img_emlines
	hdul.writeto(emlines_file, overwrite=True)
	print( ': Saved: ', emlines_file)

	hdul[1].data = continua_img
	hdul.writeto(continua_file, overwrite=True)
	print( ': Saved: ', continua_file)
	hdul.close()

	#this is added by JM #these are using the headers from flatfield_file because those have a WCS assigned
	hdul = fits.open(flatfield_file) 
	hdul[1].data = sci_img_emlines
	hdul.writeto(emlines_oldhead_file, overwrite=True)
	print( ': Saved: ', emlines_oldhead_file)
	hdul.close()

	hdul = fits.open(flatfield_file)
	hdul[1].data = continua_img
	hdul.writeto(continua_oldhead_file, overwrite=True)
	print( ': Saved: ', continua_oldhead_file)
	hdul.close()




def median_filtering(rate_file,output_dir,run_image2=False,masking=False,kx=51,ky=1,kx_gap=9): ##Depreciated 2 May
	#kx ky and kx_gap must be a odd numbers
	#FIRST DEFINE SOME FILENAMES:
	updated_filename=rate_file.replace('rate','rate_uhead')
	cal_filename=updated_filename.replace('rate_uhead','rate_uhead_cal')
	globalMed_filename=cal_filename.replace('rate_uhead_cal','globalMed')
	globalMedSubt_filename=cal_filename.replace('rate_uhead_cal','globalMedSubt')

	filter_continuum_filename=cal_filename.replace('rate_uhead_cal','filter_continuum')	
	filter_continuum_res_filename=cal_filename.replace('rate_uhead_cal','filter_continuum_residual')	

	emlines_file=cal_filename.replace('rate_uhead_cal','emline')
	continua_file=cal_filename.replace('rate_uhead_cal','continua')
	flatfield_file=cal_filename.replace('rate_uhead_cal','flatfieldstep')
	emlines_oldhead_file=cal_filename.replace('rate_uhead_cal','emline_oldhead')
	continua_oldhead_file=cal_filename.replace('rate_uhead_cal','continua_oldhead')

	mask_file=cal_filename.replace('rate_uhead_cal','mask')

	#STEP 1: Modify header such that a grism image can be treated as an image
	hdul = fits.open(rate_file)
	print('PUPIL   : ', hdul[0].header['pupil'])
	print('EXP_TYPE: ', hdul[0].header['exp_type'])

	print('Modigy header...')
	hdul[0].header['pupil'] = 'CLEAR'
	hdul[0].header['exp_type'] = 'NRC_IMAGE'
	print('PUPIL   : ', hdul[0].header['pupil'])
	print('EXP_TYPE: ', hdul[0].header['exp_type'])

	hdul.writeto(updated_filename, overwrite=True)
	hdul.close()	

	#RUN IMAGE2
	if run_image2==True:
		image2 = calwebb_image2.Image2Pipeline()
		image2.output_dir = output_dir
		image2.save_results = True
		image2.run(updated_filename)

	hdul = fits.open(cal_filename) #### JM 31 march
	imheader = hdul[1].header
	sci_img = np.copy(hdul[1].data)
	err_img = np.copy(hdul[2].data)
	dq_img = np.copy(hdul['DQ'].data)

	n_y, n_x = sci_img.shape[0], sci_img.shape[1]
	sci_img_original = np.copy(sci_img)
	sci_img_bg = np.copy(sci_img)

	if masking==True:
		print('Read mask_fil: ', mask_file)
		mask = fits.getdata(mask_file,0)
		idx_maskedout = np.where(mask>0)  ## Pixels of bright emission lines
		idx_unmasked = np.where(mask==0)   
		print('# idx_maskedout: ', idx_maskedout[0].size)
		print('# idx_unmasked : ', idx_unmasked[0].size)
		sci_img_bg[idx_maskedout]=np.nan
	else:
		idx_unmasked = np.where(np.isfinite(sci_img_bg))
	idx_negative_err = np.where(err_img <= 0)
	print('Pixels negative error: ', idx_negative_err[0].size, '(', updated_filename, ')', flush=True)
	sci_img_bg[idx_negative_err]=np.nan

	err_ptiles = np.nanpercentile(err_img[idx_unmasked], [0.15, 1., 2.5, 16., 50., 84., 97.5, 99., 99.85, 99.99])
	sci_ptiles = np.nanpercentile(sci_img[idx_unmasked], [0.15, 1., 2.5, 16., 50., 84., 97.5, 99., 99.85, 99.99])
	print('ERR Ptiles: ', err_ptiles)
	print('SCI Ptiles: ', sci_ptiles)
	# Bad pixels if outside the 1-99th percentile
	idx_bad = np.where((sci_img<sci_ptiles[1])|(sci_img>sci_ptiles[7])|
	                   (err_img<err_ptiles[1])|(err_img>err_ptiles[7]))

	sci_img_bg[idx_bad]=np.nan

	sci_img_medians_at_x = np.zeros((n_y, n_x))  # Medians at each "X", thus "vertically" striped
	sci_img_medians_at_y = np.zeros((n_y, n_x))  # Medians at each "Y", thus "horizontally" striped

	print('start global xy median')
	sci_img_bg_tmp = sci_img_bg - sci_img_medians_at_y

	for i_x in range(n_x):
	    sci_img_medians_at_x[:,i_x] = np.nanmedian(sci_img_bg_tmp[:,i_x])


	sci_img_cleaned = sci_img_original - sci_img_medians_at_x # - sci_img_medians_at_y
	medians_img =  sci_img_medians_at_x # + sci_img_medians_at_y



	hdu_tmp = fits.ImageHDU(data = sci_img_cleaned, header = imheader)
	hdu_tmp.writeto(globalMedSubt_filename, overwrite=True)

	hdul[1].data = medians_img
	hdul[1].writeto(globalMed_filename, overwrite=True)

	#SECONDARY CONTINUUM SUBTRACTION
	#sci_img_new_bg = np.copy(sci_img_cleaned)  #JM 28 April

	sci_img_new_bg = fits.getdata(filter_continuum_res_filename)
	sci_img_new_bg[idx_negative_err]=np.nan

	# Masking
	if masking==True:
		sci_img_new_bg[idx_maskedout]=np.nan
	#####sci_img_new_bg[idx_maskedout]=np.nan ##not done yet

	idx_bad_new = np.where((sci_img<sci_ptiles[1])|(sci_img>sci_ptiles[9])|
	                       (err_img<err_ptiles[1])|(err_img>err_ptiles[9]))
	sci_img_new_bg[idx_bad_new] = np.nan

	# Empty array for continua image
	continua_img = np.zeros((n_y, n_x))

	# Set the kernel
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
	    med_tmp = np.nanmedian(sci_img_new_bg[idx_y[iy,0]:idx_y[iy,-1]+1,idx_x], axis=[0,2])
	    continua_img[iy,:] = med_tmp[:]
	print( ': Loop end (', flush=True)

	# Emission-line image (continuum extracted)
	sci_img_emlines = fits.getdata(filter_continuum_res_filename) - continua_img

	continua_img+=fits.getdata(filter_continuum_filename)

	# Save files
	hdul[1].data = sci_img_emlines
	hdul.writeto(emlines_file, overwrite=True)
	print( ': Saved: ', emlines_file)

	hdul[1].data = continua_img
	hdul[1].writeto(continua_file, overwrite=True)
	print( ': Saved: ', continua_file)
	hdul.close()

	#this is added by JM #these are using the headers from flatfield_file because those have a WCS assigned
	hdul = fits.open(flatfield_file) 
	hdul[1].data = sci_img_emlines
	hdul.writeto(emlines_oldhead_file, overwrite=True)
	print( ': Saved: ', emlines_oldhead_file)
	hdul.close()

	hdul = fits.open(flatfield_file)
	hdul[1].data = continua_img
	hdul.writeto(continua_oldhead_file, overwrite=True)
	print( ': Saved: ', continua_oldhead_file)
	hdul.close()


def create_emission_line_masks(stacked_emission_image_file,filelist_maskrequired,name_to_replace_with_mask='rate_uhead_cal'):
	#requires a swarp installation and a config.swarp file
	hdu= fits.open(stacked_emission_image_file)
	hd=hdu['SCI'].header
	data=hdu['SCI'].data
	data_err=hdu['ERR'].data**0.5

	mean_noise=np.nanmedian(data_err)
	std_noise=np.nanstd(data_err)
	print(std_noise)

	sel_use=(data_err>mean_noise-0.5*std_noise)*(data_err<mean_noise+0.5*std_noise)

	sel_traces=data_err>mean_noise+0.8*std_noise


	rescale=np.nanstd(data[sel_use]) / np.nanmedian(data_err[sel_use])
	data_err=data_err*rescale


	data_err[sel_traces]=data_err[sel_traces]*10

	fits.writeto('detection_image.fits',data,hd,overwrite=True)
	fits.writeto('weight.fits',data_err,hd,overwrite=True)
	hdu.close()
	os.system('source-extractor detection_image.fits -c detect_wfss.sex')

	#create a mask for each rate file:
	for cal_filename in filelist_maskrequired:
		print('mask for',cal_filename)
		maskname=cal_filename.replace(name_to_replace_with_mask,'mask')
		hdul = fits.open(cal_filename)
		imheader = hdul[1].header
		wcs=WCS(imheader) #see above

		hdu2=fits.open('check.fits')
		wcs_mask=WCS(hdu2['SCI'].header)
		mask_data=hdu2['SCI'].data
		mask=np.zeros((2048,2048))
		for j in range(2038):
				#ra_pix,dec_pix=wcs.all_pix2world(j,k,0)
				#x_mask,y_mask=wcs_mask.all_world2pix(ra_pix,dec_pix,0)
				ra_pix,dec_pix=wcs.all_pix2world(np.arange(2038)+5,np.zeros(2038)+j+5,0) ##Changed 2048 to 2038 to not mask borders
				
				x_mask,y_mask=wcs_mask.all_world2pix(ra_pix,dec_pix,0)
				#print(np.round(x_mask))
				#print(mask_data[np.array(np.round(x_mask),dtype='int'),0])

				#print(x_mask,np.shape(mask_data))
				mask[j+5,5:-5]=mask_data[np.array(np.round(y_mask),dtype='int'),np.array(np.round(x_mask),dtype='int')]


		fits.writeto(maskname,mask,overwrite=True)
		hdul.close()


def DK_create_emission_line_masks(stacked_emission_image_file,filelist_maskrequired,name_to_replace_with_mask='rate_uhead_cal'):
	#requires a swarp installation and a config.swarp file
	hdu= fits.open(stacked_emission_image_file)
	hd=hdu['SCI'].header
	data=hdu['SCI'].data
        data_wht=hdu['WHT'].data
        data_err=hdu['ERR'].data
        data_errwht = data_err * data_wht**0.5

	mean_noise=np.nanmedian(data_errwht)
	std_noise=np.nanstd(data_errwht)
	print(std_noise)

	#sel_use=(data_err>mean_noise-0.5*std_noise)*(data_err<mean_noise+0.5*std_noise)

	sel_traces = data_errwht > mean_noise+0.8*std_noise

        

	#rescale=np.nanstd(data[sel_use]) / np.nanmedian(data_err[sel_use])
	#data_err=data_err*rescale

	data_err[sel_traces]=data_err[sel_traces]*10

	fits.writeto('detection_image.fits',data,hd,overwrite=True)
	fits.writeto('weight.fits',data_err,hd,overwrite=True)
	hdu.close()
	os.system('source-extractor detection_image.fits -c detect_wfss.sex')

	#create a mask for each rate file:
	for cal_filename in filelist_maskrequired:
		print('mask for',cal_filename)
		maskname=cal_filename.replace(name_to_replace_with_mask,'mask')
		hdul = fits.open(cal_filename)
		imheader = hdul[1].header
		wcs=WCS(imheader) #see above

		hdu2=fits.open('check.fits')
		wcs_mask=WCS(hdu2['SCI'].header)
		mask_data=hdu2['SCI'].data
		mask=np.zeros((2048,2048))
		for j in range(2038):
				#ra_pix,dec_pix=wcs.all_pix2world(j,k,0)
				#x_mask,y_mask=wcs_mask.all_world2pix(ra_pix,dec_pix,0)
				ra_pix,dec_pix=wcs.all_pix2world(np.arange(2038)+5,np.zeros(2038)+j+5,0) ##Changed 2048 to 2038 to not mask borders
				
				x_mask,y_mask=wcs_mask.all_world2pix(ra_pix,dec_pix,0)
				#print(np.round(x_mask))
				#print(mask_data[np.array(np.round(x_mask),dtype='int'),0])

				#print(x_mask,np.shape(mask_data))
				mask[j+5,5:-5]=mask_data[np.array(np.round(y_mask),dtype='int'),np.array(np.round(x_mask),dtype='int')]


		fits.writeto(maskname,mask,overwrite=True)
		hdul.close()



def xy_drizzled_to_xy_rate(emfile,drizzlefile,x_drizzled,y_drizzled):
	#the emfile is the single emission-line fine (from a _cal.fits) corresponding to a rate.fits 
	#the outout are the pixel positions in the original rate file
	hdul = fits.open(emfile)
	wcs_em=WCS(hdul['SCI'].header) 

	hdu2=fits.open(drizzlefile)
	wcs_drizzle=WCS(hdu2['SCI'].header)

	ra_pix,dec_pix=wcs_drizzle.all_pix2world(x_drizzled,y_drizzled,0)
	x_rate,y_rate=wcs_em.all_world2pix(ra_pix,dec_pix,0)

	return x_rate,y_rate
	#The rate file needs to have a WCS assigned.
	#If you want to do this:
	#from jwst.assign_wcs.assign_wcs_step import AssignWcsStep			
	#step = AssignWcsStep()
	#step.save_results = False#True
	#step.output_dir = output_dir
	#result = step.run(rate_file)

def split_list_by_primarydither(filelist):
	#Separate 3 primary dithers
	PD_1=[]
	PD_2=[]
	PD_3=[]
	for thisfile in filelist:
		hd=fits.getheader(thisfile)
		number=hd['PATT_NUM']
		if number==1:
			PD_1.append(thisfile)
		if number==2:
			PD_2.append(thisfile)
		if number==3:
			PD_3.append(thisfile)
	return PD_1,PD_2,PD_3



def continuum_finding_for_pool(rate_file):
	print('Running Continuum Finding for file',rate_file)

	configfile='../simulate_mirage/MIRAGE_DATA/mirage_data/nircam/GRISM_NIRCAM/current/NIRCAM_F356W_modA_R.conf'
	
	module=rate_file.split('nrc')[1][0]

	updated_filename=rate_file.replace('rate','rate_uhead')
	cal_filename=updated_filename.replace('rate_uhead','rate_uhead_cal')
	globalMedSubt_filename=cal_filename.replace('rate_uhead_cal','globalMedSubt')	

	filter_continuum_filename=cal_filename.replace('rate_uhead_cal','filter_continuum')	
	filter_continuum_res_filename=cal_filename.replace('rate_uhead_cal','filter_continuum_residual')	

	wavelength_filename=cal_filename.replace('rate_uhead_cal','continuum_wavs')	

	#Load the sensitivity curve (first order)

	C = grismconf.Config(configfile)
	dLdx=C.DDISPL('+1',1000,1000,0.5)/C.DDISPX('+1',1000,1000,0.5) #in the current V3 config this doesnt depend on the position, unit is micron. V4 it changes but this is some sort of average
	print(globalMedSubt_filename,dLdx)
	lam=np.arange(3.06,4.08,abs(dLdx))
	n_curve=len(lam)
	print(len(lam))

	s = C.SENS["+1"](lam)
	s=s/np.max(s)



	hdul_rate = fits.open(rate_file) 
	dq_img=np.copy(hdul_rate['DQ'].data)

	hdul = fits.open(globalMedSubt_filename) 
	imheader = hdul[1].header
	sci_img = np.copy(hdul[1].data)



	sel_mask=dq_img>0
	#sci_img[sel_mask]=np.nan
	n_y, n_x = sci_img.shape[0], sci_img.shape[1]

	print(n_x)
	x_array=np.arange(0,n_x,1)

	newimage=np.zeros(np.shape(sci_img))
	chisqimage=np.zeros(np.shape(sci_img))
	wavimage=np.zeros(np.shape(sci_img))

	for i_y in range(n_y): 
		thisrow=sci_img[i_y,:]
		sel=~np.isfinite(thisrow)
		thisrow[sel]=0.#np.nan#np.nanmedian(thisrow[~sel])

		sel_dq=(dq_img[i_y,:]>0)*(dq_img[i_y,:]<1E7) #dont change the borders

		med,std=np.nanmedian(thisrow),np.nanstd(thisrow)
		#thisrow[thisrow-med > 5*std]=np.nan

		maxrow=np.nanmax(thisrow[~sel_dq])
		sel_thresh=thisrow>0.08

		if len(thisrow[sel_thresh])<10:
			#print('ignore',maxrow)
			continue
		print(i_y,maxrow)

		thisrow[sel_dq]=np.nan
		nans, x= nan_helper(thisrow)
		thisrow[nans]= np.interp(x(nans), x(~nans), thisrow[~nans])

		if module=='a':
			model=Model(filter_times_linear_modA) #filter_times_linear(dummyvar,x_end,norm,slope)
		if module=='b':
			model=Model(filter_times_linear_modB) #filter_times_linear(dummyvar,x_end,norm,slope)

		model.set_param_hint('x_end',min=1,max=2048+1016+1)
		#model.set_param_hint('x_end',min=1770,max=1800)

		model.set_param_hint('norm',min=0.,max=10.)
		model.set_param_hint('slope',min=-0.0019,max=0.0019)

		params = model.make_params(x_end=1000.,norm=0.3,slope=0.) #a=12., x0=12, asym=0.25, d=80.,a1=12.,peaksep=200.,asym1=0.25,d1=80
		#params['norm'].vary=False
		#params['slope'].vary=False
		a=time.time()

		result = model.fit(thisrow, params,dummyvar=thisrow,nan_policy='raise',method='differential_evolution')
		b=time.time()
		chisq=result.chisqr
		#print(result.params['slope'].value)
		print(rate_file,'TIME USED',b-a)
		#print(result.fit_report())
		#pyplot.plot(x_array,thisrow)
		#pyplot.plot(x_array,model.eval(result.params,dummyvar=x_array),color='k',lw=3,ls='--')

		newimage[i_y,:]=model.eval(result.params,dummyvar=x_array)
		chisqimage[i_y,:]+=chisq

		x_end=result.params['x_end'].value
		x_start=int(x_end)-n_curve
		x_dummy=np.append(np.arange(0,x_start,1),np.arange(x_start,int(x_end),1))
		x_dummy=np.append(x_dummy,np.arange(int(x_end),n_x,1))

		y_dummy=np.zeros(len(x_dummy))
		x_curve=np.arange(0.,n_curve,1.)

		y_dummy[(x_dummy>=x_start)*(x_dummy<int(x_end))]+=lam
		ykeep=y_dummy[(x_dummy>=0)*(x_dummy<n_x)]
		wavimage[i_y,:]=ykeep
	#pyplot.show()
	
	hdu_tmp = fits.ImageHDU(data = newimage, header = imheader)
	hdu_tmp.writeto(filter_continuum_filename, overwrite=True)
	# hdu_tmp = fits.ImageHDU(data = chisqimage, header = imheader)
	# hdu_tmp.writeto('test_chisq.fits', overwrite=True)

	hdu_tmp = fits.ImageHDU(data = sci_img-newimage, header = imheader)
	hdu_tmp.writeto(filter_continuum_res_filename, overwrite=True)
	#hdu_tmp = fits.ImageHDU(data = sci_img, header = imheader)
	#hdu_tmp.writeto('test_contfit_data.fits', overwrite=True)

	hdu_tmp = fits.ImageHDU(data = wavimage, header = imheader)
	hdu_tmp.writeto(wavelength_filename, overwrite=True)


def apply_gain(rate_file,gainfile):
	hdul = fits.open(rate_file)
	gain = fits.open(gainfile)[1].data
	hdul[1].data=hdul[1].data * gain
	hdul[2].data=hdul[2].data * gain
	hdul.writeto(rate_file, overwrite=True)
	hdul.close()	
	
def mask_snowballs(ratefile):
	rate_hdu = fits.open(ratefile)

	dq_data=rate_hdu['DQ'].data
	dq_mask = (build_mask(dq_data,  '~JUMP_DET') == 0)
	segment_map = detect_sources(dq_mask, 0.5, npixels=50)
	cat = SourceCatalog(dq_mask, segment_map)
	tab =  cat.to_table()

	layers = [[50, 0.5, 5], [200, 1.0, 10], [500, 1.0, 20], [1000, 1.0, 30]]
	save_mask =np.zeros([len(layers), dq_mask.shape[0], dq_mask.shape[1]], dtype=bool)
	for il, layer in enumerate(layers):
		seg = segment_map.copy()
		cal_sel = tab[(tab['eccentricity'].value < layer[1])&(tab['area'].value > layer[0])]
		seg.keep_labels(labels=cal_sel['label'].data)
		footprint = circular_footprint(radius=layer[2])
		save_mask[il,:,:] = seg.make_source_mask(footprint=footprint)

	sball_mask=np.any(save_mask, axis=0)
	rate_hdu['DQ'].data[(sball_mask)] = 1
	rate_hdu.writeto(ratefile, overwrite=True)



def create_wfss_stack(emlinefile,contfile):
	hdul = fits.open(emlinefile)
	cont = fits.open(contfile)[1].data
	hdul[1].data=hdul[1].data + cont
	hdul.writeto('temp.fits', overwrite=True)
	hdul.close()
	
def create_masterbias(filelist,masking,mask_brightpixels,masterbias_filename):
	#FIRST DEFINE SOME FILENAMES:
	print('Creating masterflat')
	stack=[]
	for thisfile in filelist:
		mask_file=thisfile.replace('emline','mask')
		
		data=fits.getdata(thisfile)
		hdul = fits.open(thisfile) 
		data = np.copy(hdul[1].data)
		err=np.copy(hdul[2].data)
		dq = np.copy(hdul['DQ'].data)
		#data[dq>1]=np.nan
		
		idx_negative_err = np.where(err <= 0)
		data[idx_negative_err]=np.nan
		
		if mask_brightpixels==True:
			err_ptiles = np.nanpercentile(err[:,:], [0.15, 1., 2.5, 16., 50., 84., 97.5, 99., 99.85])
			sci_ptiles = np.nanpercentile(data[:,:], [0.15, 1., 2.5, 16., 50., 84., 97.5, 99., 99.85])
			idx_bad = np.where((data<sci_ptiles[1])|(data>sci_ptiles[7])| (err<err_ptiles[1])|(err>err_ptiles[7]))
			data[idx_bad]=np.nan

		if masking==True:
			mask=fits.getdata(mask_file)
			idx_maskedout = np.where(mask>0)  ## Pixels of bright emission lines
			data[idx_maskedout]=np.nan

		else:
			stack.append(data)
		hdul.close()
	median_stack=np.nanmedian(stack,axis=0)
	median_stack[np.where(np.isnan(median_stack))]=0.
	fits.writeto(masterbias_filename,median_stack,overwrite=True)


def apply_masterflat(rate_file,masterflat_file):
	hdul = fits.open(rate_file)
	masterflat = fits.getdata(masterflat_file)
	hdul[1].data=hdul[1].data / masterflat
	hdul[2].data=hdul[2].data /masterflat
	hdul.writeto(rate_file, overwrite=True)
	hdul.close()	
	
	

def apply_masterbias(rate_file,masterbias_file):
	new_filename=rate_file.replace('emline','emline_bsub')

	hdul = fits.open(rate_file)
	masterbias = fits.getdata(masterbias_file)
	hdul[1].data=hdul[1].data - masterbias
	hdul.writeto(new_filename, overwrite=True)
	hdul.close()	
	
		




