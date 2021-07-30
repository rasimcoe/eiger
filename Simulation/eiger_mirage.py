import numpy as np
from mirage.catalogs import catalog_generator
from astropy.io import fits
import h5py

def create_pointsource_catalog(target_RA,target_DEC,AB_MAGNITUDES,filename):
	#AB_MAGNITUDES is assumed to be a list with guesses for [F115W,F200W,F356W]
	ptsrc = catalog_generator.PointSourceCatalog(ra=[target_RA], dec=[target_DEC], starting_index=1)
	ptsrc.add_magnitude_column([AB_MAGNITUDES[0]], instrument='nircam', filter_name='f115w', magnitude_system='abmag')
	ptsrc.add_magnitude_column([AB_MAGNITUDES[1]], instrument='nircam', filter_name='f200w', magnitude_system='abmag')
	ptsrc.add_magnitude_column([AB_MAGNITUDES[2]], instrument='nircam', filter_name='f356w', magnitude_system='abmag')
	ptsrc.save(filename)

	print('Pointsource catalog created')


def create_galaxy_catalog(galaxycat,target_RA,target_DEC,MAGLIM,redshift_lim_low,redshift_lim_high,filename):
	fc=fits.open(galaxycat)
	data_cat=fc[1].data 
	redshift=data_cat.field('redshift')
	IDlist=data_cat.field('ID')
	RA=data_cat.field('RA')-np.nanmean(data_cat.field('RA')) + target_RA #re-center the catalog on target RA
	DEC=data_cat.field('DEC')-np.nanmean(data_cat.field('DEC')) + target_DEC #re-center the catalog on target DEC
	redshift=data_cat.field('redshift')
	F115W=8.9-2.5*np.log10(data_cat.field('NRC_F115W_fnu')*1E-9) #all need to be AB magnitudes and catalog is in nJy
	F200W=8.9-2.5*np.log10(data_cat.field('NRC_F200W_fnu')*1E-9)
	F356W=8.9-2.5*np.log10(data_cat.field('NRC_F356W_fnu')*1E-9)

	F300M=8.9-2.5*np.log10(data_cat.field('NRC_F300M_fnu')*1E-9)
	F335M=8.9-2.5*np.log10(data_cat.field('NRC_F335M_fnu')*1E-9)
	F360M=8.9-2.5*np.log10(data_cat.field('NRC_F360M_fnu')*1E-9)
	F410M=8.9-2.5*np.log10(data_cat.field('NRC_F410M_fnu')*1E-9)
	ellip=1.-data_cat.field('axis_ratio')
	radius=data_cat.field('Re_maj')
	sersic=data_cat.field('sersic_n')
	posang=data_cat.field('position_angle')

	#Select which sources to include
	sel=(redshift>redshift_lim_low)*(redshift<redshift_lim_high)
	sel2=(F115W[sel]<MAGLIM) + (F200W[sel]<MAGLIM)  + (F356W[sel]<MAGLIM) 

	print('The number of simulated galaxies is:',len(F115W[sel][sel2]))


	gal = catalog_generator.GalaxyCatalog(ra=RA[sel][sel2], dec=DEC[sel][sel2], ellipticity=ellip[sel][sel2], radius=radius[sel][sel2], sersic_index=sersic[sel][sel2],
	                                      position_angle=posang[sel][sel2], radius_units='arcsec', starting_index=1)
	gal.add_magnitude_column(F115W[sel][sel2], instrument='nircam', filter_name='f115w', magnitude_system='abmag')
	gal.add_magnitude_column(F200W[sel][sel2], instrument='nircam', filter_name='f200w', magnitude_system='abmag')
	gal.add_magnitude_column(F356W[sel][sel2], instrument='nircam', filter_name='f356w', magnitude_system='abmag')
	gal.add_magnitude_column(F300M[sel][sel2], instrument='nircam', filter_name='f300m', magnitude_system='abmag')
	gal.add_magnitude_column(F335M[sel][sel2], instrument='nircam', filter_name='f335m', magnitude_system='abmag')
	gal.add_magnitude_column(F360M[sel][sel2], instrument='nircam', filter_name='f360m', magnitude_system='abmag')
	gal.add_magnitude_column(F410M[sel][sel2], instrument='nircam', filter_name='f410m', magnitude_system='abmag')
	gal.save(filename)
	print('Galaxy catalog created')
	

def get_qso_spectrum(qso_spec,zqso,average_flamb_QSO):
	fqso=fits.open(qso_spec)
	data_qso=fqso[1].data
	l_qso=data_qso.field('restframe_wavelength')
	f_qso=data_qso.field('flambda')

	#redshift the spectrum
	l_qso=1E-4*l_qso*(1+zqso) #wavelength units are Microns
	f_qso=f_qso/(1+zqso) 
	#normalise the spectrum to the F356W magnitude
	sel_wav=(l_qso>3.)*(l_qso<4.)
	f_qso=f_qso*average_flamb_QSO/np.nanmedian(f_qso[sel_wav])
	return l_qso,f_qso
	
	
def create_galaxy_catalog_grism(galaxycat,target_RA,target_DEC,MAGLIM,redshift_lim_low,redshift_lim_high,filename):
	fc=fits.open(galaxycat)
	data_cat=fc[1].data 
	redshift=data_cat.field('redshift')
	IDlist=data_cat.field('ID')
	RA=data_cat.field('RA')-np.nanmean(data_cat.field('RA')) + target_RA
	DEC=data_cat.field('DEC')-np.nanmean(data_cat.field('DEC')) + target_DEC
	redshift=data_cat.field('redshift')
	F115W=8.9-2.5*np.log10(data_cat.field('NRC_F115W_fnu')*1E-9)
	F200W=8.9-2.5*np.log10(data_cat.field('NRC_F200W_fnu')*1E-9)
	F356W=8.9-2.5*np.log10(data_cat.field('NRC_F356W_fnu')*1E-9)

	F300M=8.9-2.5*np.log10(data_cat.field('NRC_F300M_fnu')*1E-9)
	F335M=8.9-2.5*np.log10(data_cat.field('NRC_F335M_fnu')*1E-9)
	F360M=8.9-2.5*np.log10(data_cat.field('NRC_F360M_fnu')*1E-9)
	F410M=8.9-2.5*np.log10(data_cat.field('NRC_F410M_fnu')*1E-9)
	ellip=1.-data_cat.field('axis_ratio')
	radius=data_cat.field('Re_maj')
	sersic=data_cat.field('sersic_n')
	posang=data_cat.field('position_angle')


	sel=(redshift>redshift_lim_low)*(redshift<redshift_lim_high)
	sel2=F356W[sel]<MAGLIM

	#The following is used for correct bookkeeping. The hdf5 files containing the galaxy spectra are split in redshift intervals and we need to correctly select the galaxies from those separate catalogs.
	sel_cat1=(redshift>5.)    
	sel_cat1_2=(F356W[sel_cat1]<MAGLIM) * (redshift[sel_cat1]>redshift_lim_low)*(redshift[sel_cat1]<redshift_lim_high)
	
	

	#this last part     "* (redshift[sel_cat1]>redshift_lim_low)*(redshift[sel_cat1]<redshift_lim_high)"     is to make sure that only galaxies are included within the determined z_low and z_high. As a result some of these sel_cat? arrays may select zero galaxies, but this should run fine

	sel_cat2=(redshift>4.)*(redshift<5.0)
	sel_cat2_2=(F356W[sel_cat2]<MAGLIM) * (redshift[sel_cat2]>redshift_lim_low)*(redshift[sel_cat2]<redshift_lim_high)

	sel_cat3=(redshift>3.)*(redshift<4.0) 
	sel_cat3_2=(F356W[sel_cat3]<MAGLIM)* (redshift[sel_cat3]>redshift_lim_low)*(redshift[sel_cat3]<redshift_lim_high)

	sel_cat4=(redshift>2.)*(redshift<3.0) 
	sel_cat4_2=(F356W[sel_cat4]<MAGLIM)* (redshift[sel_cat4]>redshift_lim_low)*(redshift[sel_cat4]<redshift_lim_high)

	sel_cat5=(redshift>1.5)*(redshift<2.0) 
	sel_cat5_2=(F356W[sel_cat5]<MAGLIM)* (redshift[sel_cat5]>redshift_lim_low)*(redshift[sel_cat5]<redshift_lim_high)

	sel_cat6=(redshift>1.)*(redshift<1.5)
	sel_cat6_2=(F356W[sel_cat6]<MAGLIM) * (redshift[sel_cat6]>redshift_lim_low)*(redshift[sel_cat6]<redshift_lim_high)

	sel_cat7=(redshift>0.2)*(redshift<1.) 
	sel_cat7_2=(F356W[sel_cat7]<MAGLIM)* (redshift[sel_cat7]>redshift_lim_low)*(redshift[sel_cat7]<redshift_lim_high)

	sel_spectral_cats=[sel_cat1_2,sel_cat2_2,sel_cat3_2,sel_cat4_2,sel_cat5_2,sel_cat6_2,sel_cat7_2]

	print('The number of simulated galaxies is:',len(F115W[sel][sel2]))


	gal = catalog_generator.GalaxyCatalog(ra=RA[sel][sel2], dec=DEC[sel][sel2], ellipticity=ellip[sel][sel2], radius=radius[sel][sel2], sersic_index=sersic[sel][sel2],
	                                      position_angle=posang[sel][sel2], radius_units='arcsec', starting_index=1)
	gal.add_magnitude_column(F115W[sel][sel2], instrument='nircam', filter_name='f115w', magnitude_system='abmag')
	gal.add_magnitude_column(F200W[sel][sel2], instrument='nircam', filter_name='f200w', magnitude_system='abmag')
	gal.add_magnitude_column(F356W[sel][sel2], instrument='nircam', filter_name='f356w', magnitude_system='abmag')
	gal.add_magnitude_column(F300M[sel][sel2], instrument='nircam', filter_name='f300m', magnitude_system='abmag')
	gal.add_magnitude_column(F335M[sel][sel2], instrument='nircam', filter_name='f335m', magnitude_system='abmag')
	gal.add_magnitude_column(F360M[sel][sel2], instrument='nircam', filter_name='f360m', magnitude_system='abmag')
	gal.add_magnitude_column(F410M[sel][sel2], instrument='nircam', filter_name='f410m', magnitude_system='abmag')
	gal.save(filename)
	print('Galaxy catalog created')

	return sel_spectral_cats




def create_catalog_grism_spectra(spectra_cat_root,FLUX_MULTIPLY,sel_spectral_cats,lambda_qso,flux_qso,sed_file):

	wavelength_units = 'microns' 
	flux_units = 'flam'

	#GOING TO LOAD IN ALL THE SPECTRA FROM DIFFERENT CATALOGS. ONLY SELECT THE GALAXIES THAT SATISFY CONSTRAINTS WITH sel_spectral_cats
	f=fits.open(spectra_cat_root+'0p2_1.fits')
	data_lib=f[3].data 
	redshift_lib_7=data_lib.field('redshift')[sel_spectral_cats[6]]
	spec_lib=f[1].data *FLUX_MULTIPLY #Rest-frame erg/s/cm2/A ; so divide by (1+z) as well when redshifting the spectra
	spectra_7=spec_lib[sel_spectral_cats[6],:]


	f=fits.open(spectra_cat_root+'1_1p5.fits')
	data_lib=f[3].data 
	redshift_lib_6=data_lib.field('redshift')[sel_spectral_cats[5]]
	spec_lib=f[1].data *FLUX_MULTIPLY #Rest-frame erg/s/cm2/A ; so divide by (1+z) as well when redshifting the spectra
	spectra_6=spec_lib[sel_spectral_cats[5],:]


	f=fits.open(spectra_cat_root+'1p5_2.fits')
	data_lib=f[3].data 
	redshift_lib_5=data_lib.field('redshift')[sel_spectral_cats[4]]
	spec_lib=f[1].data *FLUX_MULTIPLY #Rest-frame erg/s/cm2/A ; so divide by (1+z) as well when redshifting the spectra
	spectra_5=spec_lib[sel_spectral_cats[4],:]


	f=fits.open(spectra_cat_root+'2_3.fits')
	data_lib=f[3].data 
	redshift_lib_4=data_lib.field('redshift')[sel_spectral_cats[3]]
	spec_lib=f[1].data *FLUX_MULTIPLY #Rest-frame erg/s/cm2/A ; so divide by (1+z) as well when redshifting the spectra
	spectra_4=spec_lib[sel_spectral_cats[3],:]


	f=fits.open(spectra_cat_root+'3_4.fits')
	data_lib=f[3].data 
	redshift_lib_3=data_lib.field('redshift')[sel_spectral_cats[2]]
	spec_lib=f[1].data *FLUX_MULTIPLY #Rest-frame erg/s/cm2/A ; so divide by (1+z) as well when redshifting the spectra
	spectra_3=spec_lib[sel_spectral_cats[2],:]


	f=fits.open(spectra_cat_root+'4_5.fits')
	data_lib=f[3].data 
	redshift_lib_2=data_lib.field('redshift')[sel_spectral_cats[1]]
	spec_lib=f[1].data *FLUX_MULTIPLY #Rest-frame erg/s/cm2/A ; so divide by (1+z) as well when redshifting the spectra
	spectra_2=spec_lib[sel_spectral_cats[1],:]


	f=fits.open(spectra_cat_root+'5_15.fits')
	data_lib=f[3].data 

	redshift_lib_1=data_lib.field('redshift')[sel_spectral_cats[0]]
	spec_lib=f[1].data *FLUX_MULTIPLY #Rest-frame erg/s/cm2/A ; so divide by (1+z) as well when redshifting the spectra  
	spectra_1=spec_lib[sel_spectral_cats[0],:]

	lamb_lib=f[2].data * 1E-4 #Rest-frame wavelength in Angstrom  so multiply by 1E-4 to transofrm to micron  #It's the same for all cats


	#Append them all
	all_spectra=np.append(spectra_7,spectra_6,axis=0)
	all_spectra=np.append(all_spectra,spectra_5,axis=0)
	all_spectra=np.append(all_spectra,spectra_4,axis=0)
	all_spectra=np.append(all_spectra,spectra_3,axis=0)
	all_spectra=np.append(all_spectra,spectra_2,axis=0)
	all_spectra=np.append(all_spectra,spectra_1,axis=0)


	all_redshifts=np.append(redshift_lib_7,redshift_lib_6)
	all_redshifts=np.append(all_redshifts,redshift_lib_5)
	all_redshifts=np.append(all_redshifts,redshift_lib_4)
	all_redshifts=np.append(all_redshifts,redshift_lib_3)
	all_redshifts=np.append(all_redshifts,redshift_lib_2)
	all_redshifts=np.append(all_redshifts,redshift_lib_1)


	print('Number of galaxies in the Spectral Catalog:',np.shape(all_redshifts))



	#NOW GOING TO WRITE THE HDF5 CATALOG
	with h5py.File(sed_file, "w") as file_obj:
		#ADD THE QSO
		sel=(lambda_qso>2.0)*(lambda_qso<4.5) ##We only need these wavelengths for simulating F356W LW Grism data - Otherwise it crashes
		dset = file_obj.create_dataset(str(1), data=[lambda_qso[sel], flux_qso[sel]], dtype='f',
                                       compression="gzip", compression_opts=9)
		dset.attrs[u'wavelength_units'] = wavelength_units
		dset.attrs[u'flux_units'] = flux_units

		#ADD THE GALAXIES
		for i in range(len(all_redshifts)):
			dset = file_obj.create_dataset(str(i+2), data=[lamb_lib*(1+all_redshifts[i]), all_spectra[i]/(1+all_redshifts[i])], dtype='f',
                                       compression="gzip", compression_opts=9)
			dset.attrs[u'wavelength_units'] = wavelength_units
			dset.attrs[u'flux_units'] = flux_units



