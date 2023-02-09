import numpy
import numpy
import matplotlib.pyplot as plt
from jwst import datamodels
from scipy.stats import binned_statistic_2d
import os
from astropy.visualization import simple_norm
import h5py
import numpy as np
from astropy.io import fits
import grismconf
import time
import copy
from mirage.catalogs import catalog_generator

CATALOG='/scratch/EIGER/identification/J1148+5251_photcat_v1_noisemodel_short.fits'
NEW_CATALOG='/scratch/EIGER/identification/J1148+5251_photcat_v1_noisemodel_short_withCandidateDoublets.fits' #SAVE it like this
FOLDER='/scratch/EIGER/identification/SPECTRA_J1148/' #FOLDER WITH SPECTRA

field='J1148'
with fits.open(CATALOG) as hdul:
    orig_table = hdul[1].data
    orig_cols = orig_table.columns
    

IDlist=orig_table.field('NUMBER')
Nlines=[]
YPOS=[]
YPOS_STD=[]
ELLIP=[]
N_O3_CANDIDATES=[]
REDSHIFT_O3_CANDIDATES=[]
SN_O3_BLUE=[]






BLUELIST=[4960.295]#,4862.69,9071.1]
REDLIST=[5008.24]#,5008.24,9533.2]

MINRATIO=[1.5]#,0.3,1.0]
MAXRATIO=[6.0]#,8.0,6.0]
PAIRS=['O3doublet']#,'HbO3','S3doublet']


N=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
SNlist=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
FLUXRATIO=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
EXP_SEP=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
MEAS_SEP=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
REDSHIFT=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
YVAL=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
YVAL_STD=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
LBLUE=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]


"""
#VARIOUS OPTIONS HERE
BLUELIST=[3728,4102.9]
REDLIST=[3870,4341]
MINRATIO=[0.3,1.2]
MAXRATIO=[3,5]

PAIRS=['O2Ne3','HgHd']


N=[np.zeros(len(IDlist)),np.zeros(len(IDlist))]#,np.zeros(len(IDlist))]
SNlist=[np.zeros(len(IDlist)),np.zeros(len(IDlist))]#,np.zeros(len(IDlist))]
FLUXRATIO=[np.zeros(len(IDlist)),np.zeros(len(IDlist))]#,np.zeros(len(IDlist))]
EXP_SEP=[np.zeros(len(IDlist)),np.zeros(len(IDlist))]#,np.zeros(len(IDlist))]
MEAS_SEP=[np.zeros(len(IDlist)),np.zeros(len(IDlist))]#,np.zeros(len(IDlist))]
REDSHIFT=[np.zeros(len(IDlist)),np.zeros(len(IDlist))]#,np.zeros(len(IDlist))]
YVAL=[np.zeros(len(IDlist)),np.zeros(len(IDlist))]#,np.zeros(len(IDlist))]
YVAL_STD=[np.zeros(len(IDlist)),np.zeros(len(IDlist))]#,np.zeros(len(IDlist))]
LBLUE=[np.zeros(len(IDlist)),np.zeros(len(IDlist))]#,np.zeros(len(IDlist))]




BLUELIST=[4960.295,4862.69,9071.1]
REDLIST=[5008.24,5008.24,9533.2]

MINRATIO=[1.5,0.3,1.0]
MAXRATIO=[6.0,8.0,6.0]
PAIRS=['O3doublet','HbO3','S3doublet']

N=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
SNlist=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
FLUXRATIO=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
EXP_SEP=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
MEAS_SEP=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
REDSHIFT=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
YVAL=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
YVAL_STD=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
LBLUE=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]


BLUELIST=[4960,4860,6564,6564,3728,9069,10830]
REDLIST=[5008,5008,6585,6718,3870,9531,10938]

MINRATIO=[1.5,0.3,0.0001,0.0001,0.5,1.0,0.2]
MAXRATIO=[6.0,8.0,1000,1000,10,6.0,10]
#Experiment with deblending, NII

PAIRS=['O3doublet','HbO3','HaN2','HaS2','O2Ne3','S3doublet','PagHeI']



#ZZZZ
N=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
SNlist=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
FLUXRATIO=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
EXP_SEP=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
MEAS_SEP=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
REDSHIFT=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
YVAL=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
YVAL_STD=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]
LBLUE=[np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist)),np.zeros(len(IDlist))]

"""


##NEW
for jjj in range(100):#range(len(IDlist)):
	thisID=IDlist[jjj]+1000
	#thisID=17242
	print(thisID)
	try:
		file=FOLDER+'stacked_2D_%s_%s.fits'%(field,thisID)

		hdu= fits.open(file)
		hd=hdu['EMLINE'].header
		data=hdu['EMLINE'].data

		#RESCALING NOISE iteratively masking residuals and lines
		standard=numpy.nanstd(data[:,200:980])
		errdata=hdu['ERR'].data**0.5  #THIS IS BECAUSE THE ERR extension is actually VARIANCE -- needs to be fixed

		#print('standard deviation of emline data before modification',standard)

		use_err=copy.deepcopy(errdata[:,200:980])
		sel_ignore=use_err==0.
		use_err[sel_ignore]=numpy.nan

		use_dat=copy.deepcopy(data[:,200:980])
		sel_ignore=numpy.abs(use_dat)>5*standard
		use_dat[sel_ignore]=numpy.nan
		standard=numpy.nanstd(use_dat)
		sel_ignore=numpy.abs(use_dat)>3*standard
		use_dat[sel_ignore]=numpy.nan
		standard=numpy.nanstd(use_dat)		
		sel_ignore=numpy.abs(use_dat)>3*standard
		use_dat[sel_ignore]=numpy.nan
		standard=numpy.nanstd(use_dat)	
		sel_ignore=numpy.abs(use_dat)>3*standard
		use_dat[sel_ignore]=numpy.nan
		standard=numpy.nanstd(use_dat)	
		#print('standard deviation of emline data after modification',standard)


		errdata=errdata * standard/numpy.nanmedian(use_err) #Renormalising
		fits.writeto('weight_rms.fits',errdata,hd,overwrite=True)

		fits.writeto('2D_detection_image.fits',data,hd,overwrite=True)

		os.system('source-extractor 2D_detection_image.fits -c sex_O3_j1148.config')

		#Create Catalog
		fc=fits.open('2D_detection_catalog.fits')
		data_cat=fc[1].data 

		x=data_cat.field('X_IMAGE')
		y=data_cat.field('Y_IMAGE')
		flux=data_cat.field('FLUX_APER')
		flux_err=data_cat.field('FLUXERR_APER')
		SN=flux/flux_err

		wav=30000+x*9.75 #in angstrom

		sort=np.argsort(wav)
		wav=wav[sort]
		flux=flux[sort]
		y=y[sort]
		SN=SN[sort]
		x=x[sort]

		sel=(y>25-3)*(y<25+3)*(SN>2.)

		#print(x[sel],SN[sel])
		#print(x[sel][np.argsort(x[sel])])

		#print(wav[sel],wav[sel]/4960. -1,wav[sel]/5008. -1, SN[sel])
		###sel=(y>31)*(y<34)*(SN>3)

		if len(x[sel])==0:
			#print(thisID,'The number of detected lines is:',len(x[sel]))
			#No lines here, so skip

			continue


		for thispair in range(len(PAIRS)):
			PAIRNAME=PAIRS[thispair]
			l1=BLUELIST[thispair]
			l2=REDLIST[thispair]
			min_redshift=-1+31500/l1
			max_redshift=-1+39500/l2
			min_separation=(l2-l1)*(1+min_redshift)
			max_separation=(l2-l1)*(1+max_redshift)
			minflux=MINRATIO[thispair]
			maxflux=MAXRATIO[thispair]

			for q in range(len(x[sel])):
				thiswav=wav[sel][q]
				thisx=x[sel][q]
				thisy=y[sel][q]
				thisSN=SN[sel][q]
				thisflux=flux[sel][q]
				distances=wav[sel]-thiswav
				flux_differences=flux[sel]/thisflux

				othery=np.abs(thisy-y[sel])

				#print(flux_differences)
				#print(numpy.nanmin(SN[sel]))

				#sel_candidates=(flux_differences>1.5)*(flux_differences<6)*(distances>0)*(distances>min_separation-20.)*(distances<max_separation+20.)

				#Let's try this:
				expected_z=-1+thiswav/l1
				expected_separation=(l2-l1)*(1+expected_z)

				#print(thisID,PAIRNAME,distances,expected_separation,'fluxdif',flux_differences)

				sel_candidates=(flux_differences>minflux)*(flux_differences<maxflux)*(distances>0)*(distances>expected_separation-20.)*(distances<expected_separation+20.)*(othery<2)


				if len(distances[sel_candidates])>0.1:
					Z_thispair=np.nanmean(-1+thiswav/l1)
					SN_thispair=np.nanmin(SN[sel][sel_candidates])
					print('SN blue,red',thisSN,SN_thispair)
					FLUXRATIO_thispair=np.nanmin(flux_differences[sel_candidates])
					meas_dist=np.nanmin(distances[sel_candidates])
					stdy=np.nanmean(othery[sel_candidates])

					print('ID',thisID,'Found candidate pair!',PAIRNAME,'redshift',Z_thispair,'Fluxratio',FLUXRATIO_thispair,'Sep (Exp)',meas_dist,'(%s)'%expected_separation,'Y',thisy,'stdY',stdy)
					N[thispair][jjj]+=1
					SNlist[thispair][jjj]=SN_thispair
					REDSHIFT[thispair][jjj]=Z_thispair
					FLUXRATIO[thispair][jjj]=FLUXRATIO_thispair
					EXP_SEP[thispair][jjj]=expected_separation
					MEAS_SEP[thispair][jjj]=meas_dist
					LBLUE[thispair][jjj]=thiswav
					YVAL[thispair][jjj]=thisy
					YVAL_STD[thispair][jjj]=stdy


	except:
		continue

COLS=[]
#NOW SAVE RESULTS
for thispair in range(len(PAIRS)):
	PAIRNAME=PAIRS[thispair]
	Ng=N[thispair]
	SNg=SNlist[thispair]
	Zg=REDSHIFT[thispair]
	Fg=FLUXRATIO[thispair]
	Es=EXP_SEP[thispair]
	Ms=MEAS_SEP[thispair]
	Lg=LBLUE[thispair]
	Yg=YVAL[thispair]
	stdYg=YVAL_STD[thispair]



	COLS.append(fits.Column(name='N_%s'%PAIRNAME,format='D',array=Ng))
	COLS.append(fits.Column(name='SN_%s'%PAIRNAME,format='D',array=SNg))
	COLS.append(fits.Column(name='z_%s'%PAIRNAME,format='D',array=Zg))
	COLS.append(fits.Column(name='lamb_%s'%PAIRNAME,format='D',array=Lg))
	COLS.append(fits.Column(name='x_%s'%PAIRNAME,format='D',array=(Lg-3E4)/9.75))

	COLS.append(fits.Column(name='fluxratio_%s'%PAIRNAME,format='D',array=Fg))
	COLS.append(fits.Column(name='meas_sep_%s'%PAIRNAME,format='D',array=Ms))
	COLS.append(fits.Column(name='exp_sep_%s'%PAIRNAME,format='D',array=Es))
	COLS.append(fits.Column(name='Yval_%s'%PAIRNAME,format='D',array=Yg))
	COLS.append(fits.Column(name='YvalSTD_%s'%PAIRNAME,format='D',array=stdYg))

new_cols = fits.ColDefs(COLS)
hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)

hdu.writeto(NEW_CATALOG, overwrite=True)
