#EIGER NIRCam/HST imaging reduction  \
These scripts download the raw files, build the directory structure, run the reduction and produce mosaics.

Requirements:  \
Needs a lot of disk space (800 GB+)  \
Current version of the jwst pipeline and all it's dependencies  \
A MAST astroquery token (https://auth.mast.stsci.edu/tokens)  \
source-extractor

The basic scripts are:  \
step0__MAST_query_organise.py \
Fetches the data from MAST using astroquery, it then checks the files are complete, builds the directory structure and organises the uncal.fits files.
Check the files are all there, astroquery sometimes misses some. Currently the list of files is quite manual, you need to set the sequence ids of the different parts of the visits.
You can rerun the script and it should find them. \
There are three command line options for this script, but only this one. \
    -q or --qso quasar name on MAST  e.g. "2MASS J01001301+2802257" "QSO J1120+0641" \
    -d or --di base directory e.g. "/scratch/mruari/EIGER/imaging/J0100+2802/"" \
    -p or --pmap CRDS pmap e.g. "jwst_0988.pmap"        #use the latest one \
These options are writen to a params yaml file which the other steps read.

step1__imaging_pipe1_basic_findgrid.py \
Does the basic reduction on each individual frame (Image1 and Image2),
then it jointly finds a wcs grid which will fit all files (which is saved and read in later steps).
It then produces the basic mosaic with Image3.
Note that if you skip reject_outliers in Image3 at this point it will break the next step. Without this the relative wcs tweaks are not saved.
Also for smaller machines reduce n_procs, this is the number of parallel jobs that will be run. Set to 40.
It is slow to start because the first group of exposures are run one by one, this is to let the crds downlad the reference files to the cache before running in parallel.

step2__imaging_pipe2_wcs.py \
Does the absolute wcs calibration. The F356W image is aligned to Gaia and restacked.
The F200W and F115W mosaics are then aligned to the new F356W. These adjustments are propagated into every pipe2 file, so that all the individual images are aligned.

step3__imaging_pipe3_wisp_skyfix.py \
This step creates median sky flats to remove wisp features for the short wavelength channels. It makes a source mask on the fly for each image, this is quite slow.
It writes these files to the mycals directory in each filter directory, these are masked median stacks.
It then applies these to the individual image crf.fits files, along with sky subtraction, snowball masking. It then makes new stacks.

step4__imaging_cal2_sourcemask_wisp.py \
Uses the new stacks to make a deep sourcemask. This sourcemask is then used to make second generation wisp templates.

step5__imaging_pipe4_filter_stack.py \
Fixes the crf.fits files with the better wisp templates, as well as snowball masking. It also uses the deep source mask to filter the 1/f noise and conduct sky subtraction.
Finally it stacks the data into final mosaics.

Notes: \
It is not fully sequential. You cannot delete all the files in the last step after you run the next. step5 reads the crf.fits files from step2/pipe2 for example. \
You can delete useless files like the \*outlier_i2d.fits made in outlier rejection, or the blank \*trapsfilled.fits and intermediate \*rate.fits \*rateints.fits. \
If you are running out of space delete: \
 the redundant mast download dir (download/mastDownload), which are copied into download/organised_output. \
 the FXXXW/pipe1_basic/jw\*.fits files after running step2. \
 the FXXXW/pipe3_skyfix/jw\*.fits files after running step4.



Finally reducing HST data: \
Requirements: \
stenv, which includes drizzlepac. \
Also source-extractor.

HST_pipe__dev2.1_gaia_nrctweak.py \
Is convoluted and unstable. It needs the F356W step2 stack to exist, to align the images to.
Basically it searches MAST for the object names of exposures taken near the quasar, and the filters.
It then loops through the filters, downloads the flc.fits from MAST. 
It then breaks them into epochs, and aligned them to Gaia (due to proper motion errors). It then stacks these with cosmic ray rejection.
In the end it stacks again, using the cosmic ray masks from before.
Then lastly a catalog is run to align it to the NIRCam reference image carefully. This is propagated to all the files and it runs one last stack.
It also saves the zeropoints to the file headers.

Warning: \
Sometimes there will be a blind pointing near the quasar with the target name "ANY" or "DARK". 
The script will then try to download every MAST file with that filter and target name. This is a lot of files. 
Manually screen the names, and remove any filters or targets you aren't interested in. It was set up this way due to a limitation in astroquery.

