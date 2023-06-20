#EIGER multiband photometry
These scripts create PSF matched images and perform source detection and photometry.

Requirements:
webbpsf
photutils
source-extractor
TinyTim (HST only)

You need to copy the final images to a mosaic directory, by default basedir/junk/current_best/

You also need to copy the NIRCam source masks to current_best/src_mask/

Before running HST data you need to calculate TinyTim model PSFs. Otherwise comment out HST section of photometry1.

The scripts do the following:
photometry1__PSF_convolve.py
This step finds the files in the mosaic directory, you need to specify which HST filters you are using.
It writes out a parameter file with the data of each filter mosaic.
It PSF matches the image, convolving the images with the matching filter. For JWST these are from webbpsf, for HST they need to be calculated manually with TinyTim.

photometry2__sextractor.py
Does the multi-band runs with source-extractor. It writes the images to it's working directory, and runs each band with the direction image (F356W).
It writes a combined catalog.

photometry3__noise_model
Builds a noise model based on random apertures in blank sections of the image. It uses the source mask and random locations. It fits the function as a power law vs aperture area.
These model parameters are saved for the next step.

photometry4__apply_noise_model.py
For each source in the sextractor-catalog the code calculates the aperture area, and scales the error model to the local image WHT.
These are the errors we use in the photometry. It writes a full and short version of the catalog.



TinyTim:
HST's psf model is archaic. Use the following steps with TinyTim:
1) Get the date and time from the image.
2) Find the quasar postion in the image (roughly). X and Y, CCD1 or 2 for ACS or WFC3UVIS.
3) Get the focus position (http://focustool.stsci.edu/cgi-bin/control.py)

4) run "tiny1 J0000_FXXXW_ACS.in"
       select instrument:  15 ACS - Wide Field Channel
       specify the CCD and filter
    SED: Power law : F(nu) = nu^i
    alpha: 0.0
    PSF diameter arcsec: 2.0
    rootname: FXXXW_tinytim
5) run the tiny2 command printed at the end

6) run "tiny3 J0000_FXXXW_ACS.in SUB=5"
7) copy the files to basename/HST/PSF/FXXXW/
    the imporant file is FXXXW_tinytim00.fits, which is the subsampled PSF (relative to HST pixels) and the charge diffusion kernel



