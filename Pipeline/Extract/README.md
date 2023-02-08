This directory contains script to:

1) Extract and combine 2D spectra for direct-image detected sources (or at any RA,DEC)
extract_J0100.py -- J0100+2802 field
eiger_tracing_dk.py -- contains various relevant functions



Requirements:
1. grismconf V4 https://github.com/npirzkal/GRISM_NIRCAM
2. eiger reference files, trace cor  "'yoffset_polyreg_F356W.R*.npy" and wavelength cor "lambda_offset*.npy"
3. WCS cor files for each visit & module 
	- direct images   "img3_F356W_visit[1,2,3,4][a,b]_pdit123_i2d.fits"     ***Daichi, could you add how those were created?
	- the function "eiger_tracing_dk.radec_in_this_vismod" has a library with WCS offsets per field
4. EIGER sensitivity curves "NIRCam.F356W.R.A.1st.sensitivity.EIGER.fits"


