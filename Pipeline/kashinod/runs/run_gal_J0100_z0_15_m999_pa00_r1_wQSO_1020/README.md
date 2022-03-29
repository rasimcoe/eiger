# Imaging: 
* Detector1
  ** exe_pipe_det1.py --> pipeline_Detector1.py 
  ** Input: simulated_data/jwuncal.fits
  ** Output: calibrated_det1/rate.fits
## Image2
  * qsub_exe_pipe_img2.sh --> exe_pipe_img2.sh --> pipeline_Image2.py
  * Input: calibrated_det1/rate.fits
  * Output: calibrated_img2/cal.fits
## Image3 (standard; no manual processing)
  * exe_pipe_img3.sh --> pipeline_Image3.py
  * Input: calibrated_img2/cal.fits; asn.json
  * Output: calibrated_img3/cat.ecsv, i2d.fits 
## Create cal.fits lists for globalsky
  * create_list_for_globalsky.py
## Get global-sky images
  * exe_get_globalsky.sh --> get_globalsky.py
  * Input: list created above, calibrated_img2/cal.fits
  * Output: globalsky_nrca1_f115w.fits etc.
## Subtract global-sky & global-xy-medians from Image2/cal.fits
* qsub_exe_median_filter_img2cal.sh --> 

