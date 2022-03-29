# Imaging: 
## Detector1
 - exe_pipe_det1.py --> pipeline_Detector1.py 
 - Input: simulated_data/uncal.fits
 - Output: calibrated_det1/rate.fits
## Image2
 - qsub_exe_pipe_img2.sh --> exe_pipe_img2.sh --> pipeline_Image2.py
 - Input: calibrated_det1/rate.fits
 - Output: calibrated_img2/cal.fits
## Image3 (standard; no manual processing)
 - qsub_exe_pipe_img3.sh --> exe_pipe_img3.sh --> pipeline_Image3.py
 - Input: calibrated_img2/cal.fits; asn.json
 - Output: calibrated_img3/cat.ecsv, i2d.fits 
## Create cal.fits lists for globalsky
 - create_list_for_globalsky.py
## Get global-sky images
 - exe_get_globalsky.sh --> get_globalsky.py
 - Input: list created above, calibrated_img2/cal.fits
 - Output: globalsky_nrca1_f115w.fits etc.
## Subtract global-sky & global-xy-medians from Image2/cal.fits
 - qsub_exe_median_filter_img2cal.sh --> exe_median_filter_img2cal.sh --> median_filter_img2cal_fits_v20211029.py
 - Input: calibrated_img2/cal.fits
 - Output: calibrated_img2_medSubt/cal.fits
## Create asn files for Image3
 - create_asn_img3.ipynb
 - Output: list_nrcb5_f115w.txt, list_nrcb5_f356w_imaging.txt etc.
## Image 3 from the processed images
 - qsub_exe_pipe_img3.sh --> exe_pipe_img3.sh --> pipeline_Image3.sh
 - Input: asn file e.g., jw01243_nrc_img3_f115w_visit1_medSubt_asn.json
 - Output: calibrated_img3_medSubt/i2d.fits, cat.ecsv
## Masking
 - Prepare a ds9-format region file from Image3/cat.ecsv, your own Sextractor result, or input catalog.
 - qsub_exe_mask_objects_img2cal.sh --> exe_mask_objects_img2cal.sh --> exe_mask_objects_img2cal.sh
 - Input: calibrated_img2/cal.fits
 - Output: calibrated_img2_mask/mask.fits
## Re-create globalsky with mask
* exe_get_globalsky.sh --> get_globalsky.py
* Input: list created above, mask file list, calibrated_img2/cal.fits
* Output: globalsky_nrca1_f115w_masked.fits etc
## 

# WFSS
## Detector1
 * Same as Imaging
## Modify WFSS fits header
 * modify_wfss_fits_header.py
 * Input dir: calibrated_det1
 * Output dir: calibrated_det1_wfss_hdr_corr
## WFSS Image2:
 * qsub_exe_pipe_img2_wfss.sh --> exe_pipe_img2_wfss.sh --> pipeline_Image2_wfss.sh
 * Output dir: calibrated_img2_wfss
## Subtract global medians (at each x) and sources' continua
 * qsub_exe_median_filter_wfss_img2cal.sh --> exe_median_filter_wfss_img2cal.sh --> median_filter_wfss_img2cal_fits_v*.py
 * Kernel: kx=51, ky=1, kx_gap=9
 * Input: calibrated_img2_wfss/cal.fits
 * Output: calibrated_img2_wfss_medFiltered/cal.fits
## Create lists for global sky images
 * create_lists_for_globalsky.ipynb
 * Output: list_nrca5_f356w_wfss.txt etc 
## Globalsky from median-filtered emission-line images
  *
