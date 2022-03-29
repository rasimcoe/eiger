# Directory tree
.
├── simulated_data
│   ├── jw01243001001_01101_00001_nrca1_linear.fits
│   ├── jw01243001001_01101_00001_nrca1_uncal_cosmicrays.list
│   ├── jw01243001001_01101_00001_nrca1_uncal_F115W_CLEAR_final_seed_image.fits
│   ├── jw01243001001_01101_00001_nrca1_uncal_F115W_CLEAR_galaxy_seed_image.fits
│   ├── jw01243001001_01101_00001_nrca1_uncal_F115W_CLEAR_ptsrc_seed_image.fits
│   ├── jw01243001001_01101_00001_nrca1_uncal.fits
│   ├── jw01243001001_01101_00001_nrca1_uncal_galaxySources.list
│   ├── jw01243001001_01101_00001_nrca1_uncal_linear_dark_prep_object.fits
│   ├── jw01243001001_01101_00001_nrca1_uncal_pointsources.list
│
│
├── calibrated_det1
│   ├── jw01243001001_01101_00001_nrca1_ramp.fits
│   ├── jw01243001001_01101_00001_nrca1_rate.fits
│   ├── jw01243001001_01101_00001_nrca1_rateints.fits
│   ├── jw01243001001_01101_00001_nrca1_trapsfilled.fits
│
│
├── calibrated_img2
│   ├── jw01243001001_01101_00001_nrca1_cal.fits
│   ├── jw01243001001_01101_00001_nrca1_i2d.fits
│

# Imaging: 
## Detector1
  * exe_pipe_det1.py --> pipeline_Detector1.py 
  * Input: simulated_data/jwuncal.fits
  * Output: calibrated_det1/rate.fits
## Image2
  * qsub_exe_pipe_img2.sh --> exe_pipe_img2.sh --> pipeline_Image2.py
  * Input: calibrated_det1/rate.fits
  * Output: calibrated_img2/cal.fits
## Create cal.fits lists for globalsky
  * create_list_for_globalsky.py
