# Imaging: 
## Detector1
  * exe_pipe_det1.py --> pipeline_Detector1.py 
  * Input: simulated_data/jwuncal.fits
  * Output: calibrated_det1/rate.fits
# Image2
  * qsub_exe_pipe_img2.sh > exe_pipe_img2.sh > pipeline_Image2.py
  * Input: calibrated_det1/rate.fits
  * Output: calibrated_img2/cal.fits
