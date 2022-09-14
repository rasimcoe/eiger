Last update: Daichi Kashino, 2022-09-14

These codes assume that uncal.fits files are stored in the directories:
`../data/uncal_F115W`
`../data/uncal_F200W`
`../data/uncal_F356W`
`../data/uncal_F356W_GRISM`

The bash scripts `exe_*.sh` and `qsub_exe_*.sh` are prepared to run the process for each of the all selected files using multi cores.  Usually it's enough to edit `qsub_exe_*.sh` where you can select files and specify output directory name etc.


## WFSS:
- Detector1
  - `qsub_exe_pipe_det1.sh` >> `pipeline_Detector1.py`
  - Input: `uncal.fits`
  - Output: `calibrated_det1/rate.fits`

### Imaging-mode reduction
- Modify WFSS fits header
    - `modify_wfss_fits_header.py`
    - Input dir: `calibrated_det1`
    - Output dir: `calibrated_det1_wfss_hdr_corr`

- Image2:
    - `qsub_exe_pipe_img2_wfss.sh` >> `pipeline_Image2.py`
    - Output: `calibrated_img2_wfss/cal.fits`

- "global" median-filtering - subtract median values at each column of the images 
    - `qsub_exe_subtract_globalMed_wfss.sh` >> `subtract_globalMed_wfss.py`
    - Input: e.g., `calibrated_img2_wfss/cal.fits`
    - Output: e.g., `calibrated_img2_wfss_globalMedSubt/cal_.fits`

- Create lists for global-sky (or master bias) image for each detector
    - Output: e.g., `list_nrcalong_wfss.txt`
      - calibrated_img2_wfss_medSubt/jw01243001002_02101_00001_nrcalong_cal.fits
      - calibrated_img2_wfss_medSubt/jw01243001002_02101_00002_nrcalong_cal.fits

- Get global sky from median-filtered (iter1) images
    - python ../scripts/get_globalsky.py list_nrca5_f356w_wfss.txt nrca5_f356w_wfss
    - python ../scripts/get_globalsky.py list_nrcb5_f356w_wfss.txt nrcb5_f356w_wfss
    - Output: globalsky_nrcb5_f356w_wfss.fits

- Subtract global sky from median-filtered (iter1) images
    - scripts: ./qsub_exe_subtract_globalsky.sh > exe_subtract_globalsky.sh > python ../scripts/subtract_globalsky.py $fil $output_dir $global_sky
    - Input: calibrated_img2_wfss.fits
    - Output: calibrated_img2_wfss_medSubt2/cal.fits

- WFSS Image3
    - create_asn_img3_wfss.ipynb
        - asn_files
    - scripts: exe_pipe_img3_wfss.sh > pipeline_Image3_wfss.py
    - Input: calibrated_img2_wfss_medSubt2/cal.fits
    - Output dir: calibrated_img3_wfss_medSubt2

### imaging-mode reduction with masking
- Prepare masks
    - Sextractor
    - ./sextractor/results/sex_jw01243_nrc_img3_wfss_visit3a_i2d_selected.reg
- Masking 
    - scripts: qsub_exe_mask_objects_wfss_img2cal.sh
    - Input: calibrated_img2_wfss/cal.fits
    - Output: calibrated_img2_wfss_mask/cal.fits
- Median-filtering (iter2) with mask
    - qsub_exe_median_filter_wfss_img2cal_with_mask.sh > exe_median_filter_wfss_img2cal_with_mask.sh > 
python ../scripts/median_filter_wfss_img2cal_with_mask.py $fil $output_dir $kx $ky $kx_gap $mask_fil
        - kernel (x, y, x_gap): 51 1 9
        - ==> maskをしているので、x_gap = 1にする（自分自身で引くのを避けるため、必ずx_gap>=1とする）
    - Input: calibrated_img2_wfss_mask/cal.fits
    - Output: calibrated_img2_wfss_medSubt3/cal.fits
- WFSS Image3
    - Input: calibrated_img2_wfss_medSubt3/cal.fits
    - Output: calibrated_img3_wfss_medSubt3/cal.fits
- Get global sky from median-filtered (iter2, with mask) images
    - Create lists
        - list_nrcb5_f356w_wfss_iter2.txt
    - python ../scripts/get_globalsky.py list_nrca5_f356w_wfss_iter2.txt
    - Output: globalsky_nrca5_f356w_wfss_iter2.fits
- Subtract global sky from median-filtered (iter2) images
    - scripts: qsub_exe_subtract_globalsky.sh
    - Input: calibrated_img2_wfss_medSubt3/cal.fits
    - Output: calibrated_img2_wfss_medSubt4/cal.fits
- WFSS Image3 from fully-processed images
    - create_asn_img3_wfss.ipynb
    - asn files: jw01243_nrc_img3_wfss_visit1a_medSubt41_asn.json
    - Output: calibrated_img3_wfss_medSubt41



## Imaging
### Standard
- Detector1
  - scripts: exe_pipe_det1.py > pipeline_Detector1.py
  - Input: simulated_data/uncal.fits
  - Output: calibrated_det1/rate.fits

- Image2
  - scripts: qsub_exe_pipe_img2.sh > exe_pipe_img2.sh > pipeline_Image2.py
  - Input: calibrated_det1/rate.fits
  - Output: calibrated_img2/cal.fits

- Create asn file for Image3

- Image3 
  - scripts: exe_pipe_img3.sh > pipeline_Image3.py
  - Input: calibrated_img2/cal.fits and asn.json
  - Output: calibrated_img3/cat.ecsv, i2d.fits etc.

### with global sky + median subtraction
- Create lists for creating global sky images
  - scripts: create_lists_for_global_sky.ipynb
  - Output: e.g., list_nrca1_f200w.txt, list_nrcb5_f356w_imaging.txt

- Get global sky images
  - scripts: exe_get_global_sky_image.sh; python ../scripts/get_global_sky_image.py [list file] [outname]
  - Input: calibrated_img2/cal.fits
  - Output: globalsky_nrca1_f115w.fits etc.

- Subtract Global sky (and then) median values (in horizontal and vertical directions)
  - scripts: exe_median_filter_img2cal.sh; python ../scripts/median_filter_img2cal_fits_v20211029.py
  - Input: calibrated_img2/cal.fits
  - Output: calibrated_img2_medSubt/cal.fits

- Create asn file for Image3
 
- Image3 from the processed images
  - Input: 
  - Output: calibrated_img3_medSubt/cat.ecsv, i2d.fits etc.

### with masking detected sources
- Create masks (.reg files) from input images
  - scripts: create_mask_for_median_filtering.ipynb
  - Output: mask_f115w_m280.reg, mask_f200w_m280.reg, mask_f356w_m270.reg

- Apply the masks
  - scripts: ../scripts/mask_objects_img2cal.py
  - Input: calibrated_img2/cal.fits
  - Output: calibrated_img2_mask_f115w_m280/jw01243001003_01101_00063_nrcb2_mask.fits

- Re-creating global sky images
  - Copy and modify the mask lists; list_nrca1_f115w.txt --> mask_list_nrca1_f115w.txt
  - execute ./exe_get_global_sky_image.sh
  - Input: e.g., calibrated_img2_mask_f200w_m280/jw01243001003_01101_00077_nrca1_mask.fits
  - Output: e.g., globalsky_nrca1_f115w_masked.fits

- Global sky + median subtraction
  - execute python ../scripts/median_filter_img2cal_fits_v20211104.py $fil $output_dir
  - Input: 
  - Output: calibrated_img2_medSubt
 
 - Image3 from processed image (fully processed)
  - create_asn_img3.ipynb
  - Input: 
    - asn file: jw01243_nrc_img3_f356w_medSubt_v3_asn.json
  - Output: calibrated_img3_medSubt_v3


