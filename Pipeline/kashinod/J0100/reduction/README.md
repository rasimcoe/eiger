Last update: Daichi Kashino, 2022-09-14

These codes assume that uncal.fits files are stored in the directories:
`../data/uncal_F115W`
`../data/uncal_F200W`
`../data/uncal_F356W`
`../data/uncal_F356W_GRISM`

The bash scripts `exe_*.sh` and `qsub_exe_*.sh` are prepared to run the process for each of the all selected files using multi cores.  Usually it's enough to edit `qsub_exe_*.sh` where you can select files and specify output directory name etc.


# WFSS
- Detector1
  - `bash qsub_exe_pipe_det1.sh` to run `pipeline_Detector1.py`
  - Input files : e.g., `../data/uncal_F356W_GRISM/jw01243001001_02101_00001_nrcalong_uncal.fits`
  - Output files: e.g., `calibrated_det1/jw01243001001_02101_00001_nrcalong_rate.fits`

### Image 2 for WFSS

- Modify WFSS fits header
  - `python modify_wfss_fits_header.py`
  - Input dir : `calibrated_det1`
  - Output dir: `calibrated_det1_wfss_hdr_corr`

- Execute Image2
  - `bash qsub_exe_pipe_img2_wfss.sh` to run `pipeline_Image2.py`
  - Output files: `calibrated_img2_wfss/jw01243001001_02101_00001_nrcalong_cal.fits`

- "Global" (column) median subtraction -- subtract median values at each column of the images
  - `bash qsub_exe_subtract_globalMed_wfss.sh` to run `subtract_globalMed_wfss.py`
  - Input: `calibrated_img2_wfss/jw01243001001_02101_00001_nrcalong_cal.fits`
  - Output: `calibrated_img2_wfss_globalMedSubt/`
    - `jw01243001001_02101_00001_nrcalong_cal_globalMed.fits` -- median image
    - `jw01243001001_02101_00001_nrcalong_cal_globalMedSubt.fits` -- median-subtracted image
    
- Continuum subtraction to get EMLINE and CONTINUA images
  - `bash qsub_exe_subtract_continua_wfss.sh` to run `subtract_continua_wfss.py`
  - Input dir: `calibrated_img2_wfss_globalMedSubt`
  - Output dir: `calibrated_img2_wfss_globalMedSubt_contSubt_kx51_9_Skx21_5`
    - `jw01243001001_02101_00001_nrcalong_cal_globalMedSubt_continua.fits` - CONTINUA image
    - `jw01243001001_02101_00001_nrcalong_cal_globalMedSubt_continua_lk.fits`
    - `jw01243001001_02101_00001_nrcalong_cal_globalMedSubt_continua_sk.fits`
    - `jw01243001001_02101_00001_nrcalong_cal_globalMedSubt_continua_wht.fits`
    - `jw01243001001_02101_00001_nrcalong_cal_globalMedSubt_emline.fits` - EMLINE image

- Create lists for globalsky (or master bias) image for each detector (Mod A and B)
  - `list_nrcalong_wfss_cal_globalMedSubt_contSubt_kx51_9_Skx21_5.txt`
  - `list_nrcblong_wfss_cal_globalMedSubt_contSubt_kx51_9_Skx21_5.txt`

- Get globalsky (master bias) from EMLINE images
    - `python get_globalsky.py list_fil outname --mask_brightpixels=True`
      - list_fil = `list_nrca[b]long_wfss_cal_globalMedSubt_contSubt_kx51_9_Skx21_5.txt`
      - outname = `nrca[b]long_wfss_cal_globalMedSubt_contSubt_kx51_9_Skx21_5`
      - Output:
      	- `globalsky_nrcalong_wfss_cal_globalMedSubt_contSubt_kx51_9_Skx21_5.fits`
      	- `globalsky_nrcblong_wfss_cal_globalMedSubt_contSubt_kx51_9_Skx21_5.fits`


### Go back to the "raw" cal.fits with the globalsky (master-bias) in hand

- Subtract globalsky from raw cal.fits files
    - `bash qsub_exe_subtract_globalsky_wfss.sh` to run `subtract_globalsky.py`
    - globalsky_img: `globalsky_nrca[b]long_wfss_cal_globalMedSubt_contSubt_kx51_9_Skx21_5.fits`    
    - Input: e.g., `calibrated_img2_wfss/jw01243001001_02101_00001_nrcalong_cal.fits`
    - Output: e.g., `calibrated_img2_wfss_globalskySubtV3/jw01243001001_02101_00001_nrcalong_cal.fits`  ## V3 has no meaning.

- Again, "global" median subtraction for master-bias-subtracted images -- subtract median values at each column of the images
  - `bash qsub_exe_subtract_globalMed_wfss.sh` to run `subtract_globalMed_wfss.py`
  - Input: e.g., `calibrated_img2_wfss_globalskySubtV3/jw01243001001_02101_00001_nrcalong_cal.fits`
  - Output: e.g., `calibrated_img2_wfss_globalskySubtV3_globalMedSubt/`
    - `jw01243001001_02101_00001_nrcalong_cal_globalMed.fits` -- median image
    - `jw01243001001_02101_00001_nrcalong_cal_globalMedSubt.fits` -- median-subtracted image

- Continuum subtraction
  - `bash qsub_exe_subtract_continua_wfss.sh` to run `subtract_continua_wfss.py`
  - Input dir: `calibrated_img2_wfss_globalskySubtV3_globalMedSubt/_cal_globalMedSubt.fits`
  - Output dir: `calibrated_img2_wfss_globalskySubtV3_globalMedSubt_contSubt_kx51_9_Skx21_5`
    - `jw01243001001_02101_00001_nrcalong_cal_globalMedSubt_continua.fits` - continuum image
    - `jw01243001001_02101_00001_nrcalong_cal_globalMedSubt_continua_lk.fits`
    - `jw01243001001_02101_00001_nrcalong_cal_globalMedSubt_continua_sk.fits`
    - `jw01243001001_02101_00001_nrcalong_cal_globalMedSubt_continua_wht.fits`
    - `jw01243001001_02101_00001_nrcalong_cal_globalMedSubt_emline.fits` - emission-line image

- Create emission-line masks
  - at `./GrismMasking/`
  - run `execute_sextractor_for_emline_mask.ipynb`
  - and then `apply_emline_mask.ipynb`
    - Output mask files: e.g., `calibrated_img2_wfss/jw01243001001_02101_00001_nrcalong_mask.fits`

### Do everything again but now with mask!

- Global median subtraction with mask
  - `bash qsub_exe_subtract_globalMed_wfss_with_mask.sh` to run `subtract_globalMed_wfss.py` with `--mask` option
    - Input dir: `calibrated_img2_wfss`
    - Output dir: `calibrated_img2_wfss_emlineMasked_globalMedSubt`

- Continuum subtraction with mask
  - `bash qsub_exe_subtract_continua_wfss_with_mask.sh` to run `subtract_continua_wfss.py` with `--mask` option
    - Input dir: `calibrated_img2_wfss_emlineMasked_globalMedSubt`
    - Output dir: `calibrated_img2_wfss_emlineMasked_globalMedSubt_contSubt_kx51_9_Skx21_5`

- Create lists for global-sky (master bias) image for each detector
  - `list_nrcalong_wfss_cal_emlineMasked_globalMedSubt_contSubt_kx51_9_Skx21_5.txt`
  - `list_nrcblong_wfss_cal_emlineMasked_globalMedSubt_contSubt_kx51_9_Skx21_5.txt`

- Get globalsky with mask
  - `python get_globalSky.py  list_fil outname --mask_brightpixels=True --mask_list=mask_list
    - list_name = `list_nrca[b]long_wfss_cal_emlineMasked_globalMedSubt_contSubt_kx51_9_Skx21_5.txt`
    - out_name = `nrca[b]long_wfss_cal_emlineMasked_globalMedSubt_contSubt_kx51_9_Skx21_5`
    - mask_list = ``
    - Output:
      - `globalsky_nrcalong_wfss_cal_emlineMasked_globalMedSubt_contSubt_kx51_9_Skx21_5.fits`
      - `globalsky_nrcblong_wfss_cal_emlineMasked_globalMedSubt_contSubt_kx51_9_Skx21_5.fits`
    
- Subtract global sky from wfss cal.fits files
    - `bash qsub_exe_subtract_globalsky_wfss.sh` to run `subtract_globalsky.py`
    - globalsky_img: `globalsky_nrca[b]long_wfss_cal_emlineMasked_globalMedSubt_contSubt_kx51_9_Skx21_5.fits`
    - Input: `calibrated_img2_wfss/cal.fits`
    - Output: `calibrated_img2_wfss_globalskySubtV4/cal.fits`  # V4 has no meaning.

- Again "global" median subtraction for globalsky-subtracted images
  - `bash qsub_exe_subtract_globalMed_wfss.sh` to run `subtract_globalMed_wfss.py`
  - Input: e.g., `calibrated_img2_wfss_globalskySubtV4/cal.fits`
  - Output: e.g., `calibrated_img2_wfss_globalskySubtV4_globalMedSubt/cal_.fits` 

- Continuum subtraction
  - `bash qsub_exe_subtract_continua_wfss.sh` to run `subtract_continua_wfss.py`
  - Input dir: `calibrated_img2_wfss_globalskySubtV4_globalMedSubt/_cal_globalMedSubt.fits`
  - Output dir: `calibrated_img2_wfss_globalskySubtV4_globalMedSubt_contSubt_kx51_9_Skx21_5`
    - `jw01243001001_02101_00001_nrcalong_cal_globalMedSubt_continua.fits` - continuum image
    - `jw01243001001_02101_00001_nrcalong_cal_globalMedSubt_continua_lk.fits`
    - `jw01243001001_02101_00001_nrcalong_cal_globalMedSubt_continua_sk.fits`
    - `jw01243001001_02101_00001_nrcalong_cal_globalMedSubt_continua_wht.fits`
    - `jw01243001001_02101_00001_nrcalong_cal_globalMedSubt_emline.fits` - emission-line image

### Okay, now Image3 to stack exposures  (per visit and module)

- Create asn.json file for Image3
  - per visit and module
    - e.g., `img3_F356W_wfss_globalskySubtV4_emlineMasked_globalMedSubt_contSubt_kx51_9_Skx21_5_emline_visit4b_pdit123_asn.json`
    - e.g., `img3_F356W_wfss_globalskySubtV4_emlineMasked_globalMedSubt_contSubt_kx51_9_Skx21_5_continua_visit4b_pdit123_asn.json`

- Image3
  - `python pipeline_Image3noSkyMatch_wfss.py asn_file out_dir`
  - asn_file: e.g., `img3_F356W_wfss_globalskySubtV4_emlineMasked_globalMedSubt_contSubt_kx51_9_Skx21_5_emline_visit4b_pdit123_asn.json`
  - out_dir: `calibrated_img3noskymatch_wfss_globalskySubtV4_emlineMasked_globalMedSubt_contSubt_kx51_9_Skx21_5`






-------------------------------------------------------------------------------
# OLD INFORMATION BELOW (for simulated data)


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


