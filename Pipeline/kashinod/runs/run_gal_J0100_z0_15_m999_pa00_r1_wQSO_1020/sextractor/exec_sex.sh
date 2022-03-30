#!/bin/bash

outdir='results_calibrated_img3_wfss_medSubt41/'
mkdir $outdir

file_array=(\
'../calibrated_img3_wfss_medSubt41/jw01243_nrc_img3_wfss_visit1a_i2d.fits[1]' \
'../calibrated_img3_wfss_medSubt41/jw01243_nrc_img3_wfss_visit1b_i2d.fits[1]' \
'../calibrated_img3_wfss_medSubt41/jw01243_nrc_img3_wfss_visit2a_i2d.fits[1]' \
'../calibrated_img3_wfss_medSubt41/jw01243_nrc_img3_wfss_visit2b_i2d.fits[1]' \
'../calibrated_img3_wfss_medSubt41/jw01243_nrc_img3_wfss_visit3a_i2d.fits[1]' \
'../calibrated_img3_wfss_medSubt41/jw01243_nrc_img3_wfss_visit3b_i2d.fits[1]' \
'../calibrated_img3_wfss_medSubt41/jw01243_nrc_img3_wfss_visit4a_i2d.fits[1]' \
'../calibrated_img3_wfss_medSubt41/jw01243_nrc_img3_wfss_visit4b_i2d.fits[1]' \
#'../calibrated_img3_wfss_medSubt2/jw01243_nrc_img3_wfss_visit1a_i2d.fits[1]' \
#'../calibrated_img3_wfss_medSubt2/jw01243_nrc_img3_wfss_visit1b_i2d.fits[1]' \
#'../calibrated_img3_wfss_medSubt2/jw01243_nrc_img3_wfss_visit2a_i2d.fits[1]' \
#'../calibrated_img3_wfss_medSubt2/jw01243_nrc_img3_wfss_visit2b_i2d.fits[1]' \
#'../calibrated_img3_wfss_medSubt2/jw01243_nrc_img3_wfss_visit3a_i2d.fits[1]' \
#'../calibrated_img3_wfss_medSubt2/jw01243_nrc_img3_wfss_visit3b_i2d.fits[1]' \
#'../calibrated_img3_wfss_medSubt2/jw01243_nrc_img3_wfss_visit4a_i2d.fits[1]' \
#'../calibrated_img3_wfss_medSubt2/jw01243_nrc_img3_wfss_visit4b_i2d.fits[1]' \
#'../calibrated_img3/jw01243_nrc_img3_f200w_visit1_i2d.fits[1]' \
#'../calibrated_img3/jw01243_nrc_img3_f200w_visit2_i2d.fits[1]' \
#'../calibrated_img3/jw01243_nrc_img3_f200w_visit3_i2d.fits[1]' \
#'../calibrated_img3/jw01243_nrc_img3_f200w_visit4_i2d.fits[1]' \
#'../calibrated_img3/jw01243_nrc_img3_f200w_visit1_i2d.fits[1]' \
#'../calibrated_img3/jw01243_nrc_img3_f200w_visit2_i2d.fits[1]' \
#'../calibrated_img3/jw01243_nrc_img3_f200w_visit3_i2d.fits[1]' \
#'../calibrated_img3/jw01243_nrc_img3_f200w_visit4_i2d.fits[1]' \
#'../calibrated_img3/jw01243_nrc_img3_f356w_i2d.fits[1]' \
)


outname_array=(\
'sex_jw01243_nrc_img3_wfss_visit1a_i2d.cat' \
'sex_jw01243_nrc_img3_wfss_visit1b_i2d.cat' \
'sex_jw01243_nrc_img3_wfss_visit2a_i2d.cat' \
'sex_jw01243_nrc_img3_wfss_visit2b_i2d.cat' \
'sex_jw01243_nrc_img3_wfss_visit3a_i2d.cat' \
'sex_jw01243_nrc_img3_wfss_visit3b_i2d.cat' \
'sex_jw01243_nrc_img3_wfss_visit4a_i2d.cat' \
'sex_jw01243_nrc_img3_wfss_visit4b_i2d.cat' \
#'sex_jw01243_nrc_img3_wfss_visit1a_i2d.cat' \
#'sex_jw01243_nrc_img3_wfss_visit1b_i2d.cat' \
#'sex_jw01243_nrc_img3_wfss_visit2a_i2d.cat' \
#'sex_jw01243_nrc_img3_wfss_visit2b_i2d.cat' \
#'sex_jw01243_nrc_img3_wfss_visit3a_i2d.cat' \
#'sex_jw01243_nrc_img3_wfss_visit3b_i2d.cat' \
#'sex_jw01243_nrc_img3_wfss_visit4a_i2d.cat' \
#'sex_jw01243_nrc_img3_wfss_visit4b_i2d.cat' \
#'sex_jw01243_nrc_img3_f115w_visit1_i2d.cat' \
#'sex_jw01243_nrc_img3_f115w_visit2_i2d.cat' \
#'sex_jw01243_nrc_img3_f115w_visit3_i2d.cat' \
#'sex_jw01243_nrc_img3_f115w_visit4_i2d.cat' \
#'sex_jw01243_nrc_img3_f200w_visit1_i2d.cat' \
#'sex_jw01243_nrc_img3_f200w_visit2_i2d.cat' \
#'sex_jw01243_nrc_img3_f200w_visit3_i2d.cat' \
#'sex_jw01243_nrc_img3_f200w_visit4_i2d.cat' \
#'sex_jw01243_nrc_img3_f356w_i2d.cat' \
)

config_array=(\
'detect_wfss.sex' \
'detect_wfss.sex' \
'detect_wfss.sex' \
'detect_wfss.sex' \
'detect_wfss.sex' \
'detect_wfss.sex' \
'detect_wfss.sex' \
'detect_wfss.sex' \
#'detect_f115w.sex' \
#'detect_f115w.sex' \
#'detect_f115w.sex' \
#'detect_f115w.sex' \
#'detect_f200w.sex' \
#'detect_f200w.sex' \
#'detect_f200w.sex' \
#'detect_f200w.sex' \
#'detect_f356w.sex' \
)


# All
#sequence=`seq 0 16` 

# LW WFSS
#sequence=`seq 0 7` 
#sequence=`seq 0 1`
#sequence=`seq 2 7`

# SW F115W
#sequence=`seq 8 11`

# SW F200W
#sequence=`seq 12 15`

# LW F356W
#sequence=16


sequence=`seq 0 7`


for i in $sequence; do

    fil=${file_array[i]}
    wfil=${fil/"[1]"/"[2]"}
    out=$outdir${outname_array[i]}
    cfg=${config_array[i]}
    echo '---------------------------------'
    echo 'input : '$fil
    echo 'weight: '$wfil
    echo 'config: '$cfg
    echo 'output: '$out
    sed -e 's%_weight_image_%'$wfil'%g' \
	$cfg > $cfg'_tmp'
    ~/local/bin/sex -c $cfg'_tmp' $fil
    mv tmp.cat $out
done





