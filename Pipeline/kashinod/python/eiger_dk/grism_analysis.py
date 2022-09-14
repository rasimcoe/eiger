import os,copy
import numpy as np

from astropy.io.ascii import SExtractor
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs import utils
import astropy.units as u 
import numpy as np
from astropy.convolution import Tophat2DKernel,Gaussian2DKernel, convolve

import grismconf
from jwst import datamodels

#def get_conf(filt, grism, module):


def get_grism_sensitivity(filt, grism, module, directory='/scratch/kashinod/EIGER/eiger_reference_files'):
    sensfil = os.path.join(directory, 
                           'NIRCam.'+filt.upper()+'.'+grism.upper()+'.'+module.upper()+'.1st.sensitivity.EIGER.fits')
    sens = fits.getdata(sensfil,1)
    return sens

def plot_reticle(ax, x, y, 
                 xoff=5.0, yoff=5.0, 
                 s=50, alpha=1.0, color='white', linewidth=[1],**kwargs):
    xoff_list = xoff * np.array([-1,1,0,0])
    yoff_list = yoff * np.array([0,0,1,-1])
    #for m in [0,1,2,3]: 
    for m in [2,3]: 
        ax.scatter(x + xoff_list[m], 
                   y + yoff_list[m], 
                   marker=m, s=s, alpha=alpha, color=color, linewidth=linewidth, **kwargs)


def plot_reticle_c(ax, x, y, 
                 xoff=5.0, yoff=5.0, 
                 **kwargs):
    xoff_list = xoff * np.array([-1,1,0,0])
    yoff_list = yoff * np.array([0,0,1,-1])
    #for m in [0,1,2,3]: 
    for m in [2,3]: 
        ax.scatter(x + xoff_list[m], 
                   y + yoff_list[m], 
                   marker=m, **kwargs)

        
def read_sextractor_catalog(cat_fil, wcs):
    pixel_area = utils.proj_plane_pixel_area(wcs) * u.deg**2
    print('Pixel area:  ', pixel_area)
    print('Pixel area:  ', pixel_area.to(u.arcsec**2))
    print('Pixel area:  ', pixel_area.to(u.sr))
    pixel_area_sr = pixel_area.to_value(u.sr)
    #pixel_scales = utils.proj_plane_pixel_scales(wcs) * u.deg
    #print('Pixel scale: ', pixel_scales)
    #print('Pixel scale: ', pixel_scales.to(u.arcsec))
    
    MAG_ZEROPOINT = 8.90 - 2.5*(np.log10(pixel_area_sr) + 6.0) 
    print('MAG_ZEROPOINT: ', MAG_ZEROPOINT)
    
    cat = SExtractor().read(cat_fil)
    keys = cat.keys()
    
    FLUX_keys = [s for s in keys if 'FLUX_' in s]
    FLUXERR_keys = [s for s in keys if 'FLUXERR_' in s]
    MAG_keys = [s for s in keys if 'MAG_' in s]

    for k in FLUX_keys:
        cat[k]=cat[k] * pixel_area_sr * 1e15
        cat[k].unit = u.nJy
    for k in FLUXERR_keys:
        cat[k]=cat[k] * pixel_area_sr * 1e15
        cat[k].unit = u.nJy
    for k in MAG_keys:
        cat[k]=cat[k] + MAG_ZEROPOINT

    return cat



def exec_sextractor(sci_img,
                    err_img,
                    header,
                    name,
                    out_dir='./',
                    overwrite=False,
                    rectify_images=False,
                    sexcmd = 'source-extractor',
                    detect_sex='detect_wfss_for_fft_mask.sex'):
        
    if rectify_images:
        sci_img, err_img = rectify_images_for_sextractor(sci_img, err_img)
        
    fits.writeto('detection_image.fits', sci_img, header, overwrite=overwrite)
    fits.writeto('rms_weight.fits', err_img, header, overwrite=overwrite)
    
    cmd = sexcmd+' detection_image.fits -c '+detect_sex
    os.system(cmd)
                    
    # Rename
    segm_fil = os.path.join(out_dir, 'segm_'+name+'.fits')
    print('Save ', segm_fil)
    os.system('mv check_segm.fits '+segm_fil)
                    
    cat_fil = os.path.join(out_dir, 'sex_'+name+'.cat')
    print('Save ', cat_fil)
    os.system('mv tmp.cat '+cat_fil)
    return cat_fil, segm_fil, sci_img, err_img


def rectify_images_for_sextractor(sci_img, err_img):
    out_sci_img = copy.deepcopy(sci_img)
    out_err_img = copy.deepcopy(err_img)

    ### Modify data for sextractor
    ### Interpolate pixels with ERR=0
    ### This should be removerd when the pipeline 
    out_err_img[err_img==0] = np.nan
    out_err_img = convolve(out_err_img, np.ones((11,11)), 
                           preserve_nan=False)

    ### SCI ==> 0 where ERR == nan (out-of-field region) 
    out_sci_img[np.isnan(err_img)]= 0.0
    out_sci_img[np.isnan(sci_img)]= 0.0

    ### Increase the error level when ERR=NaN
    mederr = np.nanmedian(out_err_img)
    out_err_img[np.isnan(err_img)]=mederr*10
    out_err_img[np.isnan(sci_img)]=mederr*10
    
    return out_sci_img, out_err_img
    
def get_mask_for_emlines_for_FT(hdul, 
                                mask_traces, 
                                name,
                                out_dir='./',
                                overwrite=False,
                                sexcmd = 'source-extractor',
                                detect_sex='detect_wfss_for_fft_mask.sex'):
    
    sci_img = hdul['SCI'].data
    err_img = hdul['ERR'].data
    header = hdul['SCI'].header
    
    ### Rectify
    ### ERR=0 ==> ERR to be interpolated
    ### ERR=NaN ==> ERR should be large value, SCI should be 0
    ### SCI=NaN ==> ERR should be large value, SCI should be 0
    out_sci_img, out_err_img = rectify_images_for_sextractor(sci_img, err_img)
    
    ### Replace the pixels identified as traces
    out_sci_img[mask_traces==1.0] = 0.0
    mederr = np.nanmedian(out_err_img)
    out_err_img[mask_traces==1.0] = mederr*10
    
    cat_fil, segm_fil, img_use, err_used = exec_sextractor(
        out_sci_img,out_err_img,header,name,
        out_dir=out_dir,
        overwrite=overwrite,
        sexcmd=sexcmd, detect_sex=detect_sex)

    hdu_segm = fits.open(segm_fil)
    segm_data = hdu_segm['SCI'].data
    segm_data[segm_data>0]=1
    hdu_segm.close()
    
    mask_emlines = convolve(segm_data, Tophat2DKernel(5))
    mask_emlines[mask_emlines>0]=1.0
    return mask_emlines, out_sci_img, out_err_img
    


def get_mask_traces(cont_img, 
                    thresh=0.07,
                    kernel=np.ones((3,11))):

    #hdr = hdul_continua['SCI'].header
    #cont_img = hdul_continua['SCI'].data
    #err_img = hdul_continua['ERR'].data
    #wht_img = hdul_continua['WHT'].data
    
    ### Selection using CONTINUA image
    sel_traces = cont_img > thresh
    
    mask_traces = np.zeros(cont_img.shape)
    mask_traces[sel_traces]=1
    
    ## Little bit expand
    if kernel:
        print('Apply kernel to expand')
        mask_traces = convolve(mask_traces, kernel)
        mask_traces[mask_traces>0]=1

    return mask_traces


def FT_filter_grism_image(image,
                          wnlim_x=115,
                          wnlim_y=115,
                          remove_mean=True, 
                          FTfilt='Gauss'):

    ny0, nx0 = image.shape
    shifted_f_uv, shifted_ps2d, im = FT_grism_image(image, fftshift=True)
    
    xcen = shifted_f_uv.shape[1]/2-0.5
    ycen = shifted_f_uv.shape[0]/2-0.5
    yind, xind = np.indices(shifted_f_uv.shape)

    
    dy = (yind-ycen)
    dx = (xind-xcen)
    d = np.sqrt(dy**2 + dx**2)
    
    if FTfilt=='Gauss':
        gauss2dfilt = np.exp(-dx**2/(2*wnlim_x**2)-dy**2/(2*wnlim_y**2))
        filtered_shifted_f_uv = shifted_f_uv * gauss2dfilt
    else:
        raise ValueError('FTfilt=Gauss is only supported and recommended.')
        
    filtered_f_uv = np.fft.fftshift(filtered_shifted_f_uv)
    if remove_mean:
        filtered_f_uv[0,0]  = 0
        filtered_f_uv[0,-1] = 0
        filtered_f_uv[-1,0] = 0
        filtered_f_uv[-1,-1]= 0
    
    filtered_ps2d = 2 * np.log(np.absolute(filtered_f_uv))
    filtered_image = np.fft.ifft2(filtered_f_uv).real

    return filtered_image[:ny0, :nx0], filtered_f_uv, filtered_ps2d
    
    
def FT_grism_image(image, fftshift=False):
    
    ny0, nx0 = image.shape
    print('Shape of image: ', ny0, nx0)
    if np.mod(ny0, 2)==1 or np.mod(nx0, 2)==1:
        print('The shape should be (even, even). Add an emply column/row.')
             
    ny1 = ny0 + np.mod(ny0, 2)
    nx1 = nx0 + np.mod(nx0, 2)
    img_for_fft = np.zeros((ny1, nx1)) + np.nan
    img_for_fft[:ny0,:nx0] = image
    print('Shape of image for FFT: ', img_for_fft.shape)
        
    idx_nan = np.where(np.isnan(img_for_fft)|np.isinf(img_for_fft))
    idx_finite = np.where(np.isfinite(img_for_fft))
    
    if idx_nan[0].size>0:
        print('Replace {} NaN pixels: '.format(idx_nan[0].size)) 
        tmp = np.arange(idx_finite[0].size)
        np.random.seed(1)
        np.random.shuffle(tmp)
        img_for_fft[idx_nan] = img_for_fft[idx_finite][tmp[:idx_nan[0].size]]
        
    f_uv = np.fft.fft2(img_for_fft)

    if fftshift:
        f_uv = np.fft.fftshift(f_uv)
        
    ps2d = 2 * np.log(np.absolute(f_uv))

    return f_uv, ps2d, img_for_fft
    
    

def apply_mask_for_emlines():
    return 0


def save_spec2d_fits(filname,
                     list_of_images, 
                     waxis,
                     extname_list,
                     unit_list,
                     overwrite=False):
    hdu = fits.PrimaryHDU()
    hdul = [hdu]
    
    for i in range(len(list_of_images)):
        #print(i)
        hdr = fits.Header()
        hdr['EXTNAME']=extname_list[i]
        hdr['XTENSION']='IMAGE'
        hdr['BUNIT']=unit_list[i]
        hdr['NAXIS']=2
        hdr['NAXIS1']=list_of_images[i].shape[1]
        hdr['NAXIS2']=list_of_images[i].shape[0]
        hdr['CTYPE1']='WAVELENGTH'
        hdr['CUNIT1']='Angstrom'
        hdr['CRPIX1']=1.0
        hdr['CRVAL1']=waxis[0].to_value(u.angstrom)
        hdr['CDELT1']=(waxis[1]-waxis[0]).to_value(u.angstrom)
        hdu = fits.ImageHDU(list_of_images[i], header=hdr)
        hdul.append(hdu)
    hdul = fits.HDUList(hdul)
    hdul.writeto(filname, overwrite=overwrite)