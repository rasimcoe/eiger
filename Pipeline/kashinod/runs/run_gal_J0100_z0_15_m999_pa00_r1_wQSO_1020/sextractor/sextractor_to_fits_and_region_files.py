#!/usr/bin/env python
# coding: utf-8

# In[1]:


import astropy
from astropy.io.ascii import SExtractor
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs import utils
import numpy as np
import sys


# In[2]:


fits_fil_array = [
    '../calibrated_img3_wfss_medSubt41/jw01243_nrc_img3_wfss_visit1a_i2d.fits',
    '../calibrated_img3_wfss_medSubt41/jw01243_nrc_img3_wfss_visit1b_i2d.fits',
    '../calibrated_img3_wfss_medSubt41/jw01243_nrc_img3_wfss_visit2a_i2d.fits',
    '../calibrated_img3_wfss_medSubt41/jw01243_nrc_img3_wfss_visit2b_i2d.fits',
    '../calibrated_img3_wfss_medSubt41/jw01243_nrc_img3_wfss_visit3a_i2d.fits',
    '../calibrated_img3_wfss_medSubt41/jw01243_nrc_img3_wfss_visit3b_i2d.fits',
    '../calibrated_img3_wfss_medSubt41/jw01243_nrc_img3_wfss_visit4a_i2d.fits',
    '../calibrated_img3_wfss_medSubt41/jw01243_nrc_img3_wfss_visit4b_i2d.fits',
    ]
    
# Sextractor results
cat_name_array = [
    'results_calibrated_img3_wfss_medSubt41/sex_jw01243_nrc_img3_wfss_visit1a_i2d',
    'results_calibrated_img3_wfss_medSubt41/sex_jw01243_nrc_img3_wfss_visit1b_i2d',
    'results_calibrated_img3_wfss_medSubt41/sex_jw01243_nrc_img3_wfss_visit2a_i2d',
    'results_calibrated_img3_wfss_medSubt41/sex_jw01243_nrc_img3_wfss_visit2b_i2d',
    'results_calibrated_img3_wfss_medSubt41/sex_jw01243_nrc_img3_wfss_visit3a_i2d',
    'results_calibrated_img3_wfss_medSubt41/sex_jw01243_nrc_img3_wfss_visit3b_i2d',
    'results_calibrated_img3_wfss_medSubt41/sex_jw01243_nrc_img3_wfss_visit4a_i2d',
    'results_calibrated_img3_wfss_medSubt41/sex_jw01243_nrc_img3_wfss_visit4b_i2d',
    ]


# In[3]:

i_fil=int(sys.argv[1])

fits_fil=fits_fil_array[i_fil]
cat_name=cat_name_array[i_fil]
print('Fits file', fits_fil)
print('Cat name: ', cat_name)


# In[4]:


# Need to take wcs from the fits file

hdul = fits.open(fits_fil)
wcs = WCS(hdul[1].header)
hdul.close()


# In[5]:


# Pixel scale
pixel_area = utils.proj_plane_pixel_area(wcs)
pixel_scales = utils.proj_plane_pixel_scales(wcs)
print('Pixel scale: ', pixel_scales, ' deg/pixel')
print('Pixel area: ', pixel_area, ' sqdeg/pixel')
print('Pixel scale: ', pixel_scales*3600., ' arcsec/pixel')
print('Pixel area: ', pixel_area*3600.*3600., ' arcsec2/pixel')
print('Pixel area: ', pixel_area*(np.pi/180.)**2, 'stradian/pixel')


# In[6]:


# MAG_ZEROPOINT
MAG_ZEROPOINT = 8.90 - 2.5*np.log10(pixel_area*(np.pi/180.)**2*1e6)
print('MAG_ZEROPOINT: ', MAG_ZEROPOINT)


# In[7]:


# Sextractor catalog
cat_fil = cat_name+'.cat'
#sex = SExtractor()
cat = SExtractor().read(cat_fil)


# In[8]:


# Correct world coordinates using wcs
world = wcs.pixel_to_world(cat['X_IMAGE']-1,cat['Y_IMAGE']-1)
peak_world = wcs.pixel_to_world(cat['XPEAK_IMAGE']-1,cat['YPEAK_IMAGE']-1)


# In[9]:


cat['ALPHA_J2000']=world.ra.degree
cat['DELTA_J2000']=world.dec.degree
cat['ALPHAPEAK_J2000']=peak_world.ra.degree
cat['DELTAPEAK_J2000']=peak_world.dec.degree


# In[10]:


# Update magnitudes and fluxes using MAG_ZEROPOINT
# Physical units 'MJy/s'
cat['MAG_ISO'] = cat['MAG_ISO'] + MAG_ZEROPOINT
cat['MAG_ISOCOR'] = cat['MAG_ISOCOR'] + MAG_ZEROPOINT
cat['MAG_APER'] = cat['MAG_APER'] + MAG_ZEROPOINT
cat['MAG_AUTO'] = cat['MAG_AUTO'] + MAG_ZEROPOINT
cat['MAG_BEST'] = cat['MAG_BEST'] + MAG_ZEROPOINT


# In[ ]:





# In[11]:


cat


# In[12]:


cat.write(cat_name+'.fits', overwrite=True)


# In[13]:


# region
from astropy import units as u
from astropy.coordinates import SkyCoord
import regions 
from regions import PixCoord, EllipseSkyRegion, EllipsePixelRegion, TextSkyRegion


# In[14]:


skyCoord = SkyCoord(ra=cat['ALPHA_J2000'], dec=cat['DELTA_J2000'], unit='deg')


# In[15]:


ellipse_sky = []
ellipse_sky_2 = []
text_sky = []
for i in range(len(cat)):
    tmp = EllipseSkyRegion(center=skyCoord[i],
                           height=cat['A_WORLD'][i]*cat['KRON_RADIUS'][i]*u.degree,
                           width=cat['B_WORLD'][i]*cat['KRON_RADIUS'][i]*u.degree,
                           angle=(cat['THETA_WORLD'][i]+90.)*u.degree)
    ellipse_sky.append(tmp)
    tmp = EllipseSkyRegion(center=skyCoord[i],
                           height=cat['A_WORLD'][i]*cat['KRON_RADIUS'][i]*1.5*u.degree,
                           width=cat['B_WORLD'][i]*cat['KRON_RADIUS'][i]*1.5*u.degree,
                           angle=(cat['THETA_WORLD'][i]+90.)*u.degree)
    ellipse_sky_2.append(tmp)


# In[16]:


with open(cat_name+'.reg', 'w') as f:
    f.write('# Region file format: DS9 astropy/regions\n')
    f.write('icrs\n')
    for i in range(len(cat)):
        reg_str = regions.ds9_objects_to_string([ellipse_sky[i]], coordsys='icrs')
        reg_str_2 = regions.ds9_objects_to_string([ellipse_sky_2[i]], coordsys='icrs')
        #line = reg_str.split('\n')[-2] + ' \n'+\
        #    reg_str_2.split('\n')[-2] + ' # text={'+str(cat['NUMBER'][i])+'}\n'
        line = reg_str_2.split('\n')[-2]+' # text={'+str(cat['NUMBER'][i])+'}\n'
        f.write(line)


# In[ ]:





# In[ ]:





# In[ ]:




