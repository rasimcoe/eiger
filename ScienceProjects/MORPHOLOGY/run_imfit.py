import matplotlib
#matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import os
os.environ["WEBBPSF_PATH"] = "/scratch/mruari/EIGER/cache/webbpsf/webbpsf-data"
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.table import Table, Column
from astropy.wcs import WCS
from photutils import CircularAperture, aperture_photometry, CircularAnnulus, EllipticalAperture, EllipticalAnnulus, ApertureStats
from matplotlib.gridspec import GridSpec
from astropy.visualization import LogStretch,simple_norm
from astropy.stats import sigma_clipped_stats
import pyimfit
import webbpsf
from scipy.special import gammaincinv, gamma
import copy
from scipy.ndimage import rotate
from photutils.psf.matching import resize_psf
import pickle
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D

##########################################################
def get_model_JWST(FILTER, oversamp_fact=3, fov_arcsec=2, rot_angle=0.89):
    #get F356W PSF with extra fact 2 in oversamp
    instrument = webbpsf.NIRCam()
    instrument.filter=FILTER
    instrument.options['parity'] = 'odd'
    instrument.pixelscale = 0.03

    psf_model = instrument.calc_psf(oversample=oversamp_fact, fov_arcsec=fov_arcsec)

    psf_data_rot = rotate(psf_model['OVERSAMP'].data, rot_angle*u.deg, reshape=False)
    
    pix_clip = np.ceil(psf_model['OVERSAMP'].shape[0] * np.tan(rot_angle*u.deg)).astype(int)
    psf_data_rot_crop = psf_data_rot[pix_clip:-pix_clip,pix_clip:-pix_clip]

    if oversamp_fact == 1:
        return psf_data_rot_crop
    else:
        #trim array so final array is odd after resizing
        factor = np.floor(psf_data_rot_crop.shape[0]/oversamp_fact)
        trim = int((psf_data_rot_crop.shape[0]-(factor-(1 - (factor % 2)))*oversamp_fact)/2)
        psf_data_rot_crop_match = psf_data_rot_crop[trim:psf_data_rot_crop.shape[0]-trim,trim:psf_data_rot_crop.shape[1]-trim]
        det_samp = resize_psf(psf_data_rot_crop_match, 0.03/oversamp_fact, 0.03)
        return det_samp

def create_cutout(directimage,weightimage,RA,DEC,lx,ly):
    #directimage='reduced/test_10bright_mag26_z4to7/imaging_F356W/step_i2d.fits'
    hdu = fits.open(directimage)
    wcs = WCS(hdu['SCI'].header)
    positions=SkyCoord(RA,DEC,unit="deg")
    size = (ly,lx)
    cutout = Cutout2D(hdu['SCI'].data, position=positions, size=size, wcs=wcs)
    hdu['SCI'].data = cutout.data
    hdu['SCI'].header.update(cutout.wcs.to_header())
    HD=hdu['SCI'].header


    hdu = fits.open(weightimage)
    wcs = WCS(hdu['SCI'].header)
    positions=SkyCoord(RA,DEC,unit="deg")
    size = (ly,lx)
    cutout_wht = Cutout2D(hdu['SCI'].data, position=positions, size=size, wcs=wcs)

    return cutout.data,cutout_wht.data,HD

def galaxy_model(x0, y0, pa_disk, ell_disk, I_0, h):
    model = pyimfit.SimpleModelDescription()
    # define the limits on X0 and Y0 as +/-10 pixels relative to initial values
    model.x0.setValue(x0, [x0 - 10, x0 + 10])
    model.y0.setValue(y0, [y0 - 10, y0 + 10])

    disk = pyimfit.make_imfit_function('Exponential', label='disk')
    disk.PA.setValue(pa_disk, [0, 180])
    disk.ell.setValue(ell_disk, [0, 1])
    disk.I_0.setValue(I_0, [0, 10*I_0])
    disk.h.setValue(h, [0, 10*h])

    model.addFunction(disk)

    return model
##########################################################


CATALOG='/scratch/EIGER/extraction/J0100_photcat_v2_O3emitters_Systems_28092022.fits'
with fits.open(CATALOG) as hdul:
    orig_table = hdul[1].data
    orig_cols = orig_table.columns
    
rescale_noise=True #if true, rescales the mean(err_1d) to be equal to the std(data_1d) with some outlier removal


cat=fits.open(CATALOG)

data=cat[1].data
IDlist=data.field('NUMBER')
RAlist=data.field('RA_MANUAL')
DEClist=data.field('DEC_MANUAL')



psf_model = get_model_JWST('F356W', oversamp_fact=3, fov_arcsec=3, rot_angle=0.89)
print(psf_model)
st

mosaicname='/scratch/EIGER/identification/j0100_r.fits'
whtname='/scratch/EIGER/identification/j0100_r.fits'

cutout_data,wht_data,hd=create_cutout(mosaicname,whtname,RAlist[1],DEClist[1],100,100)


model_desc = galaxy_model(x0=50, y0=50, pa_disk=90.0, ell_disk=0.5, I_0=1, h=25)
imfit_fitter = pyimfit.Imfit(model_desc)

imfit_fitter.loadData(cutout_data,error=1./np.sqrt(wht_data))#, error=1./np.sqrt(whtcut), mask=total_mask)
result = imfit_fitter.doFit()#solver='DE')
img_model = imfit_fitter.getModelImage(newParameters=result.params)


print(result)
fits.writeto('cutout.fits',cutout_data,overwrite=True)
fits.writeto('model.fits',img_model,overwrite=True)
fits.writeto('residual.fits',cutout_data-img_model,overwrite=True)


stop
basedir = '/scratch/mruari/EIGER/imaging/J0100+2802/'
photdir = basedir + 'photometry/sextractor/'
workdir = basedir + 'junk/foreground_sub/'
os.chdir(workdir)

#make backup directory
bkupdir = workdir + 'results_backup/'
if not os.path.isdir(bkupdir): os.mkdir(bkupdir)

targ_cat = Table.read(workdir+'elliptical_catalog.fits')

catname = photdir + 'sourcecat_det_F356W.fits'
catalog = Table.read(catname)

inital_solver = 'DE'
twocomp_solver = 'DE'
inital_usePSF = True
save_name = inital_solver+twocomp_solver+'_withPSF'

#FILTER = 'F115W' 'F356W',
for FILTER in ['F200W','F115W']:

    mosaicname = basedir +'junk/current_best/stack_'+FILTER+'_pipe4_v2_fluxcal_20220831.fits'
    mosaic = fits.open(mosaicname)
    hd=mosaic['SCI'].header
    data=mosaic['SCI'].data
    weight=mosaic['WHT'].data
    wcs = WCS(hd)

    #psf_model = get_model_JWST(FILTER, oversamp_fact=1, fov_arcsec=0.63)
    psf_model = get_model_JWST(FILTER, oversamp_fact=3, fov_arcsec=2, rot_angle=0.89)

    img_sub = np.copy(data)
    img_mod = np.zeros(data.shape, dtype=float)
    save_params =  Table(names=('X01', 'Y01', 'c'), dtype=('f4', 'i4', 'S2'))

    for iell, elliptical in enumerate(targ_cat):
        pkl_intial  = bkupdir+'imfit_result_id%i_%s_'%(elliptical['NUMBER'], FILTER)+save_name+'_inital.pkl'
        pkl_twocomp = bkupdir+'imfit_result_id%i_%s_'%(elliptical['NUMBER'], FILTER)+save_name+'_twocomp.pkl'
        
        x_center = np.round(elliptical['X_IMAGE']).astype(int)-1
        y_center = np.round(elliptical['Y_IMAGE']).astype(int)-1

        boarder = 600
        xmin, xmax = x_center-300,x_center+300
        ymin, ymax = y_center-300,y_center+300
        cutout = np.copy(data[ymin:ymax,xmin:xmax])
        whtcut = np.copy(weight[ymin:ymax,xmin:xmax])
        cat_cut = copy.deepcopy(catalog[(catalog['X_IMAGE'] > xmin)&(catalog['X_IMAGE'] < xmax)&(catalog['Y_IMAGE'] > ymin)&(catalog['Y_IMAGE'] < ymax)])
        cat_cut['X_IMAGE'] -= xmin
        cat_cut['Y_IMAGE'] -= ymin

        #mask faint objects
        #elliptical = cat_cut[(np.argmin(cat_cut['MAG_AUTO']))]
        cat_cut = cat_cut[(cat_cut['MAG_AUTO'] < 28)&(cat_cut['MAG_AUTO'] > elliptical['MAG_AUTO'])]
        total_mask = np.zeros(cutout.shape, dtype=bool)
        for io, obj in enumerate(cat_cut):
            posn = np.array([obj['X_IMAGE']-1, obj['Y_IMAGE']-1])
            a_image = (obj['A_IMAGE']*obj['KRON_RADIUS']/1.5)
            b_image = (obj['B_IMAGE']*obj['KRON_RADIUS']/1.5 )
            theta_image = obj['THETA_IMAGE'] * u.deg
            apers = EllipticalAperture(posn.T, a_image, b_image, theta=theta_image)
            mask = apers.to_mask(method='center')
            obj_mask = mask.to_image(cutout.shape)
            total_mask = (total_mask)|(obj_mask == 1.0)


                
        n_guess = 4.0
        r_eff_guess = elliptical['B_IMAGE']*elliptical['KRON_RADIUS']/2.5
        bn = gammaincinv(2. * n_guess, 0.5)
        I_e_guess = elliptical['FLUX_AUTO']/(2.0*np.pi* n_guess * np.exp(bn)/(bn**(2.*n_guess)) * gamma(2.0*n_guess) *r_eff_guess**2)
        #mu_e_guess = elliptical['MAG_AUTO'] - 28.03 + 5.0*np.log10(r_eff_guess) + 2.5 * np.log10(2.0*np.pi* n_guess * np.exp(bn)/(bn**(2.*n_guess)) * gamma(2.0*n_guess))

        initalParamsDict = {'PA':  [130, 0,180],                               #elliptical['THETA_IMAGE']
                            'ell': [elliptical['ELLIPTICITY'], 0.0,1.0], 
                            'n':   [n_guess, 0.5,6.0],
                            'I_e': [I_e_guess, I_e_guess/10., I_e_guess*5], 
                            'r_e': [r_eff_guess,  5, 200.]}

        funcDict = {'name': "Sersic", 'label': "inital", 'parameters': initalParamsDict}

        functionSetDict = {'label':'v1',    
                            'X0': [elliptical['X_IMAGE']-xmin, elliptical['X_IMAGE']-25-xmin,elliptical['X_IMAGE']+25-xmin], 
                            'Y0': [elliptical['Y_IMAGE']-ymin, elliptical['Y_IMAGE']-25-ymin,elliptical['Y_IMAGE']+25-ymin], 
                            'function_list': [funcDict]}

        modelDict = {'function_sets': [functionSetDict]}
        model_desc = pyimfit.ModelDescription.dict_to_ModelDescription(modelDict)

        if inital_usePSF:
            fitter_psf = pyimfit.Imfit(model_desc, psf_model)
        else:
            fitter_psf = pyimfit.Imfit(model_desc)

   
        if os.path.isfile(pkl_intial):
            result_psf = pickle.load(open(pkl_intial, "rb"))
        else:
            fitter_psf.loadData(cutout, error=1./np.sqrt(whtcut), mask=total_mask)
            result_psf = fitter_psf.doFit(solver=inital_solver)
             
            pickle.dump(result_psf, open(pkl_intial, "wb"))

        img_model_psf = fitter_psf.getModelImage(newParameters=result_psf.params)

        # fig=plt.figure(figsize=(17.1, 12))
        # gs = GridSpec(1,2,hspace=0.25,wspace=0.2, left=0.02, bottom=0.06,right=0.98, top=0.95)

        # plt.subplot(gs[0])
        # norm = simple_norm(cutout-img_model_psf, 'linear', percent=98.)
        # plt.imshow(cutout, norm=norm, origin='lower', cmap='viridis')
        # plt.title('Before')
        # plt.subplot(gs[1])
        # plt.imshow(cutout-img_model_psf, norm=norm, origin='lower', cmap='viridis')
        # plt.title('with psf')
        # plt.show()


        inital_dict = fitter_psf.getModelAsDict()['function_sets'][0]
        param_dict =  inital_dict['function_list'][0]['parameters']

        diskParamsDict = {'PA': [param_dict['PA'][0], 0,180], 
                            'ell': [param_dict['ell'][0], 0.0,0.6], 
                            'n': [param_dict['n'][0], 0.5,4.5],
                            'I_e': [param_dict['I_e'][0], param_dict['I_e'][0]*0.1,param_dict['I_e'][0]*5], 
                            'r_e': [param_dict['r_e'][0], 1,60]}
        diskDict = {'name': "Sersic", 'label': "disk", 'parameters': diskParamsDict}
        diskfunctionSetDict = {'label':'ds','X0': [inital_dict['X0'][0], inital_dict['X0'][0]-1,inital_dict['X0'][0]+1], 
                               'Y0': [inital_dict['Y0'][0], inital_dict['Y0'][0]-1,inital_dict['Y0'][0]+1], 'function_list': [diskDict]}

        bulgeParamsDict = {'PA': [param_dict['PA'][0], 0,180], 
                            'ell': [0.6, 0.5,1.0], 
                            'n': [4.5, 4.0,8.0],
                            'I_e': [param_dict['I_e'][0], param_dict['I_e'][0]*0.1,param_dict['I_e'][0]*5], 
                            'r_e': [20, 1,30]}

        bulgeDict = {'name': "Sersic", 'label': "bulge", 'parameters': bulgeParamsDict}
        bulgefunctionSetDict = {'label':'bg','X0': [inital_dict['X0'][0], inital_dict['X0'][0]-3,inital_dict['X0'][0]+3], 
                                'Y0': [inital_dict['Y0'][0], inital_dict['Y0'][0]-3,inital_dict['Y0'][0]+3], 'function_list': [bulgeDict]}

        modelDict = {'function_sets': [diskfunctionSetDict, bulgefunctionSetDict]}
        model_desc = pyimfit.ModelDescription.dict_to_ModelDescription(modelDict)

        fitter_twocomp = pyimfit.Imfit(model_desc, psf_model)
        
        if os.path.isfile(pkl_twocomp):
            result2 = pickle.load(open(pkl_twocomp, "rb"))
        else:
            fitter_twocomp.loadData(cutout, error=1./np.sqrt(whtcut), mask=total_mask)
            result2 = fitter_twocomp.doFit(solver=twocomp_solver)
             
            pickle.dump(result2, open(pkl_twocomp, "wb"))
        
        img_model2 = fitter_twocomp.getModelImage(newParameters=result2.params)

        fig=plt.figure(figsize=(15, 8))
        gs = GridSpec(1,3,hspace=0.25,wspace=0.2, left=0.08, bottom=0.06,right=0.98, top=0.95)

        plt.subplot(gs[0])
        norm = simple_norm(cutout-img_model_psf, 'linear', percent=99.5)
        plt.imshow(cutout, norm=norm, origin='lower', cmap='viridis')
        plt.title('Before')
        plt.subplot(gs[1])
        plt.imshow(cutout-img_model_psf, norm=norm, origin='lower', cmap='viridis')
        plt.title('Inital')
        plt.subplot(gs[2])
        plt.imshow(cutout-img_model2, norm=norm, origin='lower', cmap='viridis')
        plt.title('2 Comp')
        plt.savefig(workdir+'imfit_result_id%i_%s_'%(elliptical['NUMBER'], FILTER)+save_name+'.pdf')
        plt.close()
        
            
        #new model to generate full model
        final_twocomp = pyimfit.Imfit(fitter_twocomp.getModelDescription(), psf_model)
        final_params = copy.deepcopy(result2.params)
        final_params[0] = final_params[0] + xmin
        final_params[1] = final_params[1] + ymin
        final_params[7] = final_params[7] + xmin
        final_params[8] = final_params[8] + ymin
        full_model = final_twocomp.getModelImage(newParameters=final_params, shape=img_sub.shape)

        if iell ==0:
            fit_sum = Table(np.hstack([np.array([elliptical['NUMBER']]),final_params]), 
                names=['NUMBER','X01','Y01','PA1','ell1','n1','I_e1','r_e1','X02','Y02','PA2','ell2','n2','I_e2','r_e2'])
        else:
            fit_sum.add_row(np.hstack([np.array([elliptical['NUMBER']]),final_params]))
        #sub main image
        img_sub = img_sub - full_model
        img_mod = img_mod + full_model

    fits.writeto( workdir+'stack_subbed_'+FILTER+'_'+save_name+'.fits' , img_sub, header=hd, overwrite=True)  
    fits.writeto( workdir+'models_'      +FILTER+'_'+save_name+'.fits'  , img_mod, header=hd, overwrite=True)
    fit_sum.write(workdir+'model_params_'+FILTER+'_'+save_name+'.fits' , overwrite=True)

#try basic fit with lens masked

# image_sub = cutout-img_model2                                                                                                        
# from scipy import ndimage, misc
# x_center = np.round(result_psf.params[0]).astype(int)-1
# y_center = np.round(result_psf.params[1]).astype(int)-1
# core = np.copy(image_sub[y_center-10:y_center+10+1,x_center-10:x_center+10+1])
# filtered = ndimage.median_filter(core, size=2)

# image_filtered = np.copy(image_sub)
# image_filtered[y_center-10:y_center+10+1,x_center-10:x_center+10+1] = filtered


# fig=plt.figure(figsize=(17.1, 12))
# gs = GridSpec(1,2,hspace=0.25,wspace=0.2, left=0.02, bottom=0.06,right=0.98, top=0.95)

# plt.subplot(gs[0])
# norm = simple_norm(image_filtered, 'linear', percent=100)
# plt.imshow(image_sub, norm=norm, origin='lower', cmap='viridis')
# plt.title('Before')
# plt.subplot(gs[1])
# plt.imshow(image_filtered, norm=norm, origin='lower', cmap='viridis')
# plt.title('2 Comp')
# plt.show()
