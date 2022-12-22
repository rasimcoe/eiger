"""
   This code will take any photometric object ID and then plot the best fit SED from eazy and the chi-square distribution.
   it will also show some potential solutions for Halpha and OIII emission line objects for F356W grism spectroscopy.

   Dependencies= requires EAZY to be installed.
                 requires EAZY output files.


"""

import eazy
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.gridspec import GridSpec
from collections import OrderedDict
import astropy.units as u
import eazy
import eazy.hdf5
import eazy.igm
import scipy.interpolate 
from scipy.integrate import cumtrapz
from glob import glob
from astropy.visualization import ZScaleInterval
from matplotlib.gridspec import GridSpec

tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)  

def show_fits(self, ix, ax, zshow=None, ndraws=100, fitter='nnls',draws = None,show_fnu = True,template_color='#1f77b4',
    show_upperlimits = True, show_components = True, show_redshift_draws = False, show_missing = False, 
    snr_thresh=2.0, draws_cmap = None, show_temp_err=False, show_all_temps=False, ylim_man=None, xlim=[0.4, 8], plot_filters=True):

    TEMPLATE_REDSHIFT_TYPE = 'nearest'
    IGM_OBJECT = eazy.igm.Inoue14()
    axes = None
    add_label=False
    maglim = None
    with_tef=True
    #load object data
    z = self.zbest[ix]      #optionally set to show redshift
    fnu_i = self.fnu[ix, :]
    efnu_i = self.efnu[ix,:]
    zspec_i = self.ZSPEC[ix]
    ra_i = self.RA[ix]
    dec_i = self.DEC[ix]
    lnp_i = self.lnp[ix,:]
    log_prior_i = self.full_logprior[ix,:].flatten()
    chi2_i = self.chi2_fit[ix,:]
    zgrid = self.zgrid
    ok_i = self.ok_data[ix,:]

    if zshow is not None:
        z = zshow

    chi2 = np.squeeze(chi2_i)
    prior = np.exp(log_prior_i)
    pz = np.exp(-(chi2-chi2.min())/2.) #*prior
    pz /= np.trapz(pz, zgrid)
    #pz = np.exp(lnp_i).flatten()

    ## SED        
    fnu_i = np.squeeze(fnu_i)*self.ext_redden*self.zp
    efnu_i = np.squeeze(efnu_i)*self.ext_redden*self.zp
    ok_band = (fnu_i/self.zp > self.param['NOT_OBS_THRESHOLD']) 
    ok_band &= (efnu_i/self.zp > 0)
    efnu_i[~ok_band] = self.param['NOT_OBS_THRESHOLD'] - 9.

    ## Evaluate coeffs at specified redshift
    tef_i = self.TEF(z)
    A = np.squeeze(self.tempfilt(z))
    chi2_i, coeffs_i, fmodel, draws = eazy.photoz.template_lsq(fnu_i, efnu_i, A, 
                                               tef_i, self.zp, 
                                               ndraws, fitter)
    if draws is None:
        efmodel = 0
    else:
        efmodel = np.percentile(np.dot(draws, A), [16,84], axis=0)
        efmodel = np.squeeze(np.diff(efmodel, axis=0)/2.)
        
    ## Full SED
    templ = self.templates[0]
    tempflux = np.zeros((self.NTEMP, templ.wave.shape[0]),
                        dtype=self.ARRAY_DTYPE)
    for i in range(self.NTEMP):
        zargs = {'z':z, 'redshift_type':TEMPLATE_REDSHIFT_TYPE}
        fnu = self.templates[i].flux_fnu(**zargs)*self.tempfilt.scale[i]
        try:
            tempflux[i, :] = fnu
        except:
            tempflux[i, :] = np.interp(templ.wave,
                                       self.templates[i].wave, fnu)
            
    templz = templ.wave*(1+z)

    if self.tempfilt.add_igm:
        igmz = templ.wave*0.+1
        lyman = templ.wave < 1300
        igmz[lyman] = IGM_OBJECT.full_IGM(z, templz[lyman])
    else:
        igmz = 1.

    templf = np.dot(coeffs_i, tempflux)*igmz
    if draws is not None:
        templf_draws = np.dot(draws, tempflux)*igmz
            
    fnu_factor = 10**(-0.4*(self.param['PRIOR_ABZP']+48.6))

    if show_fnu:
        if show_fnu == 2:
            templz_power = -1
            flam_spec = 1.e29/(templz/1.e4)
            flam_sed = 1.e29/self.ext_corr/(self.pivot/1.e4)
            ylabel = (r'$f_\nu / \lambda$ [$\mu$Jy / $\mu$m]')
            flux_unit = u.uJy / u.micron
        else:
            templz_power = 0
            flam_spec = 1.e29
            flam_sed = 1.e29
            ylabel = (r'$f_\nu$ [$\mu$Jy]')    
            flux_unit = u.uJy
        
    else:
        templz_power = -2
        flam_spec = utils.CLIGHT*1.e10/templz**2/1.e-19
        flam_sed = utils.CLIGHT*1.e10/self.pivot**2/self.ext_corr/1.e-19
        ylabel = (r'$f_\lambda [10^{-19}$ erg/s/cm$^2$]')
        
        flux_unit = 1.e-19*u.erg/u.s/u.cm**2/u.AA
                    

    data = OrderedDict(ix=ix, id=self.OBJID[ix], z=z,
                   z_spec=zspec_i, 
                   pivot=self.pivot, 
                   model=fmodel*fnu_factor*flam_sed,
                   emodel=efmodel*fnu_factor*flam_sed,
                   fobs=fnu_i*fnu_factor*flam_sed, 
                   efobs=efnu_i*fnu_factor*flam_sed,
                   valid=ok_i,
                   tef=tef_i,
                   templz=templz,
                   templf=templf*fnu_factor*flam_spec,
                   show_fnu=show_fnu*1,
                   flux_unit=flux_unit,
                   wave_unit=u.AA, 
                   chi2=chi2_i, 
                   coeffs=coeffs_i)


    ###### Make the plot

    #fig = plt.figure()
    #fig_axes = GridSpec(1,1,width_ratios=[1])
    #ax =  fig.add_subplot(fig_axes[0])
                    
    ax.scatter(self.pivot/1.e4, fmodel*fnu_factor*flam_sed, 
               color='w', label=None, zorder=1, s=120, marker='o')

    ax.scatter(self.pivot/1.e4, fmodel*fnu_factor*flam_sed, marker='o',
              color=template_color, label=None, zorder=2, s=50, 
              alpha=0.8)

    if draws is not None:
        ax.errorbar(self.pivot/1.e4, fmodel*fnu_factor*flam_sed,
                    efmodel*fnu_factor*flam_sed, alpha=0.8,
                    color=template_color, zorder=2,
                    marker='None', linestyle='None', label=None)

    # Missing data
    missing = (fnu_i < self.param['NOT_OBS_THRESHOLD']) 
    missing |= (efnu_i < 0)

    # Detection
    sn2_detection = (~missing) & (fnu_i/efnu_i > snr_thresh)

    # S/N < 2
    sn2_not = (~missing) & (fnu_i/efnu_i <= snr_thresh)

    # Uncertainty with TEF
    if with_tef:
        err_tef = np.sqrt(efnu_i**2+(tef_i*fnu_i)**2)            
    else:
        err_tef = efnu_i*1
        
    ax.errorbar(self.pivot[sn2_detection]/1.e4, 
                (fnu_i*fnu_factor*flam_sed)[sn2_detection], 
                (err_tef*fnu_factor*flam_sed)[sn2_detection], 
                color='k', marker='s', linestyle='None', label=None, 
                zorder=10)

    if show_upperlimits:
        ax.errorbar(self.pivot[sn2_not]/1.e4, 
                    (fnu_i*fnu_factor*flam_sed)[sn2_not], 
                    (efnu_i*fnu_factor*flam_sed)[sn2_not], color='k', 
                    marker='s', alpha=0.4, linestyle='None', label=None)

    if show_missing:
        ax.errorbar(self.pivot[missing]/1.e4, 
                    (fnu_i*fnu_factor*flam_sed)[missing]*0, 
                    (efnu_i*fnu_factor*flam_sed)[missing], 
                    color='0.7', marker='x', linestyle='None', 
                    alpha=0.4, label=None)

    pl = ax.plot(templz/1.e4, templf*fnu_factor*flam_spec, alpha=0.5, 
                 zorder=-1, color=template_color, 
                 label='z={0:.2f}'.format(z))

    if show_all_temps:
        colors = ['#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
                  '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        for i in range(self.NTEMP):
            i_templz = np.argmin(np.abs(templz - self.pivot[-2]))
            obs_flux = (fnu_i*fnu_factor*flam_sed)[-2]
            norm = fnu_i[-2]/(tempflux[i,i_templz]*igmz[i_templz])
            pi = ax.plot(templz/1.e4, 
                norm*tempflux[i,:]*igmz*fnu_factor*flam_spec, 
                      alpha=0.5, zorder=-1, 
                      label=self.templates[i].name.split('.dat')[0], 
                      color=colors[i % len(colors)])

    if show_components:
        colors = ['#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
                  '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        for i in range(self.NTEMP):
            if coeffs_i[i] != 0:
                pi = ax.plot(templz/1.e4, 
                    coeffs_i[i]*tempflux[i,:]*igmz*fnu_factor*flam_spec, 
                          alpha=0.5, zorder=-1, 
                          label=self.templates[i].name.split('.dat')[0], 
                          color=colors[i % len(colors)])
                

                        
    elif show_redshift_draws:
        
        if draws_cmap is None:
            draws_cmap = plt.cm.rainbow
            
        # Draw random values from p(z)
        pz = np.exp(lnp_i).flatten()
        pzcum = cumtrapz(pz, x=self.zgrid)
        
        if show_redshift_draws == 1:
            nzdraw = 100
        else:
            nzdraw = show_redshift_draws*1
        
        rvs = np.random.rand(nzdraw)
        zdraws = np.interp(rvs, pzcum, self.zgrid[1:])
        
        for zi in zdraws:
            Az = np.squeeze(self.tempfilt(zi))
            chi2_zi, coeffs_zi, fmodelz, __ = eazy.photoz.template_lsq(fnu_i, efnu_i, 
                                                   Az, 
                                                   self.TEF(zi), self.zp, 
                                                   0, fitter)
                                                   
            c_i = np.interp(zi, self.zgrid, np.arange(self.NZ)/self.NZ)
            
            templzi = templ.wave*(1+zi)
            if self.tempfilt.add_igm:
                igmz = templ.wave*0.+1
                lyman = templ.wave < 1300
                igmz[lyman] = IGM_OBJECT.full_IGM(zi, templzi[lyman])
            else:
                igmz = 1.

            templfz = np.dot(coeffs_zi, tempflux)*igmz                
            templfz *=  flam_spec * (templz / templzi)**templz_power
            
            plz = ax.plot(templzi/1.e4, templfz*fnu_factor,
                         alpha= 0.25, #np.maximum(0.1, 1./nzdraw), 
                         zorder=-1, color=draws_cmap(c_i))
    #print(draws)      
    if draws is not None:
        if show_temp_err:
            templf_width = np.percentile(templf_draws*fnu_factor*flam_spec, 
                                         [16,84], axis=0)
            ax.fill_between(templz/1.e4, templf_width[0,:], templf_width[1,:], 
                            color=pl[0].get_color(), alpha=0.1, label=None)
                                                   

                   
    ax.set_ylabel(ylabel)

    if sn2_detection.sum() > 0:
        ymax = (fmodel*fnu_factor*flam_sed)[sn2_detection].max()
    else:
        ymax = (fmodel*fnu_factor*flam_sed).max()
                
    if np.isfinite(ymax):
        ax.set_ylim(-0.1*ymax, 1.2*ymax)
        print(-0.1*ymax, 1.2*ymax)
    if ylim_man is not None:
        ax.set_ylim(ylim_man)

    ax.set_xlim(xlim)
    xt = np.array([0.1, 0.5, 1, 2, 4, 8, 24, 160, 500])*1.e4

    ax.semilogx()

    valid_ticks = (xt > xlim[0]*1.e4) & (xt < xlim[1]*1.e4)
    if valid_ticks.sum() > 0:
        xt = xt[valid_ticks]
        ax.set_xticks(xt/1.e4)
        ax.set_xticklabels(xt/1.e4)

    ax.set_xlabel(r'$\lambda_\mathrm{obs}$')
    ax.grid()

    if add_label:
        txt = '{0}\nID={1}'
        txt = txt.format(self.param['MAIN_OUTPUT_FILE'], 
                         self.OBJID[ix]) #, self.prior_mag_cat[ix])
                         
        ax.text(0.95, 0.95, txt, ha='right', va='top', fontsize=7,
                transform=ax.transAxes, 
                bbox=dict(facecolor='w', alpha=0.5), zorder=10)
        
        ax.legend(fontsize=7, loc='upper left')

    # Optional mag scaling if show_fnu = 1 for uJy
    if (maglim is not None) & (show_fnu == 1):
        
        ax.semilogy()
        # Limits
        ax.scatter(self.pivot[sn2_not]/1.e4,
                   ((3*efnu_i)*fnu_factor*flam_sed)[sn2_not], 
                   color='k', marker='v', alpha=0.4, label=None)
        
        # Mag axes
        axm = ax.twinx()
        ax.set_ylim(10**(-0.4*(np.array(maglim)-23.9)))
        axm.set_ylim(0,1)
        ytv = np.arange(maglim[0], maglim[1], -1, dtype=int)
        axm.set_yticks(np.interp(ytv, maglim[::-1], [1,0]))
        axm.set_yticklabels(ytv)

    if plot_filters:
        ax.set_ylim(-0.3*ymax, 1.2*ymax)
        for ii, filt in enumerate(self.filters):
            print('plotting filter',filt)
            ax.plot(filt.wave/1e4, (filt.throughput-1.0)*ymax*0.3, c='grey')

    #plt.show()

def chi2_dist(self, ix, ax):

    #load object data
    z = self.zbest[ix]      #optionally set to show redshift
    fnu_i = self.fnu[ix, :]
    efnu_i = self.efnu[ix,:]
    zspec_i = self.ZSPEC[ix]
    ra_i = self.RA[ix]
    dec_i = self.DEC[ix]
    lnp_i = self.lnp[ix,:]
    log_prior_i = self.full_logprior[ix,:].flatten()
    chi2_i = self.chi2_fit[ix,:]
    zgrid = self.zgrid
    ok_i = self.ok_data[ix,:]


    z_ml = self.zml[ix]
    self.Rz = self.compute_full_risk()
    z_min_risk = self.zgrid[np.argmin(self.Rz[ix,:])] 
    z_raw_chi2 = self.zchi2[ix]


    chi2 = np.squeeze(chi2_i)
    prior = np.exp(log_prior_i)
    pz = np.exp(-(chi2-chi2.min())/2.) #*prior
    pz /= np.trapz(pz, zgrid)
    #pz = np.exp(lnp_i).flatten()   

    chi2_min = np.nanmin(chi2[(chi2 > 0.)])

    #fig = plt.figure()
    #fig_axes = GridSpec(1,1,width_ratios=[1])
    #ax =  fig.add_subplot(fig_axes[0])
    ax.plot(zgrid, chi2, color=tableau20[2], label=None)
    #ax.plot(self.zgrid, prior/prior.max()*pz.max(), color='g',                label='prior')
    ax.fill_between(zgrid, chi2, chi2*0, color='yellow', alpha=0.5, label=None)

    #ax.vlines(zspec_i, chi2_min*0.9, chi2.max()*1.05, color=tableau20[12], label='zsp={0:.3f}'.format(zspec_i))

    ax.vlines(z_ml,       chi2_min*0.9, chi2.max()*1.05, color=tableau20[16], label='z_ML={0:.3f}'.format(z_ml))
    ax.vlines(z_min_risk, chi2_min*0.9, chi2.max()*1.05, color=tableau20[8], label='z_min_risk={0:.3f}'.format(z_min_risk))
    ax.vlines(z_raw_chi2, chi2_min*0.9, chi2.max()*1.05, color=tableau20[4], ls=':', label='z_raw_chi2={0:.3f}'.format(z_raw_chi2))
        
    plt.yscale("log")
    ax.set_ylim(chi2_min*0.9,chi2.max()*1.05)
        
    ax.set_xlim(0,zgrid[-1])

    #plot [OII] range
    y_OIII = 10.**(0.5*(np.log10(chi2_min)+np.log10(chi2.max())))
    col_OIII = tableau20[0]
    plt.text(6.15,y_OIII*1.1,'[OIII]',color=col_OIII, ha='center')
    plt.plot([5.3,7.0], [y_OIII,y_OIII], c=col_OIII)
    plt.plot([5.3,5.3], [y_OIII*0.9,y_OIII*1.1], c=col_OIII)
    plt.plot([7.0,7.0], [y_OIII*0.9,y_OIII*1.1], c=col_OIII)
    #plot best [OIII] z
    z_OIII_best = zgrid[(zgrid > 5.3)&(zgrid < 7.0)][((chi2[(zgrid > 5.3)&(zgrid < 7.0)]).argmin())]
    chi2_OIII_best = chi2[(zgrid > 5.3)&(zgrid < 7.0)][((chi2[(zgrid > 5.3)&(zgrid < 7.0)]).argmin())]
    plt.scatter(z_OIII_best, chi2_OIII_best, c=col_OIII,marker='o',ec=tableau20[0], label=r'[OIII] best $\Delta\chi^2=$%0.2f' % (chi2_OIII_best - chi2_min))
    #plt.errorbar(6.15, 3.0, xerr=0.85,capsize=0.5)

    y_Ha = 10.**(0.4*(np.log10(chi2_min)+np.log10(chi2.max())))
    col_Ha = tableau20[6]
    plt.text(4.4,y_Ha*1.1,r'H$\alpha$',color=col_Ha, ha='center')
    plt.plot([3.76,5.1], [y_Ha,y_Ha], c=col_Ha)
    plt.plot([3.76,3.76], [y_Ha*0.9,y_Ha*1.1], c=col_Ha)
    plt.plot([5.1,5.1], [y_Ha*0.9,y_Ha*1.1], c=col_Ha)
    #plot best [OIII] z
    z_Ha_best = zgrid[(zgrid > 3.76)&(zgrid < 5.1)][((chi2[(zgrid > 3.76)&(zgrid < 5.1)]).argmin())]
    chi2_Ha_best = chi2[(zgrid > 3.76)&(zgrid < 5.1)][((chi2[(zgrid > 3.76)&(zgrid < 5.1)]).argmin())]
    plt.scatter(z_Ha_best, chi2_Ha_best, c=col_Ha,marker='o',ec=tableau20[6], label=r'Ha best $\Delta\chi^2=$%0.2f' % (chi2_Ha_best - chi2_min))
    #plt.errorbar(6.15, 3.0, xerr=0.85,capsize=0.5)
            
    ax.set_xlabel('$z$')
    ax.set_ylabel(r'$ \chi (z)$')
    ax.grid()

    #ax.set_yticklabels([])

    ax.legend()
    return z_Ha_best, z_OIII_best

def plot_extra_temp(ix, zshow, color, alpha=0.5, ls='-'):
    z = zshow      #optionally set to show redshift
    fnu_i = self.fnu[ix, :]
    efnu_i = self.efnu[ix,:]
    zspec_i = self.ZSPEC[ix]

    lnp_i = self.lnp[ix,:]
    log_prior_i = self.full_logprior[ix,:].flatten()
    chi2_i = self.chi2_fit[ix,:]
    zgrid = self.zgrid
    ok_i = self.ok_data[ix,:]

    ndraws=100
    fitter='nnls'
    IGM_OBJECT = eazy.igm.Inoue14()

    chi2 = np.squeeze(chi2_i)
    prior = np.exp(log_prior_i)
    pz = np.exp(-(chi2-chi2.min())/2.) #*prior
    pz /= np.trapz(pz, zgrid)
    #pz = np.exp(lnp_i).flatten()

    ## SED        
    fnu_i = np.squeeze(fnu_i)*self.ext_redden*self.zp
    efnu_i = np.squeeze(efnu_i)*self.ext_redden*self.zp
    ok_band = (fnu_i/self.zp > self.param['NOT_OBS_THRESHOLD']) 
    ok_band &= (efnu_i/self.zp > 0)
    efnu_i[~ok_band] = self.param['NOT_OBS_THRESHOLD'] - 9.

    ## Evaluate coeffs at specified redshift
    tef_i = self.TEF(z)
    A = np.squeeze(self.tempfilt(z))
    chi2_i, coeffs_i, fmodel, draws = eazy.photoz.template_lsq(fnu_i, efnu_i, A, tef_i, self.zp, ndraws, fitter)

        
    ## Full SED
    templ = self.templates[0]
    tempflux = np.zeros((self.NTEMP, templ.wave.shape[0]),
                        dtype=self.ARRAY_DTYPE)

    for i in range(self.NTEMP):
        zargs = {'z':z, 'redshift_type':'nearest'}
        fnu = self.templates[i].flux_fnu(**zargs)*self.tempfilt.scale[i]
        try:
            tempflux[i, :] = fnu
        except:
            tempflux[i, :] = np.interp(templ.wave,
                                       self.templates[i].wave, fnu)
            
    templz = templ.wave*(1+z)


    igmz = templ.wave*0.+1
    lyman = templ.wave < 1300
    igmz[lyman] = IGM_OBJECT.full_IGM(z, templz[lyman])

    fnu_factor = 10**(-0.4*(self.param['PRIOR_ABZP']+48.6))
    flam_spec = 1.e29
    templf = np.dot(coeffs_i, tempflux)*igmz
    plt.plot(templz/1.e4, templf*fnu_factor*flam_spec, alpha=alpha, ls=ls,  zorder=-1, color=color, label='z={0:.2f}'.format(z))

def plot_cutouts(ix, obj, mosaics, gs, cutsize=100, boarder=30):
    x_min = np.max([0,              int(np.round(obj['X_IMAGE_det']))-cutsize])
    x_max = np.min([mosaics[0]['image'].shape[1], int(np.round(obj['X_IMAGE_det']))+cutsize])
    y_min = np.max([0,              int(np.round(obj['Y_IMAGE_det']))-cutsize])
    y_max = np.min([mosaics[0]['image'].shape[0], int(np.round(obj['Y_IMAGE_det']))+cutsize])

    posn_cut = [obj['X_IMAGE_det'] - x_min, obj['Y_IMAGE_det'] - y_min]
    zscale = ZScaleInterval()
    for ib, band in enumerate(mosaics):
        #make cutout
        data_cut =   band['image'][y_min:y_max,x_min:x_max]
        #now plot
        plt.subplot(gs[ib])
        if (np.all(np.isnan(data_cut))|np.all((data_cut == 0.0))): 
            z1 = 0.0
            z2 = 1.0
        else:
            z1, z2 = zscale.get_limits(data_cut)
        plt.imshow(data_cut, vmin=z1, vmax=z2, origin='lower', aspect='equal')
        plt.xlim(posn_cut[0] - boarder, posn_cut[0] + boarder)
        plt.ylim(posn_cut[1] - boarder, posn_cut[1] + boarder)



if __name__ == "__main__":
    """
    Main body of the routine to perform the analysis 
    """
    #set directories
    basedir = '/Users/bordoloi/Dropbox/Research/JWST/JWST_GTO/GTO/J0100/JWST/photoz/'
    mosaic_dir  = '/Users/bordoloi/Dropbox/Research/JWST/JWST_GTO/GTO/J0100/JWST/current_best/'
    save_dir    = basedir+'plot_photoz/'

    #Photometric catalog
    catalog_filename = basedir+"/J0100_photcat_v4_noisemodel_short.fits"
    basename = basedir+'J0100_photzcat_v4_EAZY_'            #basename for EAZY files, needs basename+'_output.h5'

    #workdir = basedir + 'photometry/EAZY/'
    #os.chdir(workdir)

    # Set what to save
    plot_images = False     #if false skip plotting cutouts
    save_plot   = False   

    #set v4 galaxy              #loop over objects here
    ix = 17807-1 #eyelash [OIII]
   


    #load image files so it's not done in each loop
    mosaics    = [dict(image=None, fname=mosaic_dir+'F606W_stack_F356Wref_final_drc_sci.fits',                  ext=0,      FILTER='F606W'),
                  dict(image=None, fname=mosaic_dir+'F775W_stack_F356Wref_final_drc_sci.fits',                  ext=0,      FILTER='F775W'),
                  dict(image=None, fname=mosaic_dir+'F850LP_stack_F356Wref_final_drc_sci.fits',                 ext=0,      FILTER='F850LP'),
                  dict(image=None, fname=mosaic_dir+'stack_F115W_pipe2_v2_pipeup_1.8.2_pub0988_20221028.fits',  ext='SCI',  FILTER='F115W'),
                  dict(image=None, fname=mosaic_dir+'stack_F200W_pipe2_v2.1_pipeup_1.8.2_pub0988_20221031.fits',ext='SCI',  FILTER='F200W'),
                  dict(image=None, fname=mosaic_dir+'stack_F356W_pipe2_v2.1_pipeup_1.8.2_pub0988_20221030.fits',ext='SCI',  FILTER='F356W')]

    if plot_images:
        for cur in mosaics:
            hdul = fits.open(cur['fname'])
            cur['image'] = hdul[cur['ext']].data

    #load data
    self = eazy.hdf5.initialize_from_hdf5(h5file=basename+'_output.h5')
    self.param['NOT_OBS_THRESHOLD'] = -90

    photcat = Table.read(catalog_filename) 


    plt.figure(figsize=(12,8))
    gs = GridSpec(2,6,hspace=0.25,wspace=0.2, left=0.09, bottom=0.06,right=0.97, top=0.95,height_ratios=[1,2])
    if plot_images: 
        plot_cutouts(ix, photcat[ix], mosaics, gs)
    ax = plt.subplot(gs[1,3:])
    z_Ha, z_OIII = chi2_dist(self, ix, ax)
    ax = plt.subplot(gs[1,:3])
    show_fits(self, ix, ax, zshow=None, show_components=False, show_redshift_draws=False, template_color='grey')
    plot_extra_temp(ix, z_Ha, tableau20[6], alpha=0.5, ls='-')
    plot_extra_temp(ix, z_OIII, tableau20[0], alpha=0.5, ls='-')
    #plt.ylim(-2E-3,10E-3)
    plt.show()
    if save_plot: 
        plt.savefig(save_dir+'id'+str(ix)+'_v4_chi2_dist.png')
