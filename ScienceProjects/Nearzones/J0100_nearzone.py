from   mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from   astropy.table import Table
from   astropy.cosmology import Planck18
from   astropy.coordinates import SkyCoord
import astropy.units as u
import astropy.constants as c
import numpy as np
import eiger.QuasarSpec.loadQsoSpec as spec

class Nearzone:

    def __init__(self):

        self.galtable = Table.read('../../Database/J0100_photcat_v2_O3emitters_Systems_28092022.fits')
        self.zqso     = 6.3258
        self.qso_ra   = 15.054268
        self.qso_dec  = 28.040497
        self.metal_absorbers = [6.1872,6.143476,6.143089,6.1438,6.1115,6.0116,5.7977]

        spec_current = spec.loadQsoSpec(4,revision='current',
                                        offline=True,
                                        spectrographs=['XShooter'])['XSH_VIS']


        spec_z = spec_current['wave'] / 1216.6701 - 1
        keep = np.where(np.logical_and(spec_z > 5.95,spec_z < self.zqso))

        self.qso_flux = (spec_current['flux']/12.0)[keep]
        self.qso_wave = (spec_current['wave'])[keep]
        

    ###################################################
        
    def plot3d(self):

        xs = [0]
        ys = [0]
        zs = [0]

        # Find the galaxy locations in a coordinate system where the QSO is at the origin
        for line in self.galtable:
            ra =  line['RA_MANUAL']
            dec = line['DEC_MANUAL']
            z = float(line['z_O3doublet_combined'])
            if (z<5.95 or z > self.zqso):
                continue

            distance_vector = self.physical_distance(self.zqso, z, ra, dec)
            
            ys.append(distance_vector['dR'])
            xs.append(distance_vector['dra'])
            zs.append(distance_vector['ddec'])

        # Catalog the metal absorber locations
        xs_2 = []
        ys_2 = []
        zs_2 = []

        for z in self.metal_absorbers:
            distance_vector = self.physical_distance(self.zqso, z, self.qso_ra, self.qso_dec)
            xs_2.append(0)
            ys_2.append(distance_vector['dR'])
            zs_2.append(0)
            
        # creating figure
        fig = plt.figure()
        ax = Axes3D(fig,auto_add_to_figure=False)
        fig.add_axes(ax)
        
        # Plot galaxy points
        ax.scatter(xs, ys, zs, color='green')
        # Plot a continuum
        ax.plot([0,0,0],[0,0,-25],color='red',linestyle='dashed',alpha=0.5)
        # Plot the QSO absorbers
        ax.scatter(xs_2, ys_2, zs_2, color='red')
        
        # setting title and labels
        ax.set_xlim(-2,2)
        ax.set_ylim(-22,0)
        ax.set_zlim(-2,2)
        ax.set_xlabel("dx (pMpc)")
        ax.set_zlabel("dy (pMpc)")
        ax.set_ylabel("Radius (pMpc)")

        spec_z = self.qso_wave / 1216.6701 - 1
        spec_R = np.array([self.physical_distance(self.zqso,this_z,self.qso_ra,self.qso_dec)['dR'] for this_z in spec_z])

        zz = spec_R
        yy = 2.0*self.qso_flux-2.0
        xx = np.zeros(len(zz))

        ax.plot3D(xx,zz,yy,color='r',alpha=0.6)
        
        # This sets the aspect ratio for a non-uniform box (by using the data values)
        ax.set_box_aspect([ub - lb for lb, ub in (getattr(ax, f'get_{a}lim')() for a in 'xyz')])
        plt.show()

        
    ###################################################

    def Gamma_galaxies(self):

        spec_z = self.qso_wave / 1216.6701 - 1
        Gamma_galaxies = np.zeros(len(spec_z))

        for gal in self.galtable:
            ra =  gal['RA_MANUAL']
            dec = gal['DEC_MANUAL']
            zgal = float(gal['z_O3doublet_combined'])
            F115W = (gal['fnu_F115W_AUTO_apcor'] * u.nJy).to('erg/(s*cm**2*Hz)')

            # I realize that F115W here is in erg/cm2/s/Hz and D_L is in Mpc
            # The Mpc will get divided out in ~5 lines of code, so will leave it
            # as-is to avoid floating point overflow/roundoff errors
            
            L_1600 = F115W * 4 * np.pi * (Planck18.luminosity_distance(zgal))**2

            distance = self.physical_distance(spec_z, zgal, ra, dec)['resultant'] * u.Mpc

            # ASSUMING 100% ESCAPE FRACTION - RETURN TO THIS
            fnu_here = (L_1600 / (4 * np.pi * distance**2)).value
            Gamma_galaxies += fnu_here

        fig, ax = plt.subplots(2)
        ax[0].plot(spec_z, Gamma_galaxies)
        ax[0].set_ylim(0,5e-19)
        
        # ax[0].set_yscale('log')

        ax[1].plot(spec_z, self.qso_flux)
        ax[1].plot(spec_z, Gamma_galaxies/5e-19,alpha=0.5)
        ax[1].set_ylim(0,1.1)

        plt.show()
    
    ###################################################
    
    def physical_distance(self,zreference, zgal, ra, dec):

        # This could be improved by calculating the full line element,
        # this is a linear approximation, assuming pythagorean theorem

        # This weird line is so that we can pass vectors or scalars
        if (hasattr(zreference, "__len__") == False):
            # cscale converts decimal degree offsets into physical Mpc
            cscale = 1./Planck18.arcsec_per_kpc_proper(zgal).value / 1000.0
            dra  = (ra-self.qso_ra)*np.cos(dec*np.pi/180.0) * cscale * 3600
            ddec = (dec-self.qso_dec) * cscale * 3600
            z1 = np.max([zreference,zgal])
            z2 = np.min([zreference,zgal])
            dR   = -1 * ((1+z1)/(1+z2)-1) * c.c.to("km/s") / Planck18.H(zreference)
            total_distance = np.sqrt(dra**2+ddec**2+dR.value**2)
            result = {'dra':dra, 'ddec':ddec,'dR':dR.value, 'resultant':total_distance}

        else:
            # cscale converts decimal degree offsets into physical Mpc
            cscale = 1./Planck18.arcsec_per_kpc_proper(zgal).value / 1000.0
            dra  = (ra-self.qso_ra)*np.cos(dec*np.pi/180.0) * cscale * 3600
            ddec = (dec-self.qso_dec) * cscale * 3600

            ratio = (1+np.array(zreference))/(1+zgal)

            dR   = -1 * (ratio-1) * c.c.to("km/s") / Planck18.H(zreference)
            dR[ratio < 1] = -1 * (1/ratio[ratio<1]-1) * c.c.to("km/s") / Planck18.H(zgal)
                
            total_distance = np.sqrt(dra**2+ddec**2+dR.value**2)
            result = {'dra':dra, 'ddec':ddec,'dR':dR.value, 'resultant':total_distance}

        return(result)
        
        
    #############################
    

