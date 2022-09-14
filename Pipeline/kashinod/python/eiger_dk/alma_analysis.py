from astropy.io import fits
import astropy.units as u
from astropy import constants as aconst
import numpy as np
# Some functions

def get_axis(header, axis=1):
    crval = header.get('CRVAL%1d'%(axis))
    crpix = header.get('CRPIX%1d'%(axis))
    cdelt = header.get('CDELT%1d'%(axis))
    num = header.get('NAXIS%1d'%(axis))
    return crval + (np.arange(num) - crpix + 1) * cdelt

def mygauss0(x, amp):
    return amp * np.exp( - (x - 0.)**2 / (2 * 1.**2))

def mygauss(x, amp, mu, sigma):
    return amp * np.exp( - (x - mu)**2 / (2 * sigma**2))

def mygauss_off(x, amp, mu, sigma, f0):
    return amp * np.exp( - (x - mu)**2 / (2 * sigma**2)) + f0

def mygauss_dbl(x, amp, mu, sigma, amp2, mu2, sigma2):
    return amp * np.exp( - (x - mu)**2 / (2 * sigma**2)) + amp2 * np.exp( - (x - mu2)**2 / (2 * sigma2**2))


def mygauss_lin(x, amp, mu, sigma, c0, c1):
    return mygauss(x, amp, mu, sigma) + c0 + c1*(x-mu)

    
def rms(x):
    N=x.size
    return np.sqrt(np.sum(x**2)/N)

def nanrms(x):
    N=(np.where(np.isfinite(x))[0]).size
    return np.sqrt(np.nansum(x**2)/N)

def ftov(f,f0):
    return (f0-f)/f0 * (aconst.c).to(u.km/u.s)

def correlated_noise(rms, theta_n, A, theta_maj, theta_min, phi, I):

    # Params A, sig_maj, sig_min, cen_x, cen_y, pa
    #                     A, theta_maj, theta_min,    x0,    y0,   phi
    a_maj = np.array([3./2.,     5./2.,     1./2., 5./2., 1./2., 1./2.])
    a_min = np.array([3./2.,     1./2.,     5./2., 1./2., 5./2., 5./2.])

    
    rho = A/rms * np.sqrt(theta_maj * theta_min)/(2*theta_n) * (1. + (theta_n/theta_maj)**2.)**(a_maj/2.) * (1. + (theta_n/theta_min)**2.)**(a_min/2.)
    
    sig_A, sig_Maj, sig_Min = np.sqrt(2.) * np.array([A, theta_maj, theta_min]) / rho[0:3]
    sig_x0 = np.sqrt(2) * theta_maj / (rho[3] * np.sqrt(8*np.log(2)))
    sig_y0 = np.sqrt(2) * theta_min / (rho[4] * np.sqrt(8*np.log(2)))
    sig_phi = np.degrees(2. * (theta_maj * theta_min)/(theta_maj**2. - theta_min**2) /rho[5])
    sig_I = I * np.sqrt( (sig_A/A)**2. + (theta_n**2/(theta_maj*theta_min)) * ((sig_Maj/theta_maj)**2. + (sig_Min/theta_min)**2.))
    sig_ra = np.sqrt( sig_x0**2 * np.sin(np.radians(phi))**2. + sig_y0**2. * np.cos(np.radians(phi))**2.)
    sig_dec = np.sqrt( sig_x0**2 * np.cos(np.radians(phi))**2. + sig_y0**2. * np.sin(np.radians(phi))**2.)
    
    return sig_A, sig_Maj, sig_Min, sig_ra, sig_dec, sig_phi, sig_I

def get_outline(mfs_field_mask):
    mfs_outline_x = []
    mfs_outline_y = []

    s=0
    for ix in np.arange(mfs_field_mask.shape[1]):
        tmp = np.where(np.isnan(mfs_field_mask[:,ix]))[0]
        if tmp.size>0:
            if s==0:
                s=1
                if tmp.size>2:
                    for j in tmp[0:-1]:
                        #print(ix,j)
                        mfs_outline_x.append(ix)
                        mfs_outline_y.append(j)
                elif tmp.size==2:
                    #print(ix,tmp[0])
                    mfs_outline_x.append(ix)
                    mfs_outline_y.append(tmp[0])
                elif tmp.size==1:
                    pass
                    #print('skip')
                else:
                    raise ValueError()
            else:                
                #print(ix,tmp[0])
                mfs_outline_x.append(ix)
                mfs_outline_y.append(tmp[0])
    s=0
    for ix in np.flip(np.arange(mfs_field_mask.shape[1])):
        tmp = np.flip(np.where(np.isnan(mfs_field_mask[:,ix]))[0])
        if tmp.size>0:
            if s==0:
                s=1
                if tmp.size>2:
                    for j in tmp[0:-1]:
                        #print(ix,j)
                        mfs_outline_x.append(ix)
                        mfs_outline_y.append(j)
                elif tmp.size==2:
                    #print(ix,tmp[0])
                    mfs_outline_x.append(ix)
                    mfs_outline_y.append(tmp[0])
                elif tmp.size==1:
                    pass
                    #print('skip')
                else:
                    raise ValueError()
            else:                
                #print(ix,tmp[0])
                mfs_outline_x.append(ix)
                mfs_outline_y.append(tmp[0])
    return np.array(mfs_outline_x), np.array(mfs_outline_y)


#def derive_error_spectrum(cubdata,
                          
def convert_L_solar_to_Kkmspc2(L, nu_rest):
    
    L = L.to_value(u.solLum)/(3.2e-11 * nu_rest.to_value(u.GHz)**3.)
    return L * u.K * u.km * u.pc**2 / u.s

def convert_L_Kkmspc2_to_solar(L, nu_rest):
    
    L_new = L.to_value(u.K*u.km/u.s*u.pc**2)*(3.2e-11 * nu_rest.to_value(u.GHz)**3.)
    return L_new * u.solLum

def convert_F_Jykms_to_cgs(F, nu_obs):

    
    F_new = F.to(u.Jy * u.km / u.s) / aconst.c.to(u.km/u.s) * nu_obs
    return F_new.to(u.erg * u.cm**(-2) * u.s**(-1))

def convert_F_cgs_to_Jykms(F, nu_obs):

    
    F_new = F.to(u.Jy * u.Hz) * aconst.c.to(u.km/u.s) / nu_obs
    return F_new.to(u.Jy * u.km / u.s)
