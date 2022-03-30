from eiger.QuasarSpec import loadQsoSpec
from astropy.convolution import convolve, Gaussian1DKernel
from matplotlib import pyplot as plt
import numpy as np
from scipy.signal import find_peaks, peak_widths
from astropy.table import Table
from os import getenv

def findAbsLines(spec, instrument):

    if (instrument == 'FIRE'):
        # FWHM=4 pixels (50 km/s, 12.5 km/s/pixel, 2.355 converts to sigma)
        lsf_kernel = Gaussian1DKernel(stddev=4/2.355) 

    elif (instrument == 'HIRES'):
        # 6km/s, 1.4 km/s/pixel 
        lsf_kernel= Gaussian1DKernel(stddev=4.285/2.355)
        spec[instrument]['cont'] = np.ones(len(spec[instrument]['flux']))

    elif (instrument == 'XSH_NIR'):
        # XShooter manual Table 13, 0.6" slit 2.9 pixels/FWHM, R=7770
        lsf_kernel= Gaussian1DKernel(stddev=2.9/2.355) 

    elif (instrument == 'XSH_VIS'):
        # XShooter manual Table 13, 0.7" slit 4.8 pixels/FWHM, R=11000
        lsf_kernel= Gaussian1DKernel(stddev=4.8/2.355) 

    elif (instrument == 'MOSFIRE'):
        # ? km/s, ? km/s/pixel
        print("MOSFIRE not yet supported")
        return(None)
        lsf_kernel= Gaussian1DKernel(stddev=2.9/2.355) 

        
    signal     = (1.0 - spec[instrument]['flux']/spec[instrument]['cont'])
    ivar_norm  = (spec[instrument]['ivar'] * spec[instrument]['cont']**2)

    numer = convolve(signal, lsf_kernel, boundary='extend')
    denom = 1/np.sqrt(convolve(ivar_norm, lsf_kernel, boundary='extend'))

    pks,props = find_peaks(numer/denom,height=5)

    widths, width_height, left_ips, right_ips = peak_widths(numer/denom, pks)

    ###########################

    wvmin = (1+spec['z_em']) * 1216.0
    gd    = np.where(np.logical_and((spec[instrument]['wave'][pks] > wvmin), \
                                    (right_ips > 0)))
    
    pkind         = pks[gd]
    pkwv          = spec[instrument]['wave'][pks[gd]]
    pkfx          = spec[instrument]['flux'][pks[gd]]/spec[instrument]['cont'][pks[gd]]
    pksnr         = numer[pks[gd]]/denom[pks[gd]]
    pk_leftwv     = spec[instrument]['wave'][np.floor(left_ips[gd]).astype(np.int)]
    pk_rtwv       = spec[instrument]['wave'][np.floor(right_ips[gd]).astype(np.int)]

    abslines = Table([pkwv,pksnr,pkind,pk_leftwv,pk_rtwv], masked=True, \
                     names=('PeakWave','PeakSNR','Index','LeftBound','RightBound'))

    return(abslines)

#################################################

def idAbsLines(linelist, spec):

    sorted_list = linelist[linelist.argsort(keys='PeakSNR',reverse=True)]
    nones = np.full(len(sorted_list), None)
    tt = np.zeros(len(sorted_list))
    sorted_list.add_columns([nones,tt,tt],names=('ion','restwv','redshift'))

    linelist = Table.read('atomic_data.txt', format='ascii')

    dv_max = 35.0
#    dv_max = 6.0
    
    # Scan for MgII Doublets
    for line1 in sorted_list:
        testwv = line1['PeakWave'] * 2803.5314853/2796.3542699
        for line2 in sorted_list:
            if (dv(line2['PeakWave'],testwv) < dv_max):
                z = line1['PeakWave'] / 2796.3542699 - 1.0
                if (z > spec['z_em']+0.1):
                    continue
                print(f"MgII Found: z={z}")
                line1['ion'] = 'MgII'
                line2['ion'] = 'MgII'
                line1['redshift'] = z
                line2['redshift'] = z
                line1['restwv'] = 2796.3542699
                line2['restwv'] = 2803.5314853
                break

    # Scan for CIV Doublets
    for line1 in sorted_list:
        testwv = line1['PeakWave'] * 1550.77845 / 1548.2049
        for line2 in sorted_list:
            if (dv(line2['PeakWave'],testwv) < dv_max):
                z = line1['PeakWave'] / 1548.2049 - 1.0
                if (z > spec['z_em']+0.1):
                    continue
                print(f"CIV Found: z={z}")
                line1['ion'] = 'CIV'
                line2['ion'] = 'CIV'
                line1['redshift'] = z
                line2['redshift'] = z
                line1['restwv'] = 1548.2049 
                line2['restwv'] = 1550.77845
                break
            
    # Check FeII against remaining stragglers
    for line1 in sorted_list:
        if (line1['ion'] == None):
            testwv = line1['PeakWave'] * 2600.1724835 / 2586.6495659
            for line2 in sorted_list:
                if (dv(line2['PeakWave'],testwv) < dv_max):
                    z = line1['PeakWave'] / 2586.6495659 - 1.0
                    if (z > spec['z_em']+0.1):
                        continue
                    print(f"FeII Found: z={z}")
                    line1['ion'] = 'FeII'
                    line2['ion'] = 'FeII'
                    line1['redshift'] = z
                    line2['redshift'] = z
                    line1['restwv'] = 2586.6495659
                    line2['restwv'] = 2600.1724835
                    break

    for line1 in sorted_list:
        if (line1['ion'] == None):
            testwv = line1['PeakWave'] * 2382.7641781 / 2344.2129601
            for line2 in sorted_list:
                if (dv(line2['PeakWave'],testwv) < dv_max):
                    z = line1['PeakWave'] / 2344.2129601 - 1.0
                    if (z > spec['z_em']+0.1):
                        continue
                    print(f"FeII Found: z={z}")
                    line1['ion'] = 'FeII'
                    line2['ion'] = 'FeII'
                    line1['redshift'] = z
                    line2['redshift'] = z
                    line1['restwv'] = 2344.2129601
                    line2['restwv'] = 2382.7641781
                    break
                
    redshifts = np.unique(sorted_list['redshift'])

    # Find other lines associated with previously identified doublets
    for z in redshifts:
        for linecand in linelist:
            obswave = linecand['wave'] * (1 + z)
            for foundline in sorted_list:
                if (dv(foundline['PeakWave'],obswave) < dv_max and \
                    foundline['ion'] == None):
                    print(f" Matched {linecand['ion']}")
                    foundline['ion'] = linecand['ion']
                    foundline['restwv'] = linecand['wave']
                    foundline['redshift'] = z

    return(sorted_list)


#################################################

def plotAbsLines(spec, linelist, instrument='FIRE'):

    plt.step(spec[instrument]['wave'], \
             spec[instrument]['flux']/spec[instrument]['cont'])
    plt.plot([8000,30000], [0,0])
    
    for line in linelist:

        ff = spec[instrument]['flux'][line['Index']] / spec[instrument]['cont'][line['Index']]
        plt.plot([line['PeakWave'],line['PeakWave']],[ff,1],color='r',linestyle='dotted')
        plt.plot([line['LeftBound'],line['RightBound']],[1,1],color='r',linestyle='-')

        if (line['restwv'] != 0):
            plt.text(line['PeakWave'],1.1,"{} {:4.0f} (z={:5.3f})".format(line['ion'],np.trunc(line['restwv']),line['redshift']), rotation='vertical')

        
    plt.ylim(-20,60)
    plt.ylim(-0.25,2.0)
    plt.xlim(9000,9150)
    plt.show()
    
##################################################

def dv(wave1, wave2):
    c  = 299792.45800 # km/s
    dv = (wave1-wave2)/((wave1+wave2)/2.0) * c
    return(np.abs(dv))
