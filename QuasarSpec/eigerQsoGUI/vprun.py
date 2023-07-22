#!/usr/bin/env python3

import sys

if len(sys.argv) < 3:
    print("Error: more arguments needed")
    sys.exit()

import eiger.QuasarSpec.loadQsoSpec as spec
import matplotlib.pyplot as plt
from numpy import sqrt, array
from pypeit.core.wave import airtovac
import astropy.units as u
from mcvp import model as vm
from mcvp import vfit as vf
from astropy.convolution import convolve, Gaussian1DKernel
import pickle
import corner

# FIRE: 0.6" slit, 0.15" pixels yields 4 pixels per resolution element FWHM
fire_kernel    = Gaussian1DKernel(stddev=4/2.355)

## XShooter spectral resolution:
## https://www.eso.org/sci/facilities/paranal/instruments/xshooter/doc/VLT-MAN-ESO-14650-4942_P93.pdf
## Table 13
## Headers indicate that the 0.9" slit was used (DECKER keyword)
## For this combination, the slit is 4.2 pixels in NIR, 6.0 in VIS (possibly 4.8??)
xsh_nir_kernel = Gaussian1DKernel(stddev=2.2/2.355)
#xsh_vis_kernel = Gaussian1DKernel(stddev=6.0/2.355)
xsh_vis_kernel = Gaussian1DKernel(stddev=4.8/2.355)

# HIRES data mostly take at 6 km/s resolution, R=50,000, or 0.86" slit, or 3 pixels
# See https://www2.keck.hawaii.edu/inst/hires/manual2.pdf
hires_kernel   = Gaussian1DKernel(stddev=3.0/2.355)

# MOSFIRE (0.7" Slit) per webpage:
# and our slits were 0.7" on the mask
#
# Y = 3380
# J = 3310
# H = 3660
# K = 3620

# VP fit model class
m = vm.Model()

##############

print ("Loading J0100 spectra....")

obj_index = 4
spec_current = spec.loadQsoSpec(obj_index,revision='current')

m.addspec('FIRE',spec_current['FIRE']['wave'],spec_current['FIRE']['flux']/spec_current['FIRE']['cont'],1/sqrt(spec_current['FIRE']['ivar'])/spec_current['FIRE']['cont'],fire_kernel)

m.addspec('XSH_NIR',spec_current['XSH_NIR']['wave'],spec_current['XSH_NIR']['flux']/spec_current['XSH_NIR']['cont'],1/sqrt(spec_current['XSH_NIR']['ivar'])/spec_current['XSH_NIR']['cont'],xsh_nir_kernel)

m.addspec('XSH_VIS',spec_current['XSH_VIS']['wave'],spec_current['XSH_VIS']['flux']/spec_current['XSH_VIS']['cont'],1/sqrt(spec_current['XSH_VIS']['ivar'])/spec_current['XSH_VIS']['cont'],xsh_vis_kernel)

if (obj_index == 1):
    m.addspec('HIRES',array(airtovac(spec_current['HIRES']['wave']*u.AA)),spec_current['HIRES']['flux']/spec_current['HIRES']['cont'],1/sqrt(spec_current['HIRES']['ivar'])/spec_current['HIRES']['cont'],hires_kernel)
else:
    m.addspec('HIRES',spec_current['HIRES']['wave'],spec_current['HIRES']['flux']/spec_current['HIRES']['cont'],1/sqrt(spec_current['HIRES']['ivar'])/spec_current['HIRES']['cont'],hires_kernel)

################

with open(sys.argv[1]) as f:
    # readlines returns an array, read is one long exec, need both
    exec(f.read())

with open(sys.argv[1]) as f:
    # readlines returns an array, read is one long exec, need both
    model_lines = f.readlines()
    
################

# Number of fit parameters = 2 * Ncomponents + Nions
# Each component has a redshift and b
# Each ion has a column density

ncomponents = sum(['addcomponent' in tmp for tmp in model_lines])
nions       = sum(['addion' in tmp for tmp in model_lines])
print(f"Ncomponents, Nions = {ncomponents}, {nions}")
nparams = 2*ncomponents + nions

print ("Running the sampler...")
sampler = vf.runmc(m,nwalkers=2*nparams,nruns=3000)
print("Reformatting the chain")
samples = sampler.chain[:,-1200:,:].reshape((-1,nparams))

print ("Finished, writing results to disk...")

fname = sys.argv[2]
outputdict = {'model':m, 'samples':samples}
with open(fname, "wb") as fp:
    pickle.dump(outputdict, fp, pickle.HIGHEST_PROTOCOL)

fig = corner.corner(samples)
fig.savefig('vpfit_corner.pdf')

