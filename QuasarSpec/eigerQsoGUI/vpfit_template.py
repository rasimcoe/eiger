#!/usr/bin/env python3
#
# Voigt profile fitting template for a locally-registered object.
#
# Usage:
#   python vpfit_template.py  <model_file.in>  <output.pickle>
#
# Before running:
#   1. Register the object in $EIGER_CACHE/local_objects.yaml
#   2. Set OBJ_NAME and INSTRUMENTS below to match your object
#   3. Adjust the kernel FWHM values if you are using non-standard slit widths
#   4. Write a <model_file.in> that calls m.addcomponent / m.addion /
#      m.addtransition (see any of the existing *.in files for syntax)

import sys

if len(sys.argv) < 3:
    print("Usage: vpfit_template.py <model_file.in> <output.pickle>")
    sys.exit(1)

import eiger.QuasarSpec.loadQsoSpec as spec
from mcvp import model as vm
from mcvp import vfit as vf
from astropy.convolution import Gaussian1DKernel
from numpy import sqrt
import pickle
import corner

# ── Object configuration ────────────────────────────────────────────────────

OBJ_NAME = "J0000+0000"   # must match the name in local_objects.yaml

# List only the instruments that this object actually has spectra for.
# Keys must match what loadLocalSpec() returns (FIRE, XSH_VIS, XSH_NIR,
# HIRES, MOSFIRE_Y, MOSFIRE_J, MOSFIRE_H, MOSFIRE_K, FIRE_XSH).
INSTRUMENTS = ['FIRE', 'XSH_VIS', 'XSH_NIR']

# ── Instrumental resolution kernels (stddev in pixels) ──────────────────────
# Adjust FWHM and pixel scale for your specific observations.
#
# FIRE:        0.6" slit, 0.15"/pix  → ~4 px FWHM  (R ≈ 6000)
# XShooter NIR: 0.9" slit            → ~2.2 px FWHM
# XShooter VIS: 0.9" slit            → ~4.8 px FWHM
# HIRES:       6 km/s resolution     → ~3 px FWHM   (R ≈ 50000)

kernels = {
    'FIRE':      Gaussian1DKernel(stddev=4.0   / 2.355),
    'XSH_NIR':  Gaussian1DKernel(stddev=2.2   / 2.355),
    'XSH_VIS':  Gaussian1DKernel(stddev=4.8   / 2.355),
    'HIRES':     Gaussian1DKernel(stddev=3.0   / 2.355),
    'FIRE_XSH':  Gaussian1DKernel(stddev=4.0   / 2.355),
    'MOSFIRE_Y': Gaussian1DKernel(stddev=3.5   / 2.355),
    'MOSFIRE_J': Gaussian1DKernel(stddev=3.5   / 2.355),
    'MOSFIRE_H': Gaussian1DKernel(stddev=3.5   / 2.355),
    'MOSFIRE_K': Gaussian1DKernel(stddev=3.5   / 2.355),
}

# ── Load spectra ─────────────────────────────────────────────────────────────

print(f"Loading {OBJ_NAME} spectra...")
sp = spec.loadLocalSpec(OBJ_NAME)

if sp is None:
    print("ERROR: could not load spectra. Check local_objects.yaml and EIGER_CACHE.")
    sys.exit(1)

# ── Build the VP model and add spectra ───────────────────────────────────────

m = vm.Model()

for inst in INSTRUMENTS:
    if inst not in sp:
        print(f"WARNING: {inst} not found in loaded spectra — skipping")
        continue
    s = sp[inst]
    norm_flux = s['flux'] / s['cont']
    norm_err  = 1.0 / sqrt(s['ivar']) / s['cont']
    m.addspec(inst, s['wave'], norm_flux, norm_err, kernels[inst])
    print(f"  Added {inst}")

# ── Load model definition from .in file ──────────────────────────────────────

with open(sys.argv[1]) as f:
    exec(f.read())

with open(sys.argv[1]) as f:
    model_lines = f.readlines()

# ── Count parameters and run MCMC ────────────────────────────────────────────
# nparams = 2 * Ncomponents (z + b per component) + Nions (N per ion)

ncomponents = sum('addcomponent' in line for line in model_lines)
nions       = sum('addion'       in line for line in model_lines)
nparams     = 2 * ncomponents + nions
print(f"Ncomponents={ncomponents}, Nions={nions}, Nparams={nparams}")

print("Running the sampler...")
sampler = vf.runmc(m, nwalkers=2 * nparams, nruns=3000)

print("Reformatting the chain...")
samples = sampler.chain[:, -1200:, :].reshape((-1, nparams))

print("Writing results to disk...")
fname = sys.argv[2]
with open(fname, "wb") as fp:
    pickle.dump({'model': m, 'samples': samples}, fp, pickle.HIGHEST_PROTOCOL)

fig = corner.corner(samples)
fig.savefig('vpfit_corner.pdf')
print("Done.")
