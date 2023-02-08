This directory contains script to:

1) Optimally extract 1D EMLINE / SCI spectra of the "total object" (i.e. clumps are merged). The spatial profile follows that of [OIII]5008 or Halpha. 
extract_total_1D_O3_emitters_SCI.py
extract_total_1D_O3_emitters.py
extract_total_Halpha.py
extract_stars.py

2) measure redshifts and line-fluxes  from the 1D spectrum
fit_redshift_1D_O3doublet.py
fit_flux_1D_O3shape.py

3) measure EWs based on line + photometry (currently applicable for [OIII] emitters only, but can be rewritten)
measure_EWs.py
measure_EWs_withF200W.py -- this is outdated

4) stack spectra of [OIII] emitters in 1D or 2D
stack_2D.py
stack_1D.py


