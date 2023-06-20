#EIGER EAZY photoz
These scripts run EAZY on the mutliband photometry catalog.

Requirements:
eazy-py


Codes:
EAZY_J1148_v3_noise_model.py
Is an example script to run EAZY on EIGER catalogs. It's specific to the field because you need to change the entries depending on the HST bands. This example is for J1148 which has F606W, F775W and F850LP.
You also need to change the Galactic EBV for different fields. This value is from Planck.

photz_picket_fence_line_prior.py 
Merges the photometry and photoz information into one catalog. The line prior stuff is currently not used.

plot__photoz_fit.py
Shows an example of how to plot the photoz distribution of a source.

