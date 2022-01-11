import os
import jwreftools


dir = '/str1/kashinod/JWST/mirage2/referenceFiles/mirage_data/nircam/GRISM_NIRCAM/current'
conf_fils = ['NIRCAM_F356W_modA_R.conf', 'NIRCAM_F356W_modB_R.conf']

for conf in conf_fils:
    conffile=os.path.join(dir,conf)
    outname = conf.split(".")[0]+'.asdf'
    print('conffile: ', conffile)
    print('outname: ', outname)

    pupil = "GRISM" + conffile.split(".")[0][-1]
    module = conffile.split(".")[0][-3]
    print('pupil: ', pupil)
    print('module: ', module)
    jwreftools.nircam.create_grism_specwcs(conffile=conffile, outname=outname)
