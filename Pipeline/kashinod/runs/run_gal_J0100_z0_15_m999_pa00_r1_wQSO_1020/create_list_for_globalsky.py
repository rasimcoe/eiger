import numpy as np
import glob
import os

directory='calibrated_img2'

filters=['f115w', 'f200w', 'f356w_imaging', ]#'f356w_wfss']
modules = ['a','b']
detectors_dict = {'f115w':[1,2,3,4],
                  'f200w':[1,2,3,4],
                  'f356w_wfss': [5],
                  'f356w_imaging': [5],}

pointings_f115w = np.ravel([np.arange( 1,13), np.arange(29,41), np.arange(57,69), np.arange(85, 97)])
pointings_f200w = np.ravel([np.arange(14,26), np.arange(42,54), np.arange(70,82), np.arange(98, 110)])
pointings_f356w_wfss = np.ravel([np.arange( 1,13), np.arange(14,26), 
                                 np.arange(29,41), np.arange(42,54), 
                                 np.arange(57,69), np.arange(70,82), 
                                 np.arange(85, 97), np.arange(98, 110)])
pointings_f356w_imaging = np.array([26,27,28,54,55,56,82,83,84,110,111,112])
pointings_dict = {'f115w': pointings_f115w, 
                  'f200w': pointings_f200w, 
                  'f356w_wfss': pointings_f356w_wfss, 
                  'f356w_imaging': pointings_f356w_imaging}

for filter in filters:
    pointings = pointings_dict[filter]
    detectors = detectors_dict[filter]
    print(filter)
    print(pointings)
    for module in modules:
        for detector in detectors:
            list_fil = 'list_nrc'+module+str(detector)+'_'+filter+'.txt'
            print(list_fil)
            ls=[]
            for pointing in pointings:
                fil_search=os.path.join(directory, 'jw0124300100?_01101_'+str(pointing).zfill(5)+
                                           '_nrc'+module+str(detector)+'_cal.fits')
                fil=glob.glob(fil_search)
                #print(fil_search)
                ls.append(fil[0])
            #print(ls)
            np.savetxt(list_fil, ls, fmt="%s")
