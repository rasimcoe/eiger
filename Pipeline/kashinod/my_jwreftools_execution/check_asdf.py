
import asdf
import os

crds_dir = '../crds_cache/references/jwst/nircam'
spwcs_modA = 'jwst_nircam_specwcs_0010.asdf'
spwcs_modB = 'jwst_nircam_specwcs_0009.asdf'


asdf_list = [os.path.join(crds_dir, spwcs_modA), 
             os.path.join(crds_dir, spwcs_modB)]
asdf_list.append('NIRCAM_F356W_modA_R.asdf')
asdf_list.append('NIRCAM_F356W_modB_R.asdf')

for asdf_fil in asdf_list:
    print('---------------------------------')
    print('Asdf name: ', asdf_fil)
    tmp = asdf.open(asdf_fil)
    print('meta: ', tmp['meta'])
    print('displ: ', tmp['displ'])
    print('dispx: ', tmp['dispx'])
    print('dispy: ', tmp['dispy'])
    print('invdispl; ', tmp['invdispl'])
    print('invdispx; ', tmp['invdispx'])
    print('invdispy; ', tmp['invdispy'])
    try:
        print('orders: ', tmp['orders'])
    except:
        print('order: ', tmp['order'], ' #order instead of orders')
