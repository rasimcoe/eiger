import numpy
import numpy
import matplotlib.pyplot as plt
from jwst import datamodels
from scipy.stats import binned_statistic_2d
import os
from astropy.visualization import simple_norm
import h5py
import numpy as np
from astropy.io import fits
import grismconf
import time
from scipy.interpolate import interp1d
from multiprocessing import Pool

import eiger_tracing_dk as eiger_tracing


### DK added
from astropy.wcs import WCS, utils
import astropy.units as u




field='J0100'
grismconf_calib_version=4 
ysize=51



directimage='/scratch/kashinod/EIGER/J0100/current_best/stack_F356W_pipe4_v3_fluxcal_pipeupdate_20220920.fits'

### Output directory

FOLDER='/scratch/EIGER/'#BROAD/SPECTRA_COLSEL/'


#CATALOG WITH SOURCES TO EXTRACT // CAN ALSO SKIP AND HAVE IDlist,RAlist and DEClist manually
CATALOG='/scratch/EIGER/identification/J0100_photcat_v4_within20arcsec_QSO.fits'
CATALOG='/scratch/EIGER/BROAD/J0100_photcat4_BROADsel_16022023.fits'
cat=fits.open(CATALOG)

data=cat[1].data
IDlist=data.field('NUMBER')
RAlist=data.field('ALPHA_J2000_det')
DEClist=data.field('DELTA_J2000_det')
#Xlist=data.field('X_IMAGE_det')

#RAlist=data.field('RA_MANUAL')
#DEClist=data.field('DEC_MANUAL')


#This bit is a placeholder for "bookkeeping" (i.e. keeping track of lines that are seen in this spectrum, but also seen in other spectra). Implemented for simulated data, not for real data.
"""
N_in_others=data.field('Nlines_appear_in_others')
wav_duplicates=data.field('wav_duplicates')
IDs_others=data.field('IDs_others')

wav_of_the_duplicates=data.field('wav_of_the_duplicates')
RA_of_the_duplicates=data.field('RA_others')
DEC_of_the_duplicates=data.field('DEC_others')
"""


#SBE_method=False
#DIRECT_SBE=False #THis means we use direct image-based SBE models
BOOKKEEPING=False


wcscor=True ##WCS correction between visits and parent catalog
###OUR IMPROVEMENTS TO THE EXTRACTION USING GRISMCONF V4
if grismconf_calib_version==4:
	tracecor=True #Trace correction derived by DK based on commissioning data
	lambcor=True #Wavelength correction derived by DK based on commissioning data
else:
	stopping_script #make sure you use grismconf V4 calibrations, otherwise tracecor and lambcor don't work

EXPOSURES=['001','002','003','004','005','006','007','008','009','010','011','012']
PRIMDITS=['2','4']  ## <--- VISITGRP=02 or 04
VISITS=[1,2,3,4]


##Example MANUAL IDs 
# #BROAD HALPHA EMITTERS ID v4
#IDlist=[12446,14947,15157,16221]
#RAlist=[15.048243624454152,15.045550156552324,15.030262274046956,15.034037901840078]
#DEClist=[28.009717338065087,28.028277304725833,28.050176515256226,28.051578609337152]

#FOR RONGMON
IDlist=['1030800']
RAlist=[15.0528963]
DEClist=[28.0416036]

### Reference files for offset correction
### Spectral trace "y-offset" map (calibrated using PID 1076)
eiger_reference_dir = '/scratch/kashinod/EIGER/eiger_reference_files/'
yoffmap_ver = 'v20230109'

fil = eiger_reference_dir+'yoffset_polyreg_F356W.R.ModA_jw01076101_'+yoffmap_ver+'.npy'
print('Spectral trace y-offset map for Mod A: ', fil)
yoffmap_A = np.load(fil, allow_pickle=True)[()]

fil = eiger_reference_dir+'yoffset_polyreg_F356W.R.ModB_jw01076101_'+yoffmap_ver+'.npy'
print('Spectral trace y-offset map for Mod B: ', fil)
yoffmap_B = np.load(fil, allow_pickle=True)[()]

yoffset_map = {'a':yoffmap_A, 
               'b':yoffmap_B}

### lambda offset (calibrated using PID 1076)
lamboff_ver = 'v20221004'
fil = eiger_reference_dir+'lambda_offset_'+lamboff_ver+'.npy'
lambda_offset = np.load(fil, allow_pickle=True)[()]
print('Lambda-offset: ', fil, lambda_offset)

"""
### WCS offset parameters for each Visit/Module
dirimg_vm_dict = {}
for vis in [1,2,3,4]:
    for mod in ['a','b']:
        dirimg_vm_fil = ('/scratch/kashinod/EIGER/J0100/reduction_imaging/calibrated_img3noskymatch/'+
                         'img3_F356W_visit'+str(vis)+mod+'_pdit123_i2d.fits')
        print('Get WCS from ', dirimg_vm_fil)
        dirimg_vm_header = fits.getheader(dirimg_vm_fil, 1)
        dirimg_vm_wcs = WCS(dirimg_vm_header)
#         dirimg_vm_pixscale = (utils.proj_plane_pixel_scales(dirimg_vm_wcs)[0] * u.deg).to(u.arcsec)
#         print('Pixel scale: ', dirimg_vm_pixscale)
        dirimg_vm_dict[str(vis)+mod]={'wcs':dirimg_vm_wcs}

#         ### Offset is defined as the offset to the sources' positions measured in the "raw" WCS in each Visit/Module
#         ### from the positions mesured in the Gaia-calibrated WCS
#         fil = '/scratch/kashinod/EIGER/J0100/reduction_imaging/checkWCS_v2/wcs_offset_params_visit'+str(vis)+mod+'.txt'
#         params = np.loadtxt(fil)        
#         dirimg_vm_dict[str(vis)+mod]={'wcs':dirimg_vm_wcs, 
#                                       'pixscale':dirimg_vm_pixscale,
#                                       'params':params}
        """
        
### photmjsr
photmjsr={'a':0.381, 'b':0.372}

### Grism SENSITIVITY
sens_A_fil = os.path.join(eiger_reference_dir, 'NIRCam.F356W.R.A.1st.sensitivity.EIGER.fits')
sens_B_fil = os.path.join(eiger_reference_dir, 'NIRCam.F356W.R.B.1st.sensitivity.EIGER.fits')
sens = {'a':fits.getdata(sens_A_fil,1),
        'b':fits.getdata(sens_B_fil,1)}



def run(qqq):

    #DEFINE THE FINAL WAVELENGTH GRID
    dwav = 0.000975  #DWAV in Angstrom
    bins= 3.00 + dwav * np.arange(1240) 

    start = time.time()
    print('_______________________________________________________________________________')
    print('___> Start: ', start)

    ID=IDlist[qqq]
    RA0=RAlist[qqq]   ### RA0 and DEC0 should be the accurate position of the source in the coadded direct image 
    DEC0=DEClist[qqq]
    print('___> ID (NUMBER) ', ID, '  (RA0,DEC0)', (RA0, DEC0))

    dd,cc,ww,mm,ee,ll,emem,contcont,qq,nexp,visitcounter,modulecounter=[],[],[],[],[],[],[],[],[],[],[],[]

    for visitnum in range(len(VISITS)): 
        visit=VISITS[visitnum]

        for module in ['a','b']:
            vm = str(visit)+module.lower()
            print('Visit: ', visit, '  Module: ', module)

            ### Correct WCS offset
            if wcscor:
                print('___> Correct WCS')
#                 RA, DEC = eiger_tracing.radec_in_this_vismod_with_wcs(RA0, DEC0, 
#                                                              dirimg_vm_dict[vm]['wcs'],
#                                                              visit=visit, module=module, field=field)
#                 RA = RA[()]
#                 DEC = DEC[()]
#                 print('___> ', (RA0,DEC0), ' ==>', (RA,DEC), ' with wcs')
                RA, DEC = eiger_tracing.radec_in_this_vismod(RA0, DEC0, 
                                                             visit=visit, module=module, field=field)
                RA = RA[()]
                DEC = DEC[()]
                print('___> ', (RA0,DEC0), ' ==>', (RA,DEC))
                
                

            else:
                print('No correction for WCS')
                RA = RA0
                DEC = DEC0
                
            
            for repeat in PRIMDITS:
                
                for qj in EXPOSURES:
                    #print('Extracting',qq,'id',ID) 

                    scidataname = ('/scratch/kashinod/EIGER/J0100/reduction/'+
                                   'calibrated_img2_wfss_globalskySubtV4_emlineMasked_globalMedSubt/'+
                                   'jw0124300100%s_0%s101_00%s_nrc%slong_cal_globalMedSubt.fits'                            
                                   %(visit,repeat,qj,module))
                    
                    emdataname = ('/scratch/kashinod/EIGER/J0100/reduction/'+
                                   'calibrated_img2_wfss_globalskySubtV4_emlineMasked_globalMedSubt_contSubt_kx51_9_Skx21_5/'+
                                  'jw0124300100%s_0%s101_00%s_nrc%slong_cal_globalMedSubt_emline.fits'
                                  %(visit,repeat,qj,module))
                    
                    contdataname = ('/scratch/kashinod/EIGER/J0100/reduction/'+
                                   'calibrated_img2_wfss_globalskySubtV4_emlineMasked_globalMedSubt_contSubt_kx51_9_Skx21_5/'+
                                   'jw0124300100%s_0%s101_00%s_nrc%slong_cal_globalMedSubt_continua.fits'
                                   %(visit,repeat,qj,module))
              
                    ### Spec2 assign_wcs
                    ratename = ('/scratch/kashinod/EIGER/J0100/reduction/calibrated_spc2_bsub/'+
                                'jw0124300100%s_0%s101_00%s_nrc%slong_assign_wcs.fits'%(visit,repeat,qj,module))
                    
                    ratename='/scratch/EIGER/reduce_pipeline/reduced/J0100_v1/grism_F356W/jw0124300100%s_0%s101_00%s_nrc%slong_flatfieldstep.fits'%(visit,repeat,qj,module) ##JM. This is only used for WCS, The other assign_wcs doesnt work for my installation. Unclear origin.
                
                    print('___> SCI    : ', scidataname)
                    print('___> EMLINE : ', emdataname)
                    print('___> CONT   : ', contdataname)
                    print('___> For WCS: ', ratename)
                                              
                    ### Get EXT=0 header 
                    h_sci = fits.getheader(scidataname,0)
                    h_em  = fits.getheader(emdataname,0)
                    h_cont = fits.getheader(contdataname,0)

                    ### Get Grism WCS from "header"
                    grism_data_wcs = WCS(fits.getheader(scidataname,1))
                    
                    ### Get datamodel
                    grism_wcs_dmdl = datamodels.open(ratename)

                    ### PHOTOM Conversion factor
                    sci_photom_conversion_factor = 1.0
                    em_photom_conversion_factor = 1.0
                    cont_photom_conversion_factor = 1.0
                    
                    if 'S_PHOTOM' in h_sci.keys():
                        if h_sci['S_PHOTOM']=='COMPLETE':
                            sci_photom_conversion_factor = photmjsr[module]

                    if 'S_PHOTOM' in h_em.keys():
                        if h_em['S_PHOTOM']=='COMPLETE':
                            em_photom_conversion_factor = photmjsr[module]

                    if 'S_PHOTOM' in h_cont.keys():
                        if h_cont['S_PHOTOM']=='COMPLETE':
                            cont_photom_conversion_factor = photmjsr[module]

                    print('___> SCI/EMLINE/CONT photom_coversion_factor: ', 
                          sci_photom_conversion_factor,
                          em_photom_conversion_factor,
                          cont_photom_conversion_factor)
                    
                    try: 
                        ### Extract spectra
                        ### Not corrected yet for SENSITIVITY
                        ### yoffset must be 0.0 (the correction will be applied explicitly later.)
                        ### Must be use_wcs=False (if True, wcs=WCS(scidata.header) will be used.)
                        d, e, q, l, c, m, w, y,em,cont,C=eiger_tracing.extract_dk(RA,DEC,
                                                                                  scidataname,
                                                                                  emdataname,
                                                                                  contdataname,
                                                                                  ratename,
                                                                                  yoffset=0.0,
                                                                                  yhsize=ysize,
                                                                                  use_wcs=False,
                                                                                  senscorrect=False)
                        print('elllll: ', l)
                        
                      
                        ### Correct trace offset
                        if tracecor: 
                            print('___> Correct trace')
                            x_grism, y_grism = grism_data_wcs.all_world2pix(RA,DEC,0)
                            print('___> x_grism, y_grism: ', x_grism, y_grism)
                            yoffset = eiger_tracing.get_yoffset(x_grism, y_grism, 
                                                                module, 
                                                                yoffset_map=yoffset_map[module])
                            print('___> Correct trace yoffset: ', yoffset)
                            y = y - yoffset

                        ### Correct lambda offset
                        if lambcor:
                            print('___> Correct lambda offset: ', lambda_offset[module.upper()], ' [um]')
                            l = l - lambda_offset[module.upper()]
                      
                        ## MJy/SR --> count rate (if PHOTOM is COMPLETE)
                        d=d/sci_photom_conversion_factor
                        e=e/sci_photom_conversion_factor
                        em=em/em_photom_conversion_factor
                        cont=cont/cont_photom_conversion_factor

                        ## Sensitivity correction (using our own slightly-modified version)
                        ## count rate --> 1e-18 erg/s/cm2/A
                        s = np.interp(l, sens[module]['WAVELENGTH'], sens[module]['SENSITIVITY']) * 1e-18
                        s[s==0]=np.nan
                        d = d/s
                        e = e/s
                        em = em/s
                        cont = cont/s
                      
                        #creating some dummies
                        fakeq=np.zeros(np.shape(q))
                        fakeerr=np.zeros(np.shape(e))+1.

                        #masking bad regions
                        ok = (np.isfinite(d)) * (np.isfinite(l)) * (q==0.)
                        d[~ok]=np.nan
                        e[~ok]=np.nan
                        em[~ok]=np.nan
                        cont[~ok]=np.nan
                      
                        # SCI
                        shifted,shifted_var,shifted_lamb=eiger_tracing.shift_columns_dk(y,d,e**2,l,q,ysize)
                        scrunchd,scrunchd_var,scrunchd_lamb=eiger_tracing.scrunch_columns_dk(bins,
                                                                                             shifted_lamb,
                                                                                             shifted,
                                                                                             shifted_var,
                                                                                             module,
                                                                                             flambda=True)

                        print('___> nanmedian(scrunchd): ',np.nanmedian(scrunchd))
                        #print('___> scrunchd)
                        dd.append(scrunchd)
                        n=np.zeros(np.shape(scrunchd))
                        n[scrunchd>-1E9]+=1
                        nexp.append(n)

                        ee.append(scrunchd_var)
                        qq.append((-1+scrunchd_var/scrunchd_var)) #just a placeholder


                        #EMLINE
                        shifted,shifted_var,shifted_lamb=eiger_tracing.shift_columns_dk(y,em,e**2,l,q,ysize)
                        scrunchd,scrunchd_var,scrunchd_lamb=eiger_tracing.scrunch_columns_dk(bins,
                                                                                             shifted_lamb,
                                                                                             shifted,
                                                                                             shifted_var,
                                                                                             module,flambda=True)
                        emem.append(scrunchd)


                        #CONT
                        shifted,shifted_var,shifted_lamb=eiger_tracing.shift_columns_dk(y,cont,e**2,l,q,ysize)
                        scrunchd,scrunchd_var,scrunchd_lamb=eiger_tracing.scrunch_columns_dk(bins,
                                                                                             shifted_lamb,
                                                                                             shifted,
                                                                                             shifted_var,
                                                                                             module,
                                                                                             flambda=True)
                        contcont.append(scrunchd)

                        #OPTIMAL EXTRACTION WEIGHT (SIMPLE)
                        #the last number is the sigma of the gaussian in pixels for the opt-weight
                        simple_optweight = eiger_tracing.create_simple_optweight(scrunchd_lamb,scrunchd,1.8) 
                        ww.append(simple_optweight)

                        ll.append(scrunchd_lamb)
                        visitcounter.append(visit)
                        modulecounter.append(module)
                    except:
                        continue


                      
    for attempt in [0]:
       # try: ###remove commented try, except bit in case you are extracting the full catalog and some sources may have no spectral coverage
            dd=numpy.array(dd)

            print('____',dd)
            ww=numpy.array(ww)
            ee=numpy.array(ee) 
            ll=numpy.array(ll)
            emem=numpy.array(emem)
            contcont=numpy.array(contcont)
            qq=numpy.array(qq)
            nn=numpy.array(nexp)
            visitcounter=numpy.array(visitcounter)
            modulecounter=numpy.array(modulecounter)


            hdu = fits.PrimaryHDU()
            hdu.header['RA']=RA0   ##  must be RA0, DEC0
            hdu.header['DEC']=DEC0   ##  must be RA0, DEC0
            hdu.header['SOURCEID']=ID
            hdu.header['BUNIT']='1E18 erg/s/cm2/A'
            hdu.header['CTYPE1']='WAVELENGTH'
            hdu.header['CUNIT1']='Angstrom'
            hdu.header['CRPIX1']=1.0
            hdu.header['CRVAL1']=bins[0]*1E4
            hdu.header['CDELT1']=dwav*1E4
            #hdu.header['Yoff_JM']=yoff

            LONGLIST=[hdu]
            first_selections=[visitcounter>0, 
                              modulecounter=='a',
                              modulecounter=='b',
                              visitcounter==1,
                              visitcounter==2,
                              visitcounter==3,
                              visitcounter==4]
            names=['','A','B']#,'V1','V2','V3','V4']
            stack_m_outliers = 3.0
            outlier_iter=5            
            stack_method='mean'  # mean or median

            for jj in range(len(names)):
                thisname=names[jj]
                this_selection=first_selections[jj]

                #Modeled contamination, source-model and weight
                fw = np.nanmedian(ww[this_selection],axis=0)

                #Wavelength
                fl = np.nanmedian(ll[this_selection],axis=0)

                #Data quality and N of exposures
                fq=np.nansum(qq[this_selection],axis=0)
                fn=np.nansum(nn[this_selection],axis=0)

                #Emission-line image and Continuum
                fem,fe,mask_outliers=eiger_tracing.stack_with_reject_outliers(emem[this_selection], 
                                                                                  ee[this_selection], 
                                                                                  fn, 
                                                                                  m = stack_m_outliers, 
                                                                                  method=stack_method,iterations=outlier_iter)
        
                fcont = eiger_tracing.stack_with_mask_outliers(contcont[this_selection], 
                                                                   mask_outliers, 
                                                                   method=stack_method)

                #NEW, data
                    #fd = np.nansum(dd[this_selection],axis=0)/fn
                fd = eiger_tracing.stack_with_mask_outliers(dd[this_selection],
                                                                mask_outliers,
                                                                method=stack_method)


                LONGLIST.append(fits.ImageHDU(data=fem,header=hdu.header,name='EMLINE%s'%thisname))

                LONGLIST.append(fits.ImageHDU(data=fd,header=hdu.header,name='SCI%s'%thisname))
                LONGLIST.append(fits.ImageHDU(data=fe,header=hdu.header,name='ERR%s'%thisname))
                LONGLIST.append(fits.ImageHDU(data=fcont,header=hdu.header,name='CONT%s'%thisname))
                #LONGLIST.append(fits.ImageHDU(data=fq,header=hdu.header,name='QUALITY%s'%thisname))
                #LONGLIST.append(fits.ImageHDU(data=fn,header=hdu.header,name='NEXP%s'%thisname)
                #LONGLIST.append(fits.ImageHDU(data=fw,header=hdu.header,name='OPT_WEIGHT%s'%thisname))
                #LONGLIST.append(fits.ImageHDU(data=fw-fc,header=hdu.header,name='RES%s'%thisname))
                #LONGLIST.append(fits.ImageHDU(data=fl,header=hdu.header,name='WAV%s'%thisname))
                

            #ADD DIRECT IMAGE
            dirimage,dirheader=eiger_tracing.create_cutout(directimage,RA0,DEC0,2048,181) ## Must be RA0, DEC0
            hdu_direct=fits.ImageHDU(data=dirimage,header=dirheader,name='STAMP')

            LONGLIST.append(hdu_direct)
            #ADD EAZY OUTPUT. Add a function t read it in
            #chi = eazycat_chi[eazy_idx,:]\n",
            #"    pdf = np.exp(-chi/2.)\n",
            #"    pdf = pdf/np.sum(pdf)\n",
            eazycat_zaxis=np.arange(0,15,0.1)
            chi=np.zeros(len(eazycat_zaxis))
            pdf=np.zeros(len(eazycat_zaxis))
            col1 = fits.Column(name='z', format='D', array=eazycat_zaxis)
            col2 = fits.Column(name='chi', format='D', array=chi)
            col3 = fits.Column(name='pdf', format='D', array=pdf)
            coldefs = fits.ColDefs([col1, col2, col3])
            hdu_eazy = fits.BinTableHDU.from_columns(coldefs)
            hdu_eazy.header['EXTNAME']='EAZY_ZPDF'


            #ADD GALAXY DATA
            COLOLOL=[]
            COLOLOL.append(fits.Column(name='ID', format='D', array=[ID]))
            COLOLOL.append(fits.Column(name='RA', format='D', array=[RA0])) ## Must be RA0, DEC0
            COLOLOL.append(fits.Column(name='DEC', format='D', array=[DEC0])) ## Must be RA0, DEC0
            #COLOLOL.append(fits.Column(name='F356W_mag', format='D', array=[DEC])) #### NEED TO UPDATE


            coldefs = fits.ColDefs(COLOLOL)
            hdu_src = fits.BinTableHDU.from_columns(coldefs)
            hdu_src.header['EXTNAME']='SRC_INFO'

            LONGLIST.append(hdu_eazy)
            LONGLIST.append(hdu_src)

            if BOOKKEEPING==True:
                COLOL2=[]
                COLOL2.append(fits.Column(name='ID', format='D', array=[ID]))

                COLOL2.append(fits.Column(name='wavs_with_duplicates', format='A3500', array=[wav_duplicates[qqq][:-1]]))
                COLOL2.append(fits.Column(name='IDs_of_duplicates', format='A3500', array=[IDs_others[qqq][:-1]]))
                COLOL2.append(fits.Column(name='wavs_solution_of_duplicates', format='A3500', array=[wav_of_the_duplicates[qqq][:-1]]))
                COLOL2.append(fits.Column(name='RA_of_duplicates', format='A3500', array=[RA_of_the_duplicates[qqq][:-1]]))
                COLOL2.append(fits.Column(name='DEC_of_duplicates', format='A3500', array=[DEC_of_the_duplicates[qqq][:-1]]))

                coldefs = fits.ColDefs(COLOL2)
                hdu_dup = fits.BinTableHDU.from_columns(coldefs)
                hdu_dup.header['EXTNAME']='DUPLICATE_INFO'
                LONGLIST.append(hdu_dup)


            # new_hdul = fits.HDUList([hdu,hdu6,hdu2, hdu3, hdu11,hdu12, hdu4, hdu5, hdu7, hdu8, hdu9,hdu_direct,hdu_eazy])
            # new_hdul.writeto('SPECTRA/stacked_2D_%s_visit%s_mod_%s.fits'%(ID,visit,module), overwrite=True)

            #ADD THE 1D:
            #RELOAD FULL STACK
            fw = np.nanmedian(ww,axis=0)

            #Wavelength
            fl = np.nanmedian(ll,axis=0)
            fn=np.nansum(nn,axis=0)


            #Emission-line image and Continuum
            fem,fe,mask_outliers=eiger_tracing.stack_with_reject_outliers(emem[visitcounter>0], 
                                                                          ee[visitcounter>0], fn, 
                                                                          m = stack_m_outliers, 
                                                                          method=stack_method,iterations=outlier_iter)

            stack_method='mean'
            fcont=eiger_tracing.stack_with_mask_outliers(contcont[visitcounter>0], 
                                                         mask_outliers,
                                                         method=stack_method)

            #NEW, data
            #fd = np.nansum(dd[visitcounter>0],axis=0)/fn
            fd=eiger_tracing.stack_with_mask_outliers(dd[visitcounter>0], 
                                                      mask_outliers,
                                                      method=stack_method)



            # Hornes optimal extraction of model
            t1 = np.nansum((fd)*fw,axis=0)
            t2 = np.nansum(fw**2,axis=0)
            opt_extracted = t1/t2	


            t1 = np.nansum((fem)*fw,axis=0)
            emline_extracted = t1/t2	

            #find where model has max:
            argmax=int(len(fd[:,0])/2) ##+1 ?
            argmax=26


            lam = np.nanmean(fl[:,:],axis=0)
            print(lam,fl,bins)
            lam = np.nanmean(fl[argmax-5:argmax+5,:],axis=0)
            print(lam,fl,bins, '  All should be identical.') 
            #s = C.SENS["+1"](lam)
            sens_A = np.interp(lam, sens['a']['WAVELENGTH'], sens['a']['SENSITIVITY'])
            sens_B = np.interp(lam, sens['b']['WAVELENGTH'], sens['b']['SENSITIVITY'])


            ####hdu.header['CRVAL1']=lam[0]*1E4


            #get some simple extractions
            boxcar=np.nansum(fd[argmax-5:argmax+5,:],axis=0)
            boxcar_emline=np.nansum(fem[argmax-5:argmax+5,:],axis=0)
            len_box=10
            #errors:
            ivar=1./fe**2

            boxcar_err=(np.nansum(fe[argmax-5:argmax+5,:]**2,axis=0)**0.5)/len_box# this is wrong

            t1 = np.nansum(ivar*fw**2,axis=0)
            t2 = np.nansum(fw,axis=0)
            simple_extracted_err = (t1/t2)**-0.5



            COLS=[]
            COLS.append(fits.Column(name='wavelength',unit='micron',format='E',array=lam))

            COLS.append(fits.Column(name='flux_opt_ext',unit='1E18 erg/s/cm2/A',format='E',array=opt_extracted))
            COLS.append(fits.Column(name='flux_opt_ext_err',unit='1E18 erg/s/cm2/A',format='E',array=simple_extracted_err))

            COLS.append(fits.Column(name='flux_opt_emline',unit='1E18 erg/s/cm2/A',format='E',array=emline_extracted))

            COLS.append(fits.Column(name='flux_boxcar',format='E',array=boxcar))
            COLS.append(fits.Column(name='flux_boxcar_err',format='E',array=boxcar_err))
            COLS.append(fits.Column(name='flux_boxcar_emline',format='E',array=boxcar_emline))

            COLS.append(fits.Column(name='sensitivity',format='E',array=sens_A))
            COLS.append(fits.Column(name='SENSITIVITYA',format='E',array=sens_A))
            COLS.append(fits.Column(name='SENSITIVITYB',format='E',array=sens_B))

            THESE_IDS=[]
            DUPLICATE_MARKERS={}
            duplicate_marker=np.zeros(len(lam))
            DUPLICATE_THEIR_MARKERS={}


            if BOOKKEEPING==True:
                unique_ids=np.unique(IDs_dups)
                for un in unique_ids:
                    DUPLICATE_MARKERS.update({un:np.zeros(len(lam))})
                    DUPLICATE_THEIR_MARKERS.update({un:np.zeros(len(lam))})

                f_x_lamb=interp1d(lam,np.arange(0,len(lam),1))

                for jj in range(len(wav_dups)):
                    this_dup_ID=IDs_dups[jj]
                    this_dup_wav=wav_dups[jj]
                    this_dup_their_wav=wav_of_dups[jj]

                    try:
                        duplicate_marker[int(np.round(f_x_lamb(this_dup_wav)))]+=1

                        DUPLICATE_MARKERS[this_dup_ID][int(np.round(f_x_lamb(this_dup_wav)))]+=1
                        DUPLICATE_THEIR_MARKERS[this_dup_ID][int(np.round(f_x_lamb(this_dup_their_wav)))]+=1
                    except:
                        continue #this is implemented as sometimes dupliactes are outside of the wavelength range of the 1D

                COLS.append(fits.Column(name='n_duplicates',format='E',array=duplicate_marker))
                duplicate_marker[duplicate_marker>0]=1.
                COLS.append(fits.Column(name='FLAG_DUPLICATE',format='E',array=duplicate_marker))

                for bb in range(len(unique_ids)):
                    COLS.append(fits.Column(name='duplicate_id%s'%unique_ids[bb],format='E',array=DUPLICATE_MARKERS[unique_ids[bb]]))
                    COLS.append(fits.Column(name='duplicate_as_in_id%s'%unique_ids[bb],format='E',array=DUPLICATE_THEIR_MARKERS[unique_ids[bb]]))

            #print(len(COLS))
            cols=fits.ColDefs(COLS)#,col13,col14])
            hdu_1D = fits.BinTableHDU.from_columns(cols)
            #hdu=fits.PrimaryHDU(numpy.arange(100.))
            hdu_1D.writeto(FOLDER+'stacked_1D_%s.fits'%(ID),overwrite=True)

            hdu_1D.header['EXTNAME']='1D'

            #LONGLIST.append(hdu_1D)

            #SAVE THE 2D
            print('Saving to:', FOLDER+'stacked_2D_%s_%s.fits'%(field,ID))
            new_hdul = fits.HDUList(LONGLIST)
            new_hdul.writeto(FOLDER+'stacked_2D_%s_%s.fits'%(field,ID), overwrite=True)
            end = time.time()

            print('Time lapsed',end-start,'seconds')

       # except:
#            continue




n_procs = 44 # number of cores available
with Pool(n_procs) as pool:
    pool.map(run,range(len(IDlist)))


# for mynum in [0]:
#     run(mynum)
