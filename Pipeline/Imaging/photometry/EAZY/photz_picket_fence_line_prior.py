import eazy
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.gridspec import GridSpec
from astropy.table import Table, Column, join
import matplotlib.colors as colors
from matplotlib import cm
import matplotlib.patheffects as pe
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)    

##########################################################
def chi_dist_line(obs_wave, chi2_i ,zgrid, line_list, plot=False):
    if plot:
        plt.plot(zgrid, chi2_i, color=tableau20[2], label=None)
        #ax.plot(self.zgrid, prior/prior.max()*pz.max(), color='g',                label='prior')
        plt.fill_between(zgrid, chi2_i, chi2_i*0, color='yellow', alpha=0.5, label=None)

        #plt.vlines(zspec_i, chi2_i.min()*0.95, chi2_i.max()*1.05, color='black', label='true '+trueid)

    #calc wavelength from redshift and true line id
    #obs_wave = obj['Wav_brightest']

    chi2_line = np.zeros(len(line_list), dtype=float)
    #calc implied redshift and plot
    for il, curline in enumerate(line_list):
        #print(il, curline)
        linez = obs_wave/curline[1] - 1.0
        #plt.plot([linez,linez],[chi2_i.min()*0.95, chi2_i.max()*1.05],c=tableau20[il], label=curline[0])
        zid = np.argmin(np.abs(linez-zgrid))
        chi2_line[il] = chi2_i[zid]
        if curline[0].endswith('_2'): continue
        if plot:
            plt.scatter(zgrid[zid], chi2_i[zid], color=tableau20[il], marker='o')
            plt.vlines(linez, chi2_i.min()*0.95, chi2_i.max()*1.05, color=tableau20[il], label=curline[0])

    if plot:
        plt.yscale("log")
        plt.ylim(chi2_i.min()*0.95,chi2_i.max()*1.05)
        plt.xlim(0,zgrid[-1])
            
        plt.xlabel('$z$')
        plt.ylabel(r'$ \chi (z)$')
        #plt.legend()
        #plt.grid()
    return chi2_line
##########################################################
def plot_line_hist(chi2_line, line_list, trueid):
    plt.yscale("log")
    plt.ylim(chi2_line.min()*0.9,chi2_line.max()*1.2)
    plt.xlim(-2,len(line_list)+1)
    line_ticks=[]
    line_labs = []
    for il, curline in enumerate(line_list):
        if curline[0].endswith('_2'): continue
        plt.scatter(il+1, chi2_line[il], color=tableau20[il], marker='o', label=curline[0])
        line_ticks.append(il+1)
        line_labs.append(curline[0])
        if curline[0] == trueid:
            plt.scatter(il+1, chi2_line[il], color='black', marker='x')
    plt.ylabel(r'$ \chi (z)$')    
    plt.xticks(line_ticks, line_labs)
##########################################################
def picket_fence(incat, line_list):
    line_cat = incat.copy()
    #open new cols for data
    #line_cat.add_column(Column(['       ']*len(line_cat), name='true_line'))
    line_cat.add_column(Column([0.0]*len(line_cat), name='obs_wave'))

    for il, curline in enumerate(line_list):
            line_cat.add_column(Column([0.0]*len(line_cat), name='chi2_'+curline[0]))

    line_cat.add_column(Column([0.0]*len(line_cat), name='z_line_best'))
    line_cat.add_column(Column([0.0]*len(line_cat), name='delta_chi2'))
    line_cat.add_column(Column(['       ']*len(line_cat), name='best_line'))


    #loop over objs
    for obj in line_cat:
        #what is observed line
        obs_wave = obj['Wav_brightest']

        #measure chi2 at different line positions
        for il, curline in enumerate(line_list):
            #print(il, curline)
            linez = obs_wave/curline[1] - 1.0
            zid = np.argmin(np.abs(linez-ZGRID))
            obj['chi2_'+curline[0]] = obj['CHI2'][zid]
        #now pick best redshift
        chi2_array = [obj['chi2_BrB'],
                        obj['chi2_BrC'],
                        obj['chi2_PaA'],
                        obj['chi2_PaB'],
                        obj['chi2_HeI'],
                        obj['chi2_SIII_1'],
                        obj['chi2_Ha'],
                        obj['chi2_OIII_1'],
                        obj['chi2_OII']]
        bestline = np.nanargmin(chi2_array)

        #obj['true_line'] = trueline[0]
        obj['obs_wave'] = obs_wave
        obj['z_line_best'] = obs_wave/line_list[bestline][1] - 1.0
        obj['delta_chi2'] = np.sort(chi2_array)[1] - np.sort(chi2_array)[0]
        obj['best_line'] = line_list[bestline][0]
    return line_cat
##########################################################
line_list= [['BrB',    26251.40],
            ['BrC',    21655.20],
            ['PaA',    18751.00],
            ['PaB',    12818.10],
            ['HeI',    10830.34],
            ['SIII_1',  9531.10],
            ['Ha',      6562.82],
            ['OIII_1',  5006.84],
            ['OII',     3727.42]]

#load data
basedir = '/scratch/mruari/EIGER/imaging/J1148+5251/photometry/'
basename = basedir+'EAZY/J1148_photzcat_v2_EAZY_'
zout = Table.read(basename+'_output.zout.fits')
data = fits.open(basename+'_output.data.fits')

sex_cat = Table.read(basedir+"catalogs/J1148+5251_photcat_v2_noisemodel_short.fits") 

#add cols to zout
ZGRID = data['ZGRID'].data
zgrid_col = Column([data['ZGRID'].data]*len(zout),            name='ZGRID')
chi2_col = Column(data['CHI2'].data,                             name='CHI2')
zout.add_columns([zgrid_col, chi2_col])
zout.rename_column('id','NUMBER')

#merge, need to match
merge_cat = join(zout, sex_cat, keys=['NUMBER'], join_type='left', table_names=['zout', 'sex'])

merge_cat.write(basename+'_photoz_merged.fits', overwrite=True)

# #line catalog
# linecat = Table.read('/net/hyperion/scratch/EIGER/forsharing/Halpha_MetalHosts_J0100/J0100_photcat_v2_earlyHAEs_withcommentsJM.fits')
# linecat.keep_columns(['NUMBER', 'Nlines_in_both_dispersions', 'Flux_brightest', 'Wav_brightest', 'z_Halpha', 'comments_JM', 'X_line',])

# #merge
# master_cat = join(linecat, merge_cat, keys=['NUMBER'], join_type='left', table_names=['zout', 'linecat'])

# line_dets = picket_fence(master_cat, line_list)

# #save cat
# line_dets.write(basename+'_fence_line_prior.fits', overwrite=True)




# obj = line_dets[(line_dets['NUMBER'] == 9950)][0]
# chi_line = chi_dist_line(obj['Wav_brightest'],obj['CHI2'] ,obj['ZGRID'], line_list, plot=True)
# plt.show()