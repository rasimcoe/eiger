from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.Qt import QStandardItemModel, QStandardItem, QGuiApplication
from PyQt5.QtWidgets import QFileDialog, QMenu
from gui import Ui_Dialog
from SpecGui import Ui_SpectrumSelector
from LineGui import Ui_LineSelector
from vpModelGui import Ui_VoigtProfileModel
from eiger.QuasarSpec.loadQsoSpec import loadQsoSpec
from astropy.table import Table, unique
import matplotlib.pyplot as plt
from pypeit.core.wave import airtovac
import astropy.units as u
import numpy as np
import pickle

import eiger.QuasarSpec.loadQsoSpec as spec
from pypeit.core.wave import airtovac
from mcvp import model as vm
from mcvp import vfit as vf
from astropy.convolution import convolve, Gaussian1DKernel

import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)

class GuiProgram(Ui_Dialog):

    def __init__(self, dialog):

        self.xmin = 0.0
        self.xmax = 10.0
        self.ymin = 0.0
        self.ymax = 40.0
        self.cid  = -1
        self.zabs            = None
        self.linelist        = None
        self.instruments     = None
        self.doplot          = [True,True,True]
        self.plotinstruments = ['XSH_NIR','HIRES','XSH_VIS']
        self.w               = None
        self.idtable         = Table(names=('redshift','ion','restwv'),dtype=('f4','S2','f4'))
        self.components      = None
        self.regionleftbound = 0.0
        self.regionrightbound = 0.0
        self.regioninstrument = ' '
        self.prof_guess = None
        self.profs_fire = None
        self.profs_hires = None
        self.profs_xsh_vis = None
        self.profs_xsh_nir = None
        self.vp_model = None
        self.vpTree = None
        
        ''' This method gets called when the window is created. '''
        Ui_Dialog.__init__(self)              # Initialize Window
        self.setupUi(dialog)                  # Set up the UI
        # Initialize the figure in our window
        figure = Figure()                     # Prep empty figure
#        axis = figure.add_subplot(111)        # Prep empty plot
        axis = figure.subplots(3,1,sharex=True,gridspec_kw={'hspace':0})        # Prep empty plot
        self.initialize_figure(figure, axis)  # Initialize!

        self.eigerObjectSelect.addItem("1: J1030+1524")
        self.eigerObjectSelect.addItem("2: P159-02")
        self.eigerObjectSelect.addItem("3: J1120+0641")
        self.eigerObjectSelect.addItem("4: J0100+2802")
        self.eigerObjectSelect.addItem("5: J1148+5251")
        self.eigerObjectSelect.addItem("6: J0148+0600")
        
        # Connect our button with plotting function
        self.reloadPushButton.clicked.connect(self.load_newplot)
        self.SpectraSelectionPushButton.clicked.connect(self.choose_spectra)
        self.vpTreePushButton.clicked.connect(self.show_vptree)
        self.loadVPfitButton.clicked.connect(self.loadVPfit)
        self.canvas.mpl_connect('key_press_event',self.on_keypress)
        self.canvas.mpl_connect('button_press_event',self.on_buttonpress)
        self.canvas.mpl_connect('button_release_event',self.on_buttonrelease)
        self.canvas.mpl_connect('pick_event',self.onpick)
        self.redshiftEdit.returnPressed.connect(self.set_redshift)
        self.load_linelist()
        
    def change_plot(self):
        ''' Plots something new in the figure. '''

        # Clear whatever was in the plot before
        nplots = len(self.ax)
        for i in range(nplots):
            self.ax[i].remove()

        # Set up the new plot

        nplots = sum(self.doplot)
        self.ax = self.fig.subplots(nplots,sharex=True,squeeze=True,gridspec_kw={'hspace':0})
            
        for i in range(len(self.doplot)):

            if (self.doplot[i]):

                instrument = self.plotinstruments[i]

                xdata = self.spec[instrument]['wave']
                ydata = self.spec[instrument]['flux']/self.spec[instrument]['cont']
                yerr = 1.0/np.sqrt(self.spec[instrument]['ivar'])/self.spec[instrument]['cont']
                
                # Plot data, add labels, change colors, ...
                self.ax[i].set_xlabel('Wavelength (AA)')
                self.ax[i].set_ylabel('Flux')
                self.ax[i].text(self.xmin+0.05*(self.xmax-self.xmin),\
                                self.ymin+0.85*(self.ymax-self.ymin),instrument)
                self.ax[i].set_xlim(self.xmin,self.xmax)
                self.ax[i].set_ylim(self.ymin,self.ymax)
                self.ax[i].step(xdata,ydata,where='mid')
                self.ax[i].plot(xdata,yerr,color='r',alpha=0.25)
                self.ax[i].plot([min(xdata),max(xdata)],[0,0],color='k',alpha=0.33)

        if (self.zabs != None or len(self.idtable)>0):
            self.label_abslines()

        if (self.profs_fire != None):
            self.plotVPfits()

        self.plotVPGuess()
            
        # Make sure everything fits inside the canvas
        self.fig.tight_layout()
        # Show the new figure in the interface
        self.canvas.draw()

    def initialize_figure(self, fig, ax):
        ''' Initializes a matplotlib figure inside a GUI container.
            Only call this once when initializing.
        '''
        # Figure creation (self.fig and self.ax)
        self.fig = fig
        self.ax = ax
        # Canvas creation
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setFocusPolicy( QtCore.Qt.ClickFocus )
        self.canvas.setFocus()
        self.plotLayout.addWidget(self.canvas)
        self.canvas.draw()
        # Toolbar creation
        self.toolbar = NavigationToolbar(self.canvas, self.plotWindow,
                                         coordinates=True)
        self.plotLayout.addWidget(self.toolbar)
        
    def load_newplot(self):

        selected_object = self.eigerObjectSelect.currentText()
        print(selected_object)

        indx = int(selected_object.split(':')[0])

        self.spec = loadQsoSpec(indx,revision='current')
        self.instruments = list(self.spec.keys())[4:]

        self.plotinstruments = self.instruments[0:3]

        # The J1030 HIRES reduction is in air wavelengths
        if (indx == 1):
            self.spec['HIRES']['wave'] = np.array(airtovac(self.spec['HIRES']['wave']*u.AA))
            
        self.xmin = 8000
        self.xmax = 20000
        self.ymin = -1.0
        self.ymax = 3.0
        self.change_plot()

    def on_buttonpress(self,event):
        # On shift left-click, define a line fitting region
        if (event.button == 1):
            mods = QGuiApplication.queryKeyboardModifiers()
            if (mods == QtCore.Qt.ShiftModifier):
                wave = event.xdata
                instrument = event.inaxes
                self.regionleftbound = wave

    def on_buttonrelease(self,event):

        # On right click, open line ID GUI
        if (event.button == 3):
            wave = event.xdata
            self.set_abslineid(wave)

        # On shift left-click, define a line fitting region
        if (event.button == 1):
            mods = QGuiApplication.queryKeyboardModifiers()
            if (mods == QtCore.Qt.ShiftModifier):
                wave = event.xdata
                instrument = event.inaxes
                self.regionrightbound = wave
                for i in range(len(self.ax)):
                    a = self.ax[i]
                    inst = self.plotinstruments[i]
                    if (a == instrument):
                        self.regioninstrument = inst
                        self.vpTree.addFitRegion(inst,self.regionleftbound,self.regionrightbound)
                        
                        
    def on_keypress(self,event):
        
        xmnx = self.ax[0].get_xlim()
        ymnx = self.ax[0].get_ylim()
        self.xmin = xmnx[0]
        self.xmax = xmnx[1]
        self.ymin = ymnx[0]
        self.ymax = ymnx[1]

        redraw = True
        
        if (event.key == 'l'):
            self.xmin=event.xdata
        elif (event.key == 'r'):
            self.xmax=event.xdata
        elif (event.key == 't'):
            self.ymax=event.ydata
        elif (event.key == 'b'):
            self.ymin=event.ydata
        elif (event.key == 'w'):
            self.xmin = min(self.spec['FIRE']['wave'])
            self.xmax = max(self.spec['FIRE']['wave'])
            self.ymin = -5
            self.ymax = 5
        elif (event.key == '}' or event.key == ']'):
            dx = self.xmax-self.xmin
            self.xmin = self.xmax
            self.xmax += dx
            self.ymin = ymnx[0]
            self.ymax = ymnx[1]
        elif (event.key == '{' or event.key == '['):
            dx = self.xmax-self.xmin
            self.xmax = self.xmin
            self.xmin -= dx
            self.ymin = ymnx[0]
            self.ymax = ymnx[1]
        elif (event.key == 'o'):
            xcen = (self.xmin+self.xmax)/2.0
            dx   = self.xmax - xcen
            self.xmin = xcen - 1.5*dx
            self.xmax = xcen + 1.5*dx
        elif (event.key == 'i'):
            xcen = (self.xmin+self.xmax)/2.0
            dx   = self.xmax - xcen
            self.xmin = event.xdata - 0.6*dx
            self.xmax = event.xdata + 0.6*dx
        elif (event.key == 'C'):
            wave0 = event.xdata        
            self.check_lineid(wave0,'CIV')
            redraw = False
        elif (event.key == 'M'):
            wave0 = event.xdata        
            self.check_lineid(wave0,'MgII')
            redraw = False
        elif (event.key == 'S'):
            wave0 = event.xdata        
            self.check_lineid(wave0,'CII*')
            redraw = False
        elif (event.key == 'F'):
            self.idtable.write("J0100_idtable.dat",format='ascii.fixed_width')
            print("Writing out ASCII table")
            redraw = False
        elif (event.key == 'R'):
            self.idtable = Table.read("J0100_idtable.dat",format='ascii.fixed_width')
        elif (event.key == 'c'):
            wave_marked = event.xdata
            obswaves    = (1+self.idtable['redshift']) * self.idtable['restwv']
            diffs = abs(wave_marked - obswaves)
            line = self.idtable[diffs == min(diffs)]
            zz = float(wave_marked / line['restwv'] - 1.0)
            # Mark an absorption "component"
            cmd = f"m.addcomponent({zz:6.4f},bpriors=[2,100])"
            self.vpTree.add_component(round(zz,4))
            print(cmd)
            redraw = False
        elif (event.key == 'I'):
            wave_marked = event.xdata
            obswaves    = (1+self.idtable['redshift']) * self.idtable['restwv']
            ions = self.idtable['ion']
            diffs = abs(wave_marked - obswaves)
            line = self.idtable[diffs == min(diffs)]
            ionname = (ions[diffs == min(diffs)])[0]
            zmatch = (ions[diffs == min(diffs)])[0]
            # Mark an absorption "Ion"
            zz = float(wave_marked / line['restwv'] - 1.0)
            print(f"m.addion({zz:6.4f}, \'{ionname}\', N=14)")
            self.vpTree.add_ion(round(zz,4),ionname)
            redraw = False
        elif (event.key == 'T'):
            # Mark an absorption "Transition"
            wave_marked = event.xdata
            obswaves    = (1+self.idtable['redshift']) * self.idtable['restwv']
            ions = self.idtable['ion']
            restwvs = self.idtable['restwv']
            diffs = abs(wave_marked - obswaves)
            line = self.idtable[diffs == min(diffs)]
            ionname = (ions[diffs == min(diffs)])[0]
            restwv = int(np.floor((restwvs[diffs == min(diffs)])[0]))
            # Mark an absorption "Ion"
            zz = float(wave_marked / line['restwv'] - 1.0)
            print(f"m.addtransition({restwv},\'{ionname}\',{zz:6.4f})")
            self.vpTree.add_transition(round(zz,4),ionname,restwv)
            redraw = False            
        else:
            self.xmin = xmnx[0]
            self.xmax = xmnx[1]
            self.ymin = ymnx[0]
            self.ymax = ymnx[1]
            redraw = False

        if(redraw):
            self.change_plot()

    def check_lineid(self, wave0, ion):
        if (ion == 'CIV'):
            wave1 = wave0 * 1550.77845 / 1548.2049           
            z = wave0 / 1548.2049 - 1
            print(f"CIV: z = {z}")
        elif (ion == 'MgII'):
            wave1 = wave0 * 2803.5314853 / 2796.3542699
            z = wave0 / 2796.354 - 1
            print(f"MgII: z = {z}")
        elif (ion == 'CII*'):
            wave1 = wave0 * 1335.7077 / 1334.5323
            z = wave0 / 1334.5323 - 1
            print(f"CII*: z = {z}")
        nplots = len(self.ax)
        for i in range(nplots):
            self.ax[i].plot([wave0,wave0,wave1,wave1],[1,1.5,1.5,1],color='r')
        self.canvas.draw()
            
    def set_redshift(self):
        try:
            entry = float(self.redshiftEdit.text())
        except:
            print("ERROR: Redshift entry is not a valid number")
            entry = None
            self.redshiftEdit.setText('')
            return()
            
        if (entry < 0 or entry > 10):
            print("Error: redshift mst be between 0 and 10")
            entry = None
            self.redshiftEdit.setText('')
            return()
            
        self.zabs = entry
        self.change_plot()

    def load_linelist(self):
        self.linelist = Table.read('../atomic_data.txt',format='ascii')

    def label_abslines(self):

        if (self.zabs != None):
            z = self.zabs
            for line in self.linelist:
                obswave = line['wave'] * (1+z)
                ion = f"{line['ion']} {line['wave']:4.0f}"
                if (obswave > self.xmin and obswave < self.xmax):
                    nplots = len(self.ax)
                    for i in range(nplots):
                        self.ax[i].plot([obswave,obswave],[-100,100],alpha=0.5,color='k')
                        self.ax[i].text(obswave,1.6,ion,rotation='vertical',horizontalalignment='center',\
                                        picker=True, fontsize=8)
                    
        for line in self.idtable:
            obswave = line['restwv'] * (1+line['redshift'])
            ion = f"{line['ion']} {line['restwv']:4.0f} ({line['redshift']:4.3f})"
            if (obswave > self.xmin and obswave < self.xmax):
                nplots = len(self.ax)
                for i in range(nplots):
                    self.ax[i].plot([obswave,obswave],[-100,100],alpha=0.15,color='g')
                ylims = self.ax[0].get_ylim()
                self.ax[0].text(obswave,1.1*ylims[1],ion,rotation='vertical',horizontalalignment='center',\
                                picker=True, color='k', fontsize=7)
                    
                    
    def choose_spectra(self):
        self.selector = QtWidgets.QDialog()
        self.specDialog = SpecSelect(self.selector)

        choices = self.instruments        
        for instrument in choices:
            self.specDialog.spec1ComboBox.addItem(instrument)
            self.specDialog.spec2ComboBox.addItem(instrument)
            self.specDialog.spec3ComboBox.addItem(instrument)

        self.specDialog.spec1ComboBox.setCurrentIndex\
            (self.specDialog.spec1ComboBox.findText(self.plotinstruments[0]))
        self.specDialog.spec2ComboBox.setCurrentIndex\
            (self.specDialog.spec2ComboBox.findText(self.plotinstruments[1]))
        self.specDialog.spec3ComboBox.setCurrentIndex\
            (self.specDialog.spec3ComboBox.findText(self.plotinstruments[2]))

        if(self.doplot[1]):
            self.specDialog.Spec2CheckBox.setChecked(True)
        if(self.doplot[2]):
            self.specDialog.Spec3CheckBox.setChecked(True)
        
        self.selector.exec()

        if (self.specDialog.yesno == True):
            self.plotinstruments = self.specDialog.plotlist
            self.doplot[1] = self.specDialog.Spec2CheckBox.isChecked()
            self.doplot[2] = self.specDialog.Spec3CheckBox.isChecked()
            
        self.change_plot()
        
    def set_abslineid(self, wave):
        self.selector   = QtWidgets.QDialog()
        self.lineDialog = LineSelect(self.selector)

        linelist = Table.read('../atomic_data.txt',format='ascii')

        for line in linelist:
            lab = f"{line['ion']}\t{line['wave']:6.2f}"
            self.lineDialog.lineListWidget.addItem(lab)

        self.selector.exec()
        self.register_z(wave)
        #self.lineDialog.buttonBox.clicked.connect(lambda: self.register_z(wave))

    def register_z(self,wave):
        # print(f"Register: {wave}, {self.lineDialog.linewave}")
        self.zabs = float(wave) / float(self.lineDialog.linewave) - 1.0
        # print(f"Redshift: {self.zabs}")
        self.label_abslines()
        self.redshiftEdit.setText(f"{self.zabs:5.4f}")
        self.change_plot()

    def onpick(self,event):
        txt = event.artist
        obswv         = txt.get_position()[0]
        ion           = txt.get_text().split()[0]
        restwv_approx = txt.get_text().split()[1]
        redshift = self.zabs

        for testline in self.linelist['wave']:
            if(np.abs(float(testline) - float(restwv_approx)) < 2.0):
                restwv = testline
                break
        
        newrow = (redshift,ion,restwv)
        self.idtable.add_row(newrow)
        print(self.idtable)

    def show_vptree(self):
        try:
            self.tree.show()
        except:
            self.tree = QtWidgets.QDialog()
            self.vpTree = VPModelTree(self.tree, self)
            self.tree.show()

    def loadVPfit(self):

        # Get the file name from the user
        options = QFileDialog.Options()
        picklefile, _ = QFileDialog.getOpenFileName(None,"Open","","All Files (*);;Text Files (*.txt)", options=options)

        if (picklefile == None or picklefile == ''):
            return()
        
        with open(picklefile, "rb") as fp:
            vpfit = pickle.load(fp)

        m       = vpfit['model']
        samples = vpfit['samples']

        fire_kernel     = Gaussian1DKernel(stddev=4.0/2.355)
        specobj_fire    = vm.Spectrum(self.spec['FIRE']['wave'],self.spec['FIRE']['flux']/self.spec['FIRE']['cont'], \
                                      1/np.sqrt(self.spec['FIRE']['ivar'])/self.spec['FIRE']['cont'],fire_kernel,\
                                      lines=[1548,1550,2796,2803,1526,1393,1402,1334,1335,2600,2586,2382,2374,2344,1670,5891,5897,2852,1854,1862,1304,1302,1260])
            
        hires_kernel     = Gaussian1DKernel(stddev=3.0/2.355)
        specobj_hires    = vm.Spectrum(self.spec['HIRES']['wave'],self.spec['HIRES']['flux']/self.spec['HIRES']['cont'], \
                                       1/np.sqrt(self.spec['HIRES']['ivar'])/self.spec['HIRES']['cont'],hires_kernel,\
                                       lines=[1548,1550,2796,2803,1526,1393,1402,1334,1335,2600,2586,2382,2374,2344,1670,5891,5897,2852,1854,1862,1304,1302,1260])

        xsh_nir_kernel  = Gaussian1DKernel(stddev=4.2/2.355)
        xsh_nir_kernel  = Gaussian1DKernel(stddev=2.2/2.355)
        specobj_xsh_nir = vm.Spectrum(self.spec['XSH_NIR']['wave'],self.spec['XSH_NIR']['flux']/self.spec['XSH_NIR']['cont'], \
                                      1/np.sqrt(self.spec['XSH_NIR']['ivar'])/self.spec['XSH_NIR']['cont'],xsh_nir_kernel,\
                                      lines=[1548,1550,2796,2803,1526,1393,1402,1334,1335,2600,2586,2382,2374,2344,1670,5891,5897,2852,1854,1862,1304,1302,1260])
        
        xsh_vis_kernel  = Gaussian1DKernel(stddev=4.8/2.355)
        specobj_xsh_vis = vm.Spectrum(self.spec['XSH_VIS']['wave'],self.spec['XSH_VIS']['flux']/self.spec['XSH_VIS']['cont'], \
                                      1/np.sqrt(self.spec['XSH_VIS']['ivar'])/self.spec['XSH_VIS']['cont'],xsh_vis_kernel,\
                                      lines=[1548,1550,2796,2803,1526,1393,1402,1334,1335,2600,2586,2382,2374,2344,1670,5891,5897,2852,1854,1862,1304,1302,1260])
        
        self.profs_fire    = vf.sampleVPFits(m,specobj_fire,samples[1000:],50)
        self.profs_hires   = vf.sampleVPFits(m,specobj_hires,samples[1000:],50)
        self.profs_xsh_nir = vf.sampleVPFits(m,specobj_xsh_nir,samples[1000:],50)
        self.profs_xsh_vis = vf.sampleVPFits(m,specobj_xsh_vis,samples[1000:],50)

        if (False):
            voigt_profiles = {'FIRE':self.profs_fire, \
                              'HIRES': self.profs_hires, \
                              'XSH_VIS': self.profs_xsh_vis, \
                              'XSH_NIR': self.profs_xsh_nir}
            with open("SDSS1030_z5.7_highions.pickle", "wb") as fp:
                pickle.dump(voigt_profiles,fp, pickle.HIGHEST_PROTOCOL)
        
        # Now get the statistics
        with open(picklefile,'rb') as fp:
            thefit = pickle.load(fp)
            mm      = thefit['model']
            samples = thefit['samples'][1000:]

        thetable = []

        linetable = Table(names=['ion','z_median','z_16pct','z_84pct','b_median','b_16pct','b_84pct','N_median','N_16pct','N_84pct'],\
                          dtype=('U1','f8','f8','f8','f8','f8','f8','f8','f8','f8'))
        linetable['z_median'].format = '7.5f'
        linetable['z_16pct'].format = '7.5f'
        linetable['z_84pct'].format = '7.5f'
        linetable['b_median'].format = '3.1f'
        linetable['b_16pct'].format = '3.1f'
        linetable['b_84pct'].format = '3.1f'
        linetable['N_median'].format = '5.2f'
        linetable['N_16pct'].format = '5.2f'
        linetable['N_84pct'].format = '5.2f'
        
        indx = 0
        for c in mm.components:
            comp = mm.components[c]
            redshift = comp.z
            z_low,z_med,z_high = np.quantile(samples[:,indx],[0.16,0.5,0.84])
            indx += 1
            b_turb   = comp.b_turb
            b_low,b_med,b_high = np.quantile(samples[:,indx],[0.16,0.5,0.84])
            indx += 1
            for ion in comp.ions:
                thision  = comp.ions[ion]
                ionname = thision.name
                Nlow,Nmed,Nhigh = np.quantile(samples[:,indx],[0.16,0.5,0.84])
                indx += 1
                #print(f"{z_med:7.5f} [{z_low:7.5f},{z_high:7.5f}]\t{b_med:3.1f} [{b_low:3.1f},{b_high:3.1f}]\t{ionname}\t\t{Nmed:5.2f} [{Nlow:5.2f},{Nhigh:5.2f}]")
                linetable.add_row([ion,z_med,z_low,z_high,b_med,b_low,b_high,Nmed,Nlow,Nhigh])
                
        linetable.sort('z_median')
        linetable.reverse()
        print(linetable)

    def plotVPGuess(self):

        linelist = [1548,1550,2796,2803,1526,1393,1402,1334,1335,2600,2586,2382,2374,2344,1670,5891,5897,2852,1854,1862,1304,1302,1260]

        if ('FIRE' in self.instruments):
            fire_kernel     = Gaussian1DKernel(stddev=4.0/2.355)
            specobj_fire    = vm.Spectrum(self.spec['FIRE']['wave'],self.spec['FIRE']['flux']/self.spec['FIRE']['cont'], \
                                          1/np.sqrt(self.spec['FIRE']['ivar'])/self.spec['FIRE']['cont'],fire_kernel,\
                                          lines=linelist)

        if ('HIRES' in self.instruments):
            hires_kernel     = Gaussian1DKernel(stddev=3.0/2.355)
            specobj_hires    = vm.Spectrum(self.spec['HIRES']['wave'],self.spec['HIRES']['flux']/self.spec['HIRES']['cont'], \
                                           1/np.sqrt(self.spec['HIRES']['ivar'])/self.spec['HIRES']['cont'],hires_kernel,\
                                           lines=linelist)
            
        if ('XSH_NIR' in self.instruments):
            # xsh_nir_kernel  = Gaussian1DKernel(stddev=4.2/2.355)
            xsh_nir_kernel  = Gaussian1DKernel(stddev=2.2/2.355)
            specobj_xsh_nir = vm.Spectrum(self.spec['XSH_NIR']['wave'],self.spec['XSH_NIR']['flux']/self.spec['XSH_NIR']['cont'], \
                                          1/np.sqrt(self.spec['XSH_NIR']['ivar'])/self.spec['XSH_NIR']['cont'],xsh_nir_kernel,\
                                          lines=linelist)
        if ('XSH_VIS' in self.instruments):
            xsh_vis_kernel  = Gaussian1DKernel(stddev=4.8/2.355)
            specobj_xsh_vis = vm.Spectrum(self.spec['XSH_VIS']['wave'],self.spec['XSH_VIS']['flux']/self.spec['XSH_VIS']['cont'], \
                                          1/np.sqrt(self.spec['XSH_VIS']['ivar'])/self.spec['XSH_VIS']['cont'],xsh_vis_kernel,\
                                          lines=linelist)
        
        
        if (self.vpTree != None):

            for i in range(3):
                if (self.plotinstruments[i] == 'XSH_NIR'):
                    thisprof = vf.vpTau2Flux(vf.vpFromModel(self.vpTree.vp_model, specobj_xsh_nir),xsh_nir_kernel)
                    thiswave = self.spec['XSH_NIR']['wave']
                elif(self.plotinstruments[i] == 'XSH_VIS'):
                    thisprof = vf.vpTau2Flux(vf.vpFromModel(self.vpTree.vp_model, specobj_xsh_vis),xsh_vis_kernel)
                    thiswave = self.spec['XSH_VIS']['wave']
                elif(self.plotinstruments[i] == 'FIRE'):
                    thisprof = vf.vpTau2Flux(vf.vpFromModel(self.vpTree.vp_model, specobj_fire),fire_kernel)
                    thiswave = self.spec['FIRE']['wave']
                elif(self.plotinstruments[i] == 'HIRES'):
                    thisprof = vf.vpTau2Flux(vf.vpFromModel(self.vpTree.vp_model, specobj_hires),hires_kernel)
                    thiswave = self.spec['HIRES']['wave']
                self.ax[i].plot(thiswave,thisprof,color='c',alpha=1.0)

    def plotVPfits(self):

        if (self.profs_fire == None):
            return()
        
        if (self.plotinstruments[0] == 'FIRE'):
            profs0 = self.profs_fire
        elif (self.plotinstruments[0] == 'HIRES'):
            profs0 = self.profs_hires
        elif (self.plotinstruments[0] == 'XSH_VIS'):
            profs0 = self.profs_xsh_vis
        elif (self.plotinstruments[0] == 'XSH_NIR'):
            profs0 = self.profs_xsh_nir

        for thisprof in profs0:
            self.ax[0].plot(self.spec[self.plotinstruments[0]]['wave'],thisprof,color='r',alpha=0.2)

        if (self.plotinstruments[1] == 'FIRE'):
            profs1 = self.profs_fire
        elif (self.plotinstruments[1] == 'HIRES'):
            profs1 = self.profs_hires
        elif (self.plotinstruments[1] == 'XSH_VIS'):
            profs1 = self.profs_xsh_vis
        elif (self.plotinstruments[1] == 'XSH_NIR'):
            profs1 = self.profs_xsh_nir

        for thisprof in profs1:
            self.ax[1].plot(self.spec[self.plotinstruments[1]]['wave'],thisprof,color='r',alpha=0.2)


        if (self.plotinstruments[2] == 'FIRE'):
            profs2 = self.profs_fire
        elif (self.plotinstruments[2] == 'HIRES'):
            profs2 = self.profs_hires
        elif (self.plotinstruments[2] == 'XSH_VIS'):
            profs2 = self.profs_xsh_vis
        elif (self.plotinstruments[2] == 'XSH_NIR'):
            profs2 = self.profs_xsh_nir

        for thisprof in profs2:
            self.ax[2].plot(self.spec[self.plotinstruments[2]]['wave'],thisprof,color='r',alpha=0.2)

############################################################################            
               
class SpecSelect(Ui_SpectrumSelector):

    def __init__(self, selector):
        Ui_SpectrumSelector.__init__(self)
        self.setupUi(selector)
        self.buttonBox.accepted.connect(self.setplots)
        self.buttonBox.rejected.connect(self.reject)
        self.yesno = None
        
    def setplots(self):
        self.plotlist = [self.spec1ComboBox.currentText()]
        self.plotlist.append(self.spec2ComboBox.currentText())
        self.plotlist.append(self.spec3ComboBox.currentText())
        self.yesno = True
            
    def reject(self):
        print("Spectrum selection cancelled")
        self.yesno = False

############################################################################
        
class LineSelect(Ui_LineSelector):

    def __init__(self, selector):
        Ui_LineSelector.__init__(self)
        self.setupUi(selector)
        self.buttonBox.accepted.connect(self.chooseline)
        self.buttonBox.rejected.connect(self.devnull)
        self.linewave = 0
        self.yesno = None
        
    def chooseline(self):
        selected_line = self.lineListWidget.selectedItems()[0]
        print(selected_line.text().split()[1])
        self.linewave = float(selected_line.text().split()[1])
        self.yesno = True

    def devnull(self):
        print("Line selection cancelled")
        self.yesno = False

############################################################################        

class VPModelTree(Ui_VoigtProfileModel):

    def __init__(self, selector, parentobj):
        Ui_VoigtProfileModel.__init__(self)
        self.setupUi(selector)
        self.parentobj = parentobj

        try:
            print(self.components)
        except:
            self.components = []

        self.vpModelTree.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.vpModelTree.customContextMenuRequested.connect(self.openMenu)
        self.writeButton.clicked.connect(self.traverseVPTree)
        self.readButton.clicked.connect(self.loadVPTree)
        self.clearButton.clicked.connect(self.clearVPTree)
        
        try:
            self.vpModelTree.setModel(self.treeModel)
        except:
            self.treeModel = QStandardItemModel()
            self.rootNode = self.treeModel.invisibleRootItem()
            self.vpModelTree.setModel(self.treeModel)

        self.treeModel.setColumnCount(3)
        self.vpModelTree.setColumnWidth(0,300)
        self.vpModelTree.setColumnWidth(1,100)
        self.vpModelTree.setColumnWidth(2,100)
        
        self.treeModel.setHorizontalHeaderLabels(['Component','b (km/s)','N (cm-2)'])
        self.vp_model = vm.Model()        


    def updateModel(self,item):
        # print("Updating the model")
        txt = item.text()
        dat = item.data()
        row = item.row()
        col = item.column()
        if (col == 0):
            # This is changing the redshift of a component, need to set the data field
            item.setData([txt,dat[1]])
            dat2=''
        elif (col == 1):
            # This is changing the b of a component, need to set the data field
            parentitem = self.treeModel.item(row,column=0)
            dat2 = parentitem.data()
            parentitem.setData([dat2[0],float(txt)])
        elif (col == 2):
            # This is changing the column density of an ion
            parentrow = item.parent().row()
            parentitem = self.treeModel.item(parentrow,column=0)
            dat2 = parentitem.data()
            parentitem.child(row).setData([dat2[0],float(txt)])
            
        self.makeVPModel()
        GuiProgram.change_plot(self.parentobj)
        
    def add_component(self, redshift, b = 10.0):
        newcomponent = vpComponent(redshift, bparam=b)
        self.components.append(newcomponent)
        self.rootNode.appendRow(newcomponent)
        rr = newcomponent.row()
        mm = newcomponent.model()
        bparam_item = QStandardItem(str(b))
        mm.setItem(rr,1,bparam_item)

    def add_ion(self, redshift, ion, column=13.5):
        newion = vpIon(ion, column=column)

        component_redshifts = [c.redshift for c in self.components]
        dz = np.abs(np.array(component_redshifts) - redshift)
        mindz = np.min(dz)
        if (mindz > 0.005):
            print("No matching absorption component found within dz = 0.005")
        else:
            for c in self.components:
                if (abs(redshift-c.redshift) == mindz):
                    logN_item = QStandardItem(str(column))
                    blank_item = QStandardItem('')
                    c.appendRow([newion,blank_item,logN_item])
                    
    def openMenu(self, position):
        indexes = self.vpModelTree.selectedIndexes()
        if len(indexes) > 0:
            level = 0
            index = indexes[0]
            while index.parent().isValid():
                index = index.parent()
                level += 1

        self.menu = QtWidgets.QMenu()

        if (level == 0):
            tt = Table.read('../atomic_data.txt',format='ascii')
            ion_list = unique(tt,keys='ion')['ion']
            [self.menu.addAction(ion_list[i]) for i in range(len(ion_list))]
            self.menu.addSeparator()
            self.menu.addAction("Delete Component")
            self.menu.triggered.connect(self.addIon)
        elif (level == 1):
            ion = indexes[0].data()
            linelist = Table.read('../atomic_data.txt',format='ascii')
            transitions = linelist[linelist['ion'] == ion]['wave']
            for t in transitions:
                action = self.menu.addAction(f"{np.floor(t):6.1f}")
            self.menu.addSeparator()
            self.menu.addAction("Delete Ion")
            self.menu.triggered.connect(self.addTransition)
        elif (level == 2):
            self.menu.addAction("Delete Transition")
            self.menu.triggered.connect(self.removeTransition)
        elif (level == 3):
            self.menu.addAction("Delete FitRegion")
            self.menu.triggered.connect(self.removeFitRegion)
            
        self.menu.exec_(self.vpModelTree.viewport().mapToGlobal(position))

        
    def addIon(self,action):
        selected_items = self.vpModelTree.selectionModel().selectedRows()
        if (action.text() == "Delete Component"):
            if (len(selected_items) != 0):
                self.treeModel.removeRow(selected_items[0].row())
            return()
        else:
            newion = vpIon(action.text(),column=12.0)
            logN_item = QStandardItem('12.0')
            blank_item = QStandardItem('')
            indexes = self.vpModelTree.selectedIndexes()
            item = self.treeModel.itemFromIndex(indexes[0])
            item.appendRow([newion,blank_item,logN_item])
                
    def removeIon(self,action):
        if (action.text() == "Delete Ion"):
            selected_items = self.vpModelTree.selectionModel().selectedRows()
            indexes = self.vpModelTree.selectedIndexes()
            for ii in sorted(indexes):
                selected_item = self.treeModel.itemFromIndex(ii)
                parent = selected_item.parent()
                print(selected_item)
                print(f"Deleting: {selected_item.data()}, {selected_item.row()} ")
                parent.removeRow(selected_item.row())
                
    def addTransition(self,action):

        if (action.text() == "Delete Ion"):
            self.removeIon(action)
            return()
        else:
            restwv = float(action.text())

        indexes = self.vpModelTree.selectedIndexes()
        redshift = float(indexes[0].parent().data())
        ion = indexes[0].data()
        newtransition = vpTransition(restwv)
        self.treeModel.itemFromIndex(indexes[0]).appendRow(newtransition)

    def removeTransition(self,action):
        if (action.text() == "Delete Transition"):
            selected_items = self.vpModelTree.selectionModel().selectedRows()
            indexes = self.vpModelTree.selectedIndexes()
            for ii in sorted(indexes):
                selected_item = self.treeModel.itemFromIndex(ii)
                parent = selected_item.parent()
                print(selected_item)
                print(f"Deleting: {selected_item.data()}, {selected_item.row()} ")
                parent.removeRow(selected_item.row())
            
    def addFitRegion(self,instrument,minwv,maxwv):

        indexes = self.vpModelTree.selectedIndexes()
        if len(indexes) > 0:
            level = 0
            index = indexes[0]
            while index.parent().isValid():
                index = index.parent()
                level += 1

        if (level < 2):
            print("Please select a transition (by restwv) to which to append the fit region")
            return()
        else:
            highlighted = self.vpModelTree.selectedIndexes()
            
        if (len(highlighted) == 1):
            thisitem = self.treeModel.itemFromIndex(highlighted[0])
            if (level == 2):

                restwv_approx = float(thisitem.data(0))
                tt = abs(self.parentobj.linelist['wave']-restwv_approx)
                restwv = float(self.parentobj.linelist[tt == min(tt)]['wave'])

                redshift = float(thisitem.parent().parent().data(0))
                obswv = restwv * (1+redshift)
                vmin = (minwv - obswv)/obswv * 299792.4
                vmax = (maxwv - obswv)/obswv * 299792.4
                newregion = vpFitRegion(instrument,[vmin,vmax])
                self.treeModel.itemFromIndex(highlighted[0]).appendRow(newregion)
                print(f"{redshift},{minwv},{maxwv},{restwv}")
            elif (level == 3):
                restwv_approx   = float(thisitem.parent().data(0))
                tt = abs(self.parentobj.linelist['wave']-restwv_approx)
                restwv = float(self.parentobj.linelist[tt == min(tt)]['wave'])
                
                redshift = float(thisitem.parent().parent().parent().data(0))
                obswv = restwv * (1+redshift)
                vmin = (minwv - obswv)/obswv * 299792.4
                vmax = (maxwv - obswv)/obswv * 299792.4
                newregion = vpFitRegion(instrument,[vmin,vmax])
                self.treeModel.itemFromIndex(highlighted[0]).parent().appendRow(newregion)
            else:
                print("Please select a transition (by restwv) to which to append the fit region")
                return()
        else:
            print("Error: only one transition can be assigned a fit region at one time")

    def removeFitRegion(self, action):
        if (action.text() == "Delete FitRegion"):
            selected_items = self.vpModelTree.selectionModel().selectedRows()
            indexes = self.vpModelTree.selectedIndexes()
            for ii in sorted(indexes):
                selected_item = self.treeModel.itemFromIndex(ii)
                parent = selected_item.parent()
                print(f"Deleting: {selected_item.data()} ")
                parent.removeRow(selected_item.row())

    def clearVPTree(self):
        # Clear any existing tree
        self.treeModel.removeRows(0,self.treeModel.rowCount())
        self.vp_model = None
        self.vp_model = vm.Model()
        self.components = []
        
    def loadVPTree(self):                
        # Clear any existing tree
        try:
            self.treeModel.removeRows(0,self.treeModel.rowCount())
        except:
            print("Reloading model")
            
        # Get the file name from the user
        options = QFileDialog.Options()
        fileName, _ = QFileDialog.getOpenFileName(self.vpModelTree, \
                                                  "Open","","All Files (*);;Text Files (*.txt)", options=options)

        # Read in the requested filename and populate the tree
        if (fileName):
            with open(fileName, "r") as fp:
                entries = fp.readlines()

            # First loop through and populate all the components (level 1)
            for l in entries:
                if('addcomponent' in l):
                    redshift = float(l.split(',')[0].split('(')[1])
                    bparam = float(l.split('b_turb=')[1].split(')')[0])
                    self.add_component(redshift, b=bparam)
                    
            # Second, loop through and populate all of the ions
            for l in entries:
                if('addion' in l):
                    ion = l.split('\'')[1]
                    redshift = float(l.split(",")[1])
                    column = float(l.split("N=")[1].split(',')[0])
                    self.add_ion(redshift,ion,column=column)
                    
            # Third, loop through and populate all of the transitions
            for l in entries:
                if('addtransition' in l):
                    restwv = float(l.split(',')[0].split('(')[1])
                    ion = l.split("\'")[1]
                    redshift = float(l.split(',')[-1][:-2])
                    for i in range(self.treeModel.rowCount()):
                        indx = self.treeModel.index(i,0)
                        component_z = float(self.treeModel.data(indx))
                        # print(f"component_z: {component_z}, {redshift}")
                        if (redshift == component_z):
                            component = self.treeModel.itemFromIndex(indx)
                            break

                    if component.hasChildren():
                        # print(f"Adding transition: {restwv}, {component_z}")
                        for i in range(component.rowCount()):
                            if(component.child(i).data()[0] == ion):
                                component.child(i).appendRow(vpTransition(restwv))
                                
            # Finally, loop through and populate all of the fitregions
            for l in entries:
                if('fitregion' in l):
                    restwv = float(l.split(',')[0].split('(')[1])
                    ion_name = l.split("\'")[1]
                    redshift = float(l.split(',')[2])
                    instrument = l.split("\'")[-2]
                    lower_dv = float(l.split('[')[1].split(',')[0])
                    upper_dv = float(l.split(']')[0].split(',')[-1])


                    for i in range(self.treeModel.rowCount()):
                        indx = self.treeModel.index(i,0)
                        component_z = float(self.treeModel.data(indx))
                        if (redshift == component_z):
                            component = self.treeModel.itemFromIndex(indx)
                            break

                    if component.hasChildren():
                        for i in range(component.rowCount()):
                            if(component.child(i).data(0) == ion_name):
                                ion = component.child(i)
                                
                    if ion.hasChildren():
                        for i in range(ion.rowCount()):
                            if (ion.child(i).data() == restwv):
                                ion.child(i).appendRow(vpFitRegion(instrument,[lower_dv,upper_dv]))

        self.makeVPModel()
        self.treeModel.itemChanged.connect(self.updateModel)
        try:
            GuiProgram.change_plot(self.parentobj)
        except:
            print("WARNING: loading model, but no spectrum has been read in")

    # This is to write out a file for running fitting
    def traverseVPTree(self):

        options = QFileDialog.Options()
        fileName, _ = QFileDialog.getSaveFileName(self.vpModelTree, "Save","","All Files (*);;Text Files (*.txt)", options=options)
        commands = []

        if (fileName == None or fileName == ''):
            return()
        
        ncomponents = 0
        
        for i in range(self.treeModel.rowCount()):
            item = self.treeModel.item(i)
            commands.append(f'm.addcomponent({item.data()[0]},bpriors=[3,50],b_turb={item.data()[1]})')
            ncomponents += 1
            level = 0
            self.getItem(item,level,commands)

        print(f"Ncomponents = {ncomponents}")
        with open(fileName,'w') as fp:
            for cmd in commands:
                print(cmd)
                fp.write(f'{cmd}\n')
            
    # For recursive descent through the tree.
    def getItem(self, item, level, commands):
        if (item != None):
            nions = 0
            if item.hasChildren():
                level += 1
                for i in range(item.rowCount()):
                    childitem = item.child(i)
                    if (childitem != None):
                        if (level == 1):
                            commands.append(f'm.addion(\'{childitem.data()[0]}\',{item.data(0)},N={childitem.data()[1]},Npriors=[11,15])')
                            nions += 1
                        elif (level == 2):
                            commands.append(f'm.addtransition({childitem.data(0)},\'{item.data(0)}\',{item.parent().data(0)})')
                        elif (level == 3):
                            # Restwv, ion, redshift [bounds], instrument
                            # commands.append(f'm.fitregion({childitem.data(0)},\'{item.data(0)}\',{item.parent().parent().data(0)})')
                            commands.append(f'm.addfitregion({item.data(0)},\'{item.parent().data(0)}\',{item.parent().parent().data(0)},[{childitem.data()[1]:4.1f},{childitem.data()[2]:4.1f}],\'{childitem.data()[0]}\')')
                            
                    self.getItem(childitem, level, commands)

                print(f"Nions = {nions}")
                return(commands)

    def makeVPModel(self):

        # clear the model
        # print("Clearing the model")
        self.vp_model = None
        self.vp_model = vm.Model()

        for i in range(self.treeModel.rowCount()):
            item = self.treeModel.item(i)
            z = float(item.data()[0])
            b = float(item.data()[1])
            # print(f"AddComponent: {z},{b}")
            self.vp_model.addcomponent(z,b_turb=b,bpriors=[5,100])
            level = 0
            self.descendTree(item,level)

    def descendTree(self, item, level):
        if (item != None):
            nions = 0
            if item.hasChildren():
                level += 1
                for i in range(item.rowCount()):
                    childitem = item.child(i)
                    if (childitem != None):
                        if (level == 1):
                            ion_name = childitem.data(0)
                            z = float(item.data()[0])
                            N = childitem.data()[1]
                            self.vp_model.addion(ion_name, z, N=N)
                            # print(f"AddIon: {ion_name},{z},{N}")
                        elif (level == 2):
                            restwv = float(childitem.data())
                            ion_name = item.data()[0]
                            z = float(item.parent().data()[0])
                            # print(f"Adding transition {restwv},{ion_name},{z}")
                            self.vp_model.addtransition(restwv, ion_name, z)
                        elif (level == 3):
                            # Restwv, ion, redshift [bounds], instrument
                            self.vp_model.addfitregion(float(item.data(0)),item.parent().data(0),float(item.parent().parent().data(0)),[float(childitem.data()[1]),float(childitem.data()[2])],childitem.data()[0])
                            
                    self.descendTree(childitem, level)


############################################################################
                    
class vpComponent(QStandardItem):

    def __init__(self,redshift,bparam=10.0):
        super().__init__()
        self.redshift = redshift

        self.setText(str(redshift))
        self.setData([redshift,bparam])
        
class vpIon(QStandardItem):

    def __init__(self,ion,column=14.0):
        super().__init__()
        self.ion = ion
        self.setText(ion)
        self.setData([ion,column])

class vpTransition(QStandardItem):

    def __init__(self,transition):
        super().__init__()
        self.restwv = float(transition)
        self.setText(f"{transition:6.1f}")
        self.setData(transition)

class vpFitRegion(QStandardItem):

    def __init__(self,instrument,fitbounds):
        super().__init__()
        self.fitbounds = fitbounds
        self.setText(f"{instrument}:\t[{fitbounds[0]:5.2f},{fitbounds[1]:5.2f}] km/s")
        tmp = QtCore.QVariant([instrument,fitbounds[0],fitbounds[1]])
        self.setData(tmp)
        
        
