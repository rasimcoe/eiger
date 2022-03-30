from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.Qt import QStandardItemModel, QStandardItem
from gui import Ui_Dialog
from SpecGui import Ui_SpectrumSelector
from LineGui import Ui_LineSelector
from vpModelGui import Ui_VoigtProfileModel
from eiger.QuasarSpec.loadQsoSpec import loadQsoSpec
from astropy.table import Table
import matplotlib.pyplot as plt
from pypeit.core.wave import airtovac
import astropy.units as u
import numpy as np

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

        ''' This method gets called when the window is created. '''
        Ui_Dialog.__init__(self)              # Initialize Window
        self.setupUi(dialog)                  # Set up the UI
        # Initialize the figure in our window
        figure = Figure()                     # Prep empty figure
#        axis = figure.add_subplot(111)        # Prep empty plot
        axis = figure.subplots(3,1,sharex=True,gridspec_kw={'hspace':0})        # Prep empty plot
        self.initialize_figure(figure, axis)  # Initialize!
        # Connect our button with plotting function
        self.reloadPushButton.clicked.connect(self.load_newplot)
        self.SpectraSelectionPushButton.clicked.connect(self.choose_spectra)
        self.vpTreePushButton.clicked.connect(self.show_vptree)
        self.canvas.mpl_connect('key_press_event',self.on_keypress)
        self.canvas.mpl_connect('button_release_event',self.on_buttonpress)
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

                if (instrument == 'HIRES'):
                    xdata = np.array(airtovac(self.spec[instrument]['wave']*u.AA))
                else:
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

        self.spec = loadQsoSpec(1,revision='current')
        self.instruments = list(self.spec.keys())[4:]

        tmp = loadQsoSpec(1,revision='1')
        self.spec['FIRE'] = tmp['FIRE']
        
        self.xmin = min(self.spec['FIRE']['wave'])
        self.xmax = max(self.spec['FIRE']['wave'])
        self.ymin = -1.0
        self.ymax = 3.0
        self.change_plot()

    def on_buttonpress(self,event):

        # On right click, open line ID GUI
        if (event.button == 3):
            wave = event.xdata
            self.set_abslineid(wave)
        
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
            self.idtable.write("J1030_idtable.dat",format='ascii.fixed_width')
            print("Writing out ASCII table")
            redraw = False
        elif (event.key == 'R'):
            self.idtable = Table.read("J1030_idtable.dat",format='ascii.fixed_width')
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
                        self.ax[i].plot([obswave,obswave],[1,1.4],alpha=0.5,color='k')
                        self.ax[i].text(obswave,1.6,ion,rotation='vertical',horizontalalignment='center',\
                                        picker=True, fontsize=8)
                    
        for line in self.idtable:
            obswave = line['restwv'] * (1+line['redshift'])
            ion = f"{line['ion']} {line['restwv']:4.0f} ({line['redshift']:4.3f})"
            if (obswave > self.xmin and obswave < self.xmax):
                nplots = len(self.ax)
                for i in range(nplots):
                    self.ax[i].plot([obswave,obswave],[1,1.4],alpha=0.5,color='r')
                    self.ax[i].text(obswave,1.6,ion,rotation='vertical',horizontalalignment='center',\
                                    picker=True, color='r', fontsize=8)
                    
                    
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

        self.selector.show()
        self.lineDialog.buttonBox.clicked.connect(lambda: self.register_z(wave))

    def register_z(self,wave):
        self.zabs = float(wave) / float(self.lineDialog.linewave) - 1.0
        print(f"Redshift: {self.zabs}")
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
            self.vpTree = VPModelTree(self.tree)
            self.tree.show()
        
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


class VPModelTree(Ui_VoigtProfileModel):

    def __init__(self, selector):
        Ui_VoigtProfileModel.__init__(self)
        self.setupUi(selector)

        try:
            print(self.components)
        except:
            self.components = []

        try:
           self.vpModelTree.setModel(self.treeModel)
        except:
            self.treeModel = QStandardItemModel()
            self.rootNode = self.treeModel.invisibleRootItem()

            self.vpModelTree.setModel(self.treeModel)

    def add_component(self, redshift):
        newcomponent = vpComponent(redshift)
        self.components.append(newcomponent)
        self.rootNode.appendRow(newcomponent)

    def add_ion(self, redshift, ion):
        newion = vpIon(ion)
        for component in self.components:
            print(f"{component.data()},{redshift}")
            if (component.data() == redshift):
                component.appendRow(newion)

    def add_transition(self,redshift,ion,restwv):
        newtransition = vpTransition(restwv)
        
                
class vpComponent(QStandardItem):

    def __init__(self,redshift):
        super().__init__()
        self.redshift = redshift

        self.setText(str(redshift))
        self.setData(redshift)
        
class vpIon(QStandardItem):

    def __init__(self,ion):
        super().__init__()

        self.setText(ion)
        self.setData(ion)

class vpTransition(QStandardItem):

    def __init__(self,transition):
        super().__init__()

        self.setText(str(transition))
        self.setData(transition)

        
