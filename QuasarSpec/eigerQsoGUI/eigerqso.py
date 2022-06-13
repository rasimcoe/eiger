#! /usr/bin/env python

import sys
from PyQt5 import QtCore, QtGui, QtWidgets
from gui import Ui_Dialog
from SpecGui import Ui_SpectrumSelector

from gui_logic import GuiProgram, SpecSelect

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    dialog = QtWidgets.QDialog()
    program = GuiProgram(dialog)
    dialog.show()
    sys.exit(app.exec_())

