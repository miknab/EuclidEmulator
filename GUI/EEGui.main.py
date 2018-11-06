import sys
from PyQt5 import QtGui
from PyQt5.QtWidgets import QApplication, QLabel, QLineEdit, QSlider, QTextEdit, QWidget, QPushButton, QVBoxLayout, \
    QHBoxLayout, QFileDialog, QGridLayout
from PyQt5.QtCore import Qt

import CosmoPanel as CP
import PlotPanel as PP

class Window(QWidget):
    def __init__(self):
            # super().__init() is the python3-specific syntax
            super(Window, self).__init__()  # this syntax is compatible with python2 & python3
            self.init_ui()

    def init_ui(self):
        # Emulator parameter space
        parList = [("omega_b",[0.0217,0.0233]),
                   ("omega_m",[0.1326,0.1526]),
                   ("n_s",[0.9345,0.9965]),
                   ("h", [0.6251,0.7211]),
                   ("w_0", [-1.25,-0.75]),
                   ("sigma_8", [0.7684,0.8614]),
                   ( "z", [0.0,5.0])]

        # Take EuclidReference cosmology as default
        DefaultCosmoDict = {'omega_b': 0.021996,
                            'omega_m': 0.121203,
                            'n_s': 0.96,
                            'h': 0.67,
                            'w_0': -1.0,
                            'sigma_8': 0.83}

        DefaultRedshift = 0.0

        # Create grid layout
        glay = QGridLayout()

        # Create canvas
        canv = PP.PlotCanvas(DefaultCosmoDict, DefaultRedshift)
        glay.addWidget(canv)

        # Create parameter panel
#        glay.addWidget(QLabel("Cosmology"), 0,0,1,6)
#
#        for idx, par in enumerate(parList): 
#            CP.CreateCosmoField(glay, par, idx+1, canv, DefaultCosmoDict, DefaultRedshift)

        self.setWindowTitle("EuclidEmulator")
        self.setLayout(glay)

        self.show()
     
app = QApplication(sys.argv)  # create application loop
a_window = Window()
sys.exit(app.exec_())
