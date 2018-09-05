import sys
from PyQt5 import QtGui
from PyQt5.QtWidgets import QApplication, QLabel, QLineEdit, QSlider, QTextEdit, QWidget, QPushButton, QVBoxLayout, \
    QHBoxLayout, QFileDialog, QGridLayout
from PyQt5.QtCore import Qt
import CosmoPanel as CP

class Window(QWidget):
    def __init__(self):
            # super().__init() is the python3-specific syntax
            super(Window, self).__init__()  # this syntax is compatible with python2 & python3
            self.init_ui()

    def init_ui(self):
        glay = QGridLayout()
        glay.addWidget(QLabel("Cosmology"), 0,0,1,6)
        parList = [("omega_b",[0.0217,0.0233]),
                   ("omega_m",[0.1326,0.1526]),
                   ("n_s",[0.9345,0.9965]),
                   ("h", [0.6251,0.7211]),
                   ("w_0", [-1.25,-0.75]),
                   ("sigma_8", [0.7684,0.8614]),
                   ( "z", [0.0,5.0])]

        for idx, par in enumerate(parList): 
            CP.CreateCosmoField(glay, par, idx+1)

        self.setWindowTitle("EuclidEmulator")
        self.setLayout(glay)

        self.show()
     

app = QApplication(sys.argv)  # create application loop
a_window = Window()
sys.exit(app.exec_())
