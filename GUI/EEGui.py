import sys
from PyQt5 import QtGui
from PyQt5.QtWidgets import QApplication, QLabel, QLineEdit, QSlider, QTextEdit, QWidget, QPushButton, QVBoxLayout, \
    QHBoxLayout, QFileDialog, QGridLayout
from PyQt5.QtCore import Qt

class Window(QWidget):
    def __init__(self):
            # super().__init() is the python3-specific syntax
            super(Window, self).__init__()  # this syntax is compatible with python2 & python3
            self.init_ui()

    def init_ui(self):
        glay = QGridLayout()

        parList = ["Omega_b", "Omega_m", "n_s", "h", "w_0", "sigma_8", "z"]

        for idx, label in enumerate(parList): 
            CreateCosmoField(glay, label, idx)

        self.setWindowTitle("EuclidEmulator")
        self.setLayout(glay)

        self.show()
    
def CreateCosmoField(layout, label, row):
    qlab = QLabel(label)
    qle = QLineEdit()
    qsld = QSlider(Qt.Horizontal)
    qsld.setMinimum(0)
    qsld.setMaximum(1000)
    qsld.setValue(500)
    qsld.setTickInterval(100)
    qsld.setTickPosition(QSlider.TicksBelow)

    layout.addWidget(qlab, row, 0, 2, 1)
    layout.addWidget(qle, row, 2, 2, 1)
    layout.addWidget(qsld, row, 4, 2, 1)

app = QApplication(sys.argv)  # create application loop
a_window = Window()
sys.exit(app.exec_())
