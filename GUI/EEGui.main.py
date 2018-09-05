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
        glay.addWidget(QLabel("Cosmology"), 0,0,1,6)
        parList = [("omega_b",[0.0217,0.0233]),
                   ("omega_m",[0.1326,0.1526]),
                   ("n_s",[0.9345,0.9965]),
                   ("h", [0.6251,0.7211]),
                   ("w_0", [-1.25,-0.75]),
                   ("sigma_8", [0.7684,0.8614]),
                   ( "z", [0.0,5.0])]

        for idx, par in enumerate(parList): 
            CreateCosmoField(glay, par, idx+1)

        self.setWindowTitle("EuclidEmulator")
        self.setLayout(glay)

        self.show()
    
def CreateCosmoField(layout, (label, parrange), row):

    minimum = 0
    maximum = 1000

    if label != "z":
        startvalue = (parrange[1]+parrange[0])/2
        start = int((maximum-minimum)/2)
    else:
        startvalue = 0
        start = minimum

    qlab = QLabel(label)
    qle = QLineEdit(str(startvalue))
    qsld = QSlider(Qt.Horizontal)
    qsld.setMinimum(minimum)
    qsld.setMaximum(maximum)
    qsld.setValue(start)
    qsld.setTickInterval(int(maximum/10))
    qsld.setTickPosition(QSlider.TicksBelow)

    qsld.valueChanged.connect(lambda: v_change(qsld, qle, parrange, maximum))

    layout.addWidget(qlab, row, 0, 2, 1)
    layout.addWidget(qle, row, 2, 2, 1)
    layout.addWidget(qsld, row, 4, 2, 1)

def v_change(sld, le, parrange, maximum):
        my_value = str(linparmap(sld.value(),parrange, maximum))
        le.setText(my_value)

def linparmap(val,parrange,maximum):
        del_old = maximum  # del_old = value range of slider (in Qt this are
                           # only integers). We start from 0 such that
                           # del_old = maximum. 
        del_new = parrange[1]-parrange[0]
        slope = del_new/del_old
        
        return slope*val+parrange[0]
        

app = QApplication(sys.argv)  # create application loop
a_window = Window()
sys.exit(app.exec_())
