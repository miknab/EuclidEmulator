import sys
from PyQt5 import QtGui
from PyQt5.QtWidgets import QApplication, QLabel, QLineEdit, QSlider, QTextEdit, QWidget, QPushButton, QVBoxLayout, \
    QHBoxLayout, QFileDialog, QGridLayout
from PyQt5.QtCore import Qt

import PlotPanel as PP

def CreateCosmoField(layout, (label, parrange), row, canvas):

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

    layout.addWidget(qlab, row, 0, 2, 1)
    layout.addWidget(qle, row, 2, 2, 1)
    layout.addWidget(qsld, row, 4, 2, 1)

    qsld.valueChanged.connect(lambda: v_change(qsld, qle, parrange, maximum, layout, canvas))

def v_change(sld, le, parrange, maximum, layout, canvas):
        my_value = str(linparmap(sld.value(),parrange, maximum))
        le.setText(my_value)

        for sld in layout.findChildren(QSlider):
            print(sld.value())

        # Need to get CosmoDict and z
        #canvas.plot(CosmoDict, z)

def linparmap(val,parrange,maximum):
        del_old = maximum  # del_old = value range of slider (in Qt this are
                           # only integers). We start from 0 such that
                           # del_old = maximum. 
        del_new = parrange[1]-parrange[0]
        slope = del_new/del_old

        return slope*val+parrange[0]


