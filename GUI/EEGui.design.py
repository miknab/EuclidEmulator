import sys
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

app = QApplication(sys.argv)
w = QWidget()
glay = QGridLayout(w)
glay.addWidget(QLabel("Cosmology"), 0, 0, 1, 6)
glay.addWidget(QLabel("Plot"), 0, 9, 10, 6)
glay.addWidget(QLabel("Omb_lb"), 1, 0, 1, 2)
glay.addWidget(QLabel("Omb_le"), 1, 2, 1, 2)
glay.addWidget(QLabel("Omb_sld"), 1, 4, 1, 2)
glay.addWidget(QLabel("Om0_lb"), 2, 0, 1, 2)
glay.addWidget(QLabel("Om0_le"), 2, 2, 1, 2)
glay.addWidget(QLabel("Om0_sld"), 2, 4, 1, 2)
glay.addWidget(QLabel("ns_lb"), 3, 0, 1, 2)
glay.addWidget(QLabel("ns_le"), 3, 2, 1, 2)
glay.addWidget(QLabel("ns_sld"), 3, 4, 1, 2)
glay.addWidget(QLabel("h_lb"), 4, 0, 1, 2)
glay.addWidget(QLabel("h_le"), 4, 2, 1, 2)
glay.addWidget(QLabel("h_sld"), 4, 4, 1, 2)
glay.addWidget(QLabel("w0_lb"), 5, 0, 1, 2)
glay.addWidget(QLabel("w0_le"), 5, 2, 1, 2)
glay.addWidget(QLabel("w0_sld"), 5, 4, 1, 2)
glay.addWidget(QLabel("sig8_lb"), 6, 0, 1, 2)
glay.addWidget(QLabel("sig8_le"), 6, 2, 1, 2)
glay.addWidget(QLabel("sig8_sld"), 6, 4, 1, 2)
glay.addWidget(QLabel("z_lb"), 7, 0, 1, 2)
glay.addWidget(QLabel("z_le"), 7, 2, 1, 2)
glay.addWidget(QLabel("z_sld"), 7, 4, 1, 2)

glay.addWidget(QLabel("Configuration"), 9, 0, 1, 6)
glay.addWidget(QLabel("Plin(k)_chkbx"), 10, 0, 1, 3)
glay.addWidget(QLabel("Pnonlin(k)_chkbx"), 11, 0, 1, 3)
glay.addWidget(QLabel("Boost(k)_chkbx"), 12, 0, 1, 3)
glay.addWidget(QLabel("config4"), 10, 3, 1, 3)
glay.addWidget(QLabel("config5"), 11, 3, 1, 3)
glay.addWidget(QLabel("config6"), 12, 3, 1, 3)

glay.addWidget(QLabel("Save data button"), 11, 10, 1, 2)
glay.addWidget(QLabel("Save plot button"), 11, 12, 1, 2)
qsrand(QTime.currentTime().msec())

for label in w.findChildren(QLabel):
    color = QColor(qrand() % 256, qrand() % 256, qrand() % 256)
    label.setStyleSheet('.QLabel{{background: rgb({}, {}, {});}}'.format(color.red(), color.green(), color.blue()))

w.show()
sys.exit(app.exec_())
