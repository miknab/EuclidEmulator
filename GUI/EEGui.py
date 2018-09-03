# import numpy as np
import sys
from PyQt5 import QtGui
from PyQt5.QtWidgets import QApplication, QLabel, QLineEdit, QSlider, QTextEdit, QWidget, QPushButton, QVBoxLayout, \
    QHBoxLayout, QFileDialog
from PyQt5.QtCore import Qt


class Window(QWidget):

    def __init__(self):
        # super().__init() is the python3-specific syntax
        super(Window, self).__init__()  # this syntax is compatible with python2 & python3
        self.init_ui()

    def init_ui(self):
        self.btn = QPushButton('Push Me!')

        self.lbl = QLabel()
        self.lbl.setPixmap(QtGui.QPixmap('phi.jpg'))

        top_layout = QHBoxLayout()

        # Cosmology layout
        cosmo_layout = QVBoxLayout()
        self.ttllabel = QLabel("Cosmology:")
        cosmo_layout.addWidget(self.ttllabel)
        cosmo_layout.addLayout(self.cosmo_par_layout("Omega_b"))
        cosmo_layout.addLayout(self.cosmo_par_layout("Omega_0"))
        cosmo_layout.addLayout(self.cosmo_par_layout("n_s"))
        cosmo_layout.addLayout(self.cosmo_par_layout("h"))
        cosmo_layout.addLayout(self.cosmo_par_layout("w0"))
        cosmo_layout.addLayout(self.cosmo_par_layout("sigma_8"))
        cosmo_layout.addLayout(self.cosmo_par_layout("z"))

        # Plot layout
        plot_layout = QVBoxLayout()
        plot_layout.addWidget(self.lbl)

        top_layout.addLayout(cosmo_layout)
        top_layout.addLayout(plot_layout)

        bottom_layout = QHBoxLayout()
        config_layout = QVBoxLayout()
        config_layout.addWidget(self.lbl)
        menubtn_layout = QHBoxLayout()
        menubtn_layout.addWidget(self.btn)

        bottom_layout.addLayout(config_layout)
        bottom_layout.addLayout(menubtn_layout)

        main_layout = QVBoxLayout()
        main_layout.addLayout(top_layout)
        main_layout.addLayout(bottom_layout)

        self.setWindowTitle("EuclidEmulator")
        self.setLayout(main_layout)

        self.show()

    def cosmo_par_layout(self, ttl):
        lyt = QHBoxLayout()
        self.lab = QLabel(ttl)
        self.le = QLineEdit()
        self.sld = QSlider(Qt.Horizontal)
        self.sld.setMinimum(0.0)
        self.sld.setMaximum(1.0)
        self.sld.setValue(0.5)
        self.sld.setTickInterval(0.1)
        self.sld.setTickPosition(QSlider.TicksBelow)

        lyt.addWidget(self.lab)
        lyt.addWidget(self.le)
        lyt.addWidget(self.sld)

        return lyt


app = QApplication(sys.argv)  # create application loop
a_window = Window()
sys.exit(app.exec_())
