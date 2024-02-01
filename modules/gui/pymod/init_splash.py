from PyQt5 import QtCore, QtGui, QtWidgets

import os
import ost
from ost import gui

LOGO_PATH = os.path.join(ost.GetSharedDataPath(), "gui", "images", "logo-small.png")


SPLASH_TEXT=""""Welcome to <b>Openstructure</b>!<br/><br/>
You are running version %s<br /><br />If you are new to OpenStructure, we 
invite you to run the demos from the examples directory. Scripts can be 
displayed by right clicking on the file and selecting 'Show source'.<br/><br/>
Feel free visit our website at:<br /> 
<a href='https://www.openstructure.org'>https://www.openstructure.org</a>
""" % ost.VERSION

class SplashDialog(QtWidgets.QDialog):
  def __init__(self, parent=None):
    QtWidgets.QDialog.__init__(self, parent)
    layout = QtWidgets.QHBoxLayout(self)
    self.setLayout(layout)
    imageLabel = QtWidgets.QLabel();
    self.pix_map = QtGui.QPixmap(LOGO_PATH);
    imageLabel.setPixmap(self.pix_map);
    layout.addWidget(imageLabel)
    self.label = QtWidgets.QTextBrowser()
    self.label.setReadOnly(True)
    self.label.setOpenExternalLinks(True)
    self.label.setHtml(SPLASH_TEXT)
    layout.addWidget(self.label)
    
def _InitSplash():
  splash = SplashDialog(gui.GostyApp.Instance().perspective.main_area.qobject)
  splash.exec_()

  #QtCore.QTimer.singleShot(30000, splash.close);


__all__=('SplashDialog',)
