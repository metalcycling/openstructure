#------------------------------------------------------------------------------
# This file is part of the OpenStructure project <www.openstructure.org>
#
# Copyright (C) 2008-2020 by the OpenStructure authors
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3.0 of the License, or (at your option)
# any later version.
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#------------------------------------------------------------------------------

import sys
from ost import gui
from ost import gfx
import ost

from PyQt5 import QtCore, QtGui
from ost.gui import FileLoader
from ost.gui.init_splash import _InitSplash
from ost.gui.dng import termuse
class InitMenuBar(QtCore.QObject):
  def __init__(self, menu_bar=None):
    QtCore.QObject.__init__(self, menu_bar)
    
    persp=gui.GostyApp.Instance().perspective
    file=persp.GetMenu("File")
    options=persp.GetMenu("Options")
    window=persp.GetMenu("Window")
    help=persp.GetMenu("Help")
    
    load = QtGui.QAction(QtGui.QIcon('icons/open.png'), '&Open', self)
    load.setStatusTip('Load a file')
    load.setShortcut('Ctrl+O')
    self.connect(load, QtCore.SIGNAL('triggered()'), self.Load)
    file.addAction(load)

    webpage = QtGui.QAction('&Documentation', self)
    webpage.setStatusTip('Documentation')
    webpage.setShortcut('Ctrl+D')
    self.connect(webpage, QtCore.SIGNAL('triggered()'), self.OpenDocs)
    help.addAction(webpage)
    if sys.platform=='darwin':
      install_ctl=QtGui.QAction('Install Command Line Tool', self)    
      self.connect(install_ctl, QtCore.SIGNAL('triggered()'), 
                   termuse.InstallTerminalPrograms)
      help.addAction(install_ctl)
    about = QtGui.QAction('&About', self)
    about.setStatusTip('About')
    about.setShortcut('Ctrl+A')
    self.connect(about, QtCore.SIGNAL('triggered()'), self.About)
    help.addAction(about)
    
        
    window.addMenu(persp.panels.menu)
    gl_win = QtGui.QAction('&GL Window', self)
    gl_win.setStatusTip('Display gl windows')
    gl_win.setShortcut('Ctrl+G')
    self.connect(gl_win, QtCore.SIGNAL('triggered()'), self.ShowGLWin)
    window.addAction(gl_win)
    
    reset = QtGui.QAction('Reset View', self)
    reset.setStatusTip('Reset the Panels and Widgets')
    self.connect(reset, QtCore.SIGNAL('triggered()'), self.ResetView)
    window.addAction(reset)
    
  def Exit(self):
    reply = QtGui.QMessageBox()
    reply.addButton(QtGui.QMessageBox.Yes)
    

  def Load(self):
    filename = QtGui.QFileDialog.getOpenFileName(None, 'Open file','')
    if(QtCore.QFileInfo(filename).isFile()):
      FileLoader.LoadObject(str(filename))

  def OpenDocs(self):
    QtGui.QDesktopServices.openUrl(QtCore.QUrl("http://www.openstructure.org/docs/"))
    
  def About(self):
    _InitSplash()
  
  def ShowGLWin(self):
    gosty=gui.GostyApp.Instance()
    gl_win=gosty.GetGLWin()
    if gl_win and gl_win.qobject.isHidden():
      gl_win.Show()
  
  def ResetView(self):
    msg_box = QtGui.QMessageBox()
    msg_box.setWindowTitle("Reset the Panels and Widget");
    msg_box.setIcon(QtGui.QMessageBox.Question)
    msg_box.setText("Do you really want to reset the Panels and Widgets?");
    msg_box.setStandardButtons(QtGui.QMessageBox.Yes | QtGui.QMessageBox.Cancel);
    msg_box.setDefaultButton(QtGui.QMessageBox.Cancel);
    ret = msg_box.exec_();
    if(ret == QtGui.QMessageBox.Yes):
      settings = QtCore.QSettings()
      settings.setValue("restore_settings",QtCore.QVariant(False))
      info_box = QtGui.QMessageBox()
      info_box.setStandardButtons(QtGui.QMessageBox.Ok)
      info_box.setIcon(QtGui.QMessageBox.Information)
      info_box.setWindowTitle("Restart OpenStructure")
      info_box.setText("You must restart OpenStructure for the changes to take effect!");
      info_box.exec_();
            
def _InitMenuBar(app):
  #InitMenuBar(app.perspective.menubar)
  pass
## \example menubar_example.py
#
# Shows how to use PyQt to add a menu from within Python and interact 
# with the currently selected objects in the scene menu.
