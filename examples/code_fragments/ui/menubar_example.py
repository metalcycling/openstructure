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
from ost import gui
from ost.gui import FileLoader

from PyQt5 import QtCore, QtGui, QtWidgets

class InitMenuBar(QtCore.QObject):
  def __init__(self, menu_bar=None):
    QtCore.QObject.__init__(self, menu_bar)
    self.scene_selection_ = gui.SceneSelection.Instance()
    
    test_action = QtWidgets.QAction('Test Menu Point', self)
    test_action.setStatusTip('Print Hello World')
    test_action.setShortcut('Ctrl+T')
    test_action.triggered.connect(self.TestMethod)

    test = menu_bar.addMenu('&Test')
    test.addAction(test_action)

  def TestMethod(self):
    reply = QtWidgets.QMessageBox()
    
    node_count = self.scene_selection_.GetActiveNodeCount()
    if(node_count > 0):
      string = "";
      for i in range(node_count):
        entity = self.scene_selection_.GetActiveNode(i)
        string=string + entity.GetName()
      reply.setText("Oh, there are selected entities: %s" % string)
    else:
      reply.setText("This is a test!")
    reply.addButton(QtWidgets.QMessageBox.Yes)
    reply.exec_()
    
menu_bar=gui.GostyApp.Instance().perspective.GetMenuBar()
InitMenuBar(menu_bar)
