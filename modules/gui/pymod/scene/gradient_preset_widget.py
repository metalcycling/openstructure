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
# -*- coding: utf-8 -*-

import os
import ost
from ost import gui
from ost import gfx
from ost import info
from datetime import datetime

from datetime import datetime
from PyQt5 import QtCore, QtWidgets, QtGui
from .gradient_info_handler import GradientInfoHandler
from .gradient_list_model import GradientListModel


#Gradient Preset Widget
class GradientPresetWidget(QtWidgets.QWidget):
  ICONS_DIR = os.path.join(ost.GetSharedDataPath(), "gui", "icons/")

  gradientSelected = QtCore.pyqtSignal(object, name="gradientSelected")

  def __init__(self, gradient_edit, parent=None):
    QtWidgets.QWidget.__init__(self, parent)
    
    #Title
    self.text_ = "Gradient Presets"
    
    #Refrences
    self.gradient_edit_ = gradient_edit   

    #Create Ui elements
    self.list_view_ = QtWidgets.QListView()
    
    #Create Model
    self.list_model_ = GradientListModel(self)
    self.list_view_.setModel(self.list_model_)
    self.list_view_.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        
    preset_label = QtWidgets.QLabel(self.text_)
    font = preset_label.font()
    font.setBold(True)
    
    self.add_action = QtWidgets.QAction("+",self)
    self.add_action.setIcon(QtGui.QIcon(GradientPresetWidget.ICONS_DIR+"add_icon.png"))
    
    self.add_action.triggered.connect(self.Add)
    
    self.add_button_ = QtWidgets.QToolButton(self)
    self.add_button_.setIconSize(QtCore.QSize(20,20))
    self.add_button_.setDefaultAction(self.add_action)
        
    grid = QtWidgets.QGridLayout()
    grid.setContentsMargins(0,5,0,0)
    grid.addWidget(preset_label, 0, 0, 1, 1)
    qhbox = QtWidgets.QHBoxLayout()
    grid.addWidget(self.list_view_,1,0,3,3)
    grid.addWidget(self.add_button_,4,0,1,1)
    self.setLayout(grid)
    
    self.list_view_.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
    self.list_view_.customContextMenuRequested.connect(self.contextMenuEvent)
    self.CreateImmutableContextMenu()
    self.CreateContextMenu()

    self.list_view_.doubleClicked.connect(self.Load)

  def CreateImmutableContextMenu(self):  
    self.immucontextMenu_ = QtWidgets.QMenu("Context menu", self)
    self.load_ = QtWidgets.QAction("Load", self.list_view_)  
    self.immucontextMenu_.addAction(self.load_)
    #Connect Signal with Slot  
    self.load_.triggered.connect(self.LoadCurrentIndex)

  def CreateContextMenu(self):
    self.contextMenu_ = QtWidgets.QMenu("Context menu", self)
    self.remove_ = QtWidgets.QAction("Remove", self.list_view_)
    self.rename_ = QtWidgets.QAction("Rename", self.list_view_)
    self.contextMenu_.addAction(self.load_)
    self.contextMenu_.addAction(self.remove_)
    self.contextMenu_.addAction(self.rename_)
    #Connect Signals with Slots  
    self.remove_.triggered.connect(self.Remove)
    self.rename_.triggered.connect(self.Rename)
  
  def contextMenuEvent(self, pos):
    #ContextMenu
    index = self.list_view_.indexAt(pos)
    if index.isValid(): 
      if self.list_model_.IsEditable(index.row()):
        self.contextMenu_.popup(QtWidgets.QCursor.pos())
      else:
        self.immucontextMenu_.popup(QtWidgets.QCursor.pos())
  
  def Add(self):
    if(self.list_view_.currentIndex().isValid()):
      self.list_view_.closePersistentEditor(self.list_view_.currentIndex())
    row = self.list_model_.GetLastRow()
    if self.list_model_.AddItem(datetime.now().isoformat(' '), self.gradient_edit_.GetGradient(), row, True, True):
      index = self.list_model_.index(row)
      self.list_view_.setCurrentIndex(index)
      self.Rename()
    else:
      QtWidgets.QMessageBox.information(self, "Gradient not added", "The gradient could not be added!")
      
  def Remove(self):
    if(self.list_view_.currentIndex().isValid()):
      ret = QtWidgets.QMessageBox.warning(self, "Delete Gradient",
                   "Delete Gradient?",
                   QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
      if ret == QtWidgets.QMessageBox.Yes:
        self.list_model_.RemoveItem(self.list_view_.currentIndex().row())

  def LoadCurrentIndex(self):
    if(self.list_view_.currentIndex().isValid()):
      self.Load(self.list_view_.currentIndex())

  def Load(self, index):
    if(index.isValid()):
      self.gradientSelected.emit(self.list_model_.GetGradient(index))

  def Rename(self):
    if(self.list_view_.currentIndex().isValid()):
      self.list_view_.edit(self.list_view_.currentIndex())
      
  def GetText(self):
    return self.text_
      
