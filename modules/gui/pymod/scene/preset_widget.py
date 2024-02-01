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
from ost import gfx
import ost
import os
from datetime import datetime
from PyQt5 import QtCore, QtWidgets, QtGui
from .scene_selection_helper import SelHelper
from .preset_list_model import PresetListModel
from .preset_editor_widget import PresetEditor
from .preset import Preset

class PresetWidget(QtWidgets.QWidget):
  PRESET_XML_FILE = os.path.join(ost.GetSharedDataPath(), "scene", "presets.xml")
  ICONS_DIR = os.path.join(ost.GetSharedDataPath(), "gui", "icons/")
  def __init__(self, parent=None):   
    QtWidgets.QWidget.__init__(self, parent)
    
    self.text_ = "Presets"
        
    #Create Ui elements
    self.list_view_ = QtWidgets.QListView()
        
    #Create Model
    self.list_model_ = PresetListModel(self)
    self.list_view_.setModel(self.list_model_)
    self.list_view_.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
    
    self.add_action = QtWidgets.QAction("+",self)
    self.add_action.setIcon(QtGui.QIcon(PresetWidget.ICONS_DIR+"add_icon.png"))
    
    self.add_action.triggered.connect(self.Add)
    
    self.add_button_ = QtWidgets.QToolButton(self)
    self.add_button_.setIconSize(QtCore.QSize(20,20))
    self.add_button_.setDefaultAction(self.add_action)
    
    self.preset_editor_ = PresetEditor(self)
    
    grid = QtWidgets.QGridLayout()
    grid.setContentsMargins(0,5,0,0)
    grid.addWidget(self.list_view_,0,0,3,3)
    grid.addWidget(self.add_button_,3,0,1,1)
    self.setLayout(grid)
    
    self.list_view_.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
    self.list_view_.customContextMenuRequested.connect(self.contextMenuEvent)
    self.CreateImmutableContextMenu()
    self.CreateContextMenu()
  
    self.list_view_.doubleClicked.connect(self.Load)
    
    self.setMinimumSize(250,200)
    
  def CreateImmutableContextMenu(self):  
    self.immucontext_menu_ = QtWidgets.QMenu("Context menu", self)
    self.load_ = QtWidgets.QAction("Load", self.list_view_)  
    self.immucontext_menu_.addAction(self.load_)
    #Connect Signal with Slot  
    self.load_.triggered.connect(self.LoadCurrentIndex)

  def CreateContextMenu(self):
    self.context_menu_ = QtWidgets.QMenu("Context menu", self)
    self.remove_ = QtWidgets.QAction("Remove", self.list_view_)
    self.rename_ = QtWidgets.QAction("Rename", self.list_view_)
    self.edit_ = QtWidgets.QAction("Edit", self.list_view_)
    self.context_menu_.addAction(self.load_)
    self.context_menu_.addAction(self.remove_)
    self.context_menu_.addAction(self.rename_)
    self.context_menu_.addAction(self.edit_)
    #Connect Signals with Slots  
    self.remove_.triggered.connect(self.Remove)
    self.rename_.triggered.connect(self.Rename)
    self.edit_.triggered.connect(self.Edit)
  
  def contextMenuEvent(self, pos):
    #ContextMenu
    index = self.list_view_.indexAt(pos)
    if index.isValid(): 
      if self.list_model_.IsEditable(index.row()):
        self.context_menu_.popup(QtWidgets.QCursor.pos())
      else:
        self.immucontext_menu_.popup(QtWidgets.QCursor.pos())
  
  def Add(self):
    row = self.list_model_.GetLastRow()
    preset = Preset(datetime.now().isoformat(' '))
    self.preset_editor_.SetPreset(preset)
    if(self.preset_editor_.exec_()):
      if self.list_model_.AddItem(preset, row, True, True):
        index = self.list_model_.index(row)
        self.list_view_.setCurrentIndex(index)
        self.Rename()
      
  def Remove(self):
    if(self.list_view_.currentIndex().isValid()):
      ret = QtWidgets.QMessageBox.warning(self, "Delete Preset",
                   "Delete Preset?",
                   QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
      if ret == QtWidgets.QMessageBox.Yes:
        self.list_model_.RemoveItem(self.list_view_.currentIndex().row())

  def Edit(self):
    if(self.list_view_.currentIndex().isValid()):
        preset = self.list_model_.GetPreset(self.list_view_.currentIndex())
        self.preset_editor_.SetPreset(preset)
        if(self.preset_editor_.exec_()):
            row = self.list_view_.currentIndex().row()
            self.list_model_.RemoveItem(row)
            self.list_model_.AddItem(preset, row, True, True)
        
  def LoadCurrentIndex(self):
    if(self.list_view_.currentIndex().isValid()):
      self.Load(self.list_view_.currentIndex())

  def Load(self, index):
    if(index.isValid()):
      scene_selection = gui.SceneSelection.Instance()
      preset=self.list_model_.GetPreset(index)
      for i in range(0,scene_selection.GetActiveNodeCount()):
        node = scene_selection.GetActiveNode(i)
        if isinstance(node, gfx.Entity):
          node.CleanColorOps()
          preset.ApplyOn(node)
          
  def ChangeColor(self,node):
    self.LoadCurrentIndex()
  
  def Update(self):
    self.setEnabled(True)
    if SelHelper().CheckAllFlags(SelHelper.NO_SELECTION):
      self.setEnabled(False)
      return
    
    if SelHelper().CheckNotFlags(SelHelper.HAS_ENTITY | SelHelper.IS_ONE_TYPE):
      self.setEnabled(False)
      return

  def Rename(self):
    if(self.list_view_.currentIndex().isValid()):
      self.list_view_.edit(self.list_view_.currentIndex())
    
  def GetText(self):
    return self.text_
