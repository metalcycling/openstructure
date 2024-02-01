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

import sys
from ost import mol
from ost import gui
from ost import gfx
try: 
  from ost import img
  _img_present=True
except ImportError:
  _img_present=False
  pass

from PyQt5 import QtCore, QtWidgets, QtWidgets
from .scene_selection_helper import SelHelper
from .gradient_editor_widget import GradientEditor
from .uniform_color_widget import UniformColorWidget
from .combo_options_widget import ComboOptionsWidget

class ColorOptionsWidget(ComboOptionsWidget):
  def __init__(self, parent=None):
    ComboOptionsWidget.__init__(self, parent)
    
    #Title
    self.text_ = "Color Options"
    conly_label_ = QtWidgets.QLabel('carbons only')
    self.conly_box_ = QtWidgets.QCheckBox()
    
    #Add options to menu
    self.entity_widgets_ = list()
    self.entity_widgets_.append(["Color by Element", ByElementWidget("Color by Element", self)])
    self.entity_widgets_.append(["Color by Chain", ByChainWidget("Color by Chain", self)])
    self.entity_widgets_.append(["Color by Entity", ByEntityWidget("Color by Entity", self)])
    self.entity_widgets_.append(["Color by Property", GradientEditor(self)])
    self.entity_widgets_.append(["Uniform",UniformColorWidget(self)])
  
    self.img_widgets_ = list()
    self.img_widgets_.append(["Uniform",UniformColorWidget()])

    qw = QtWidgets.QWidget(self)
    gl = QtWidgets.QGridLayout(qw)
    gl.addWidget(self.conly_box_, 0, 0, 1, 1)
    gl.addWidget(conly_label_, 0, 1, 1, 4)
    self.grid_layout_.addWidget(qw, 2, 0, 1, 1)

    self.setMinimumSize(250,200)
    
  def OnComboChange(self, item):
    scene_selection = gui.SceneSelection.Instance()
    if hasattr(item,"PrepareColoring"):
      item.PrepareColoring()
    for i in range(0,scene_selection.GetActiveNodeCount()):
      node = scene_selection.GetActiveNode(i)
      item.ChangeColor(node)
    
    if(scene_selection.GetActiveViewCount() > 0):
      entity = scene_selection.GetViewEntity()
      view = scene_selection.GetViewUnion()
      item.ChangeViewColor(entity,view)

    self.DoResize()
    
    
  def Update(self):
    
    ComboOptionsWidget.setEnabled(self,True)
    
    if SelHelper().CheckAllFlags(SelHelper.NO_SELECTION):
      ComboOptionsWidget.setEnabled(self,False)
      return
    
    for w in self.entity_widgets_:
      self.RemoveWidget(w[0])
    for w in self.img_widgets_:
      self.RemoveWidget(w[0])
    
    
    if SelHelper().CheckFlags(SelHelper.HAS_IMG | SelHelper.IS_ONE_TYPE):
      for w in self.img_widgets_:
        self.AddWidget(w[0], w[1])
    elif SelHelper().CheckMinOneFlag(SelHelper.HAS_ENTITY| SelHelper.HAS_VIEW| SelHelper.HAS_SURFACE) and SelHelper().CheckNotFlags(SelHelper.HAS_IMG):
      for w in self.entity_widgets_:
        self.AddWidget(w[0], w[1])
    else:
      ComboOptionsWidget.setEnabled(self,False)
      return
    
    self.GetCurrentWidget().Update()
    
  def GetText(self):
    return self.text_

  def GetCarbonsOnly(self):
    return self.conly_box_.isChecked()
  

class ByElementWidget(QtWidgets.QWidget):
  def __init__(self, text, parent=None):
    QtWidgets.QWidget.__init__(self, parent)
    self.parent_ = parent

    
    #Title
    self.text_ = text
    
    #UI
    text_label = QtWidgets.QLabel(text)
    font = text_label.font()
    font.setBold(True)
    
    grid = QtWidgets.QGridLayout()
    grid.addWidget(text_label,0,0,1,1)
    grid.addWidget(QtWidgets.QLabel("No Settings available"), 1, 0, 1, 3)
    grid.setRowStretch(2,1)
    self.setLayout(grid)
    self.setMinimumSize(250,60)
    
  def Update(self):
   pass #Do Nothing
    
  def ChangeColor(self, node):
    if isinstance(node, gfx.Entity):
      node.CleanColorOps()
      if self.parent_.GetCarbonsOnly():
        node.ColorByElement("ele=C")
      else:
        node.ColorByElement()
      
  def ChangeViewColor(self, entity, view):
    if isinstance(entity, gfx.Entity) and isinstance(view, mol.EntityView):
      if self.parent_.GetCarbonsOnly():
        beco=gfx.ByElementColorOp(mol.QueryViewWrapper(mol.Query("ele=C"), view))
      else:
        beco=gfx.ByElementColorOp(mol.QueryViewWrapper(view))
      entity.Apply(beco)
      
  def GetText(self):
    return self.text_
  
  
class ByChainWidget(QtWidgets.QWidget):
  def __init__(self, text, parent=None):
    QtWidgets.QWidget.__init__(self, parent)
    self.parent_ = parent
    
    #Title
    self.text_ = text
    
    #UI
    text_label = QtWidgets.QLabel(text)
    font = text_label.font()
    font.setBold(True)
    
    grid = QtWidgets.QGridLayout()
    grid.addWidget(text_label,0,0,1,1)
    grid.setRowStretch(2,1)
    self.setLayout(grid)
    self.setMinimumSize(250,60)
    
  def Update(self):
   pass #Do Nothing
    
  def ChangeColor(self, node):
    if isinstance(node, gfx.Entity):
      node.CleanColorOps()
      if self.parent_.GetCarbonsOnly():
        node.ColorByChain('ele=C')
      else:
        node.ColorByChain()
      
  def ChangeViewColor(self, entity, view):
    if isinstance(entity, gfx.Entity) and isinstance(view, mol.EntityView):
      if self.parent_.GetCarbonsOnly():
        bco=gfx.ByChainColorOp(mol.QueryViewWrapper(mol.Query("ele=C"),view))
      else:
        bco=gfx.ByChainColorOp(mol.QueryViewWrapper(view))
      entity.Apply(bco)
      
  def GetText(self):
    return self.text_

class ByEntityWidget(QtWidgets.QWidget):
  def __init__(self, text, parent=None):
    QtWidgets.QWidget.__init__(self, parent)
    self.parent_ = parent
    
    #Title
    self.text_ = text
    
    #UI
    text_label = QtWidgets.QLabel(text)
    font = text_label.font()
    font.setBold(True)
    
    grid = QtWidgets.QGridLayout()
    grid.addWidget(text_label,0,0,1,1)
    grid.setRowStretch(2,1)
    self.setLayout(grid)
    self.setMinimumSize(250,60)
    
    self.gradient_ = gfx.Gradient("RAINBOW")
    
  def Update(self):
   pass #Do Nothing
    
  def PrepareColoring(self):
    scene_selection = gui.SceneSelection.Instance()
    entity_count = scene_selection.GetActiveNodeCount()
    for i in range(0,scene_selection.GetActiveNodeCount()):
      if i<=0:
        color=self.gradient_.GetColorAt(0.0)
      else:
        color=self.gradient_.GetColorAt(float(i) / entity_count)
      node = scene_selection.GetActiveNode(i)
      if self.parent_.GetCarbonsOnly():
        node.SetColor(color, 'ele=C')
      else:
        node.SetColor(color)
  
  def ChangeColor(self, node):
    pass
     
  def ChangeViewColor(self, entity, view):
    if isinstance(entity, gfx.Entity) and isinstance(view, mol.EntityView):
      if self.parent_.GetCarbonsOnly():
        bco=gfx.ByChainColorOp(mol.QueryViewWrapper(mol.Query("ele=C"),view))
      else:
        bco=gfx.ByChainColorOp(mol.QueryViewWrapper(view))
      entity.Apply(bco)
      
  def GetText(self):
    return self.text_
