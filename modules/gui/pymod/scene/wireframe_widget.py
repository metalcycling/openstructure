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

from ost import gui
from ost import gfx
from PyQt5 import QtCore, QtWidgets
try: 
  from ost import img
  _img_present=True
except ImportError:
  _img_present=False
  pass
from .scene_selection_helper import SelHelper

#Wireframe Options
class WireframeWidget(QtWidgets.QWidget):
  def __init__(self, parent=None):
    QtWidgets.QWidget.__init__(self, parent)
    
    #Title
    self.text_ = "Wireframe"
        
    #Defaults
    min_line_width = 0.01
    max_line_width = 20
    
    #Set Render Mode
    self.mode_ = gfx.RenderMode.SIMPLE
        
    #Create Ui elements
    self.aa_rendering_cb_ = QtWidgets.QCheckBox()
  
    self.radius_spinbox_ = QtWidgets.QDoubleSpinBox()
    self.radius_spinbox_.setRange(min_line_width, max_line_width)
    self.radius_spinbox_.setDecimals(2)
    self.radius_spinbox_.setSingleStep(0.1)

    simple_label = QtWidgets.QLabel("Wireframe Settings")
    font = simple_label.font()
    font.setBold(True)
    
    radius_label = QtWidgets.QLabel("Line Width")
    aa_label = QtWidgets.QLabel("AA-Lines")

    grid = QtWidgets.QGridLayout()
    grid.addWidget(simple_label,0,0,1,3)
    grid.addWidget(aa_label, 1, 0, 1, 3)
    grid.addWidget(self.aa_rendering_cb_, 1, 2, 1, 1)
    grid.addWidget(radius_label,2,0,1,3)
    grid.addWidget(self.radius_spinbox_,2,2,1,1)
    grid.setRowStretch(5,1)
    self.setLayout(grid)
    
    self.radius_spinbox_.valueChanged.connect(self.UpdateLineWidth)
    self.aa_rendering_cb_.stateChanged.connect(self.UpdateAA)  
    
    self.setMinimumSize(250,100)
    
  def UpdateAA(self, value):
    scene_selection = gui.SceneSelection.Instance()
    for i in range(0,scene_selection.GetActiveNodeCount()):
      node = scene_selection.GetActiveNode(i)
      node.SetAALines(value)
          
  def UpdateLineWidth(self, value):
    scene_selection = gui.SceneSelection.Instance()
    for i in range(0,scene_selection.GetActiveNodeCount()):
      node = scene_selection.GetActiveNode(i)
      node.SetLineWidth(value)

  def UpdateGui(self):
    scene_selection = gui.SceneSelection.Instance()
    node = scene_selection.GetActiveNode(0)
    self.radius_spinbox_.setValue(node.GetLineWidth())

  def Update(self):
    self.setEnabled(True)
    self.UpdateGui()
    if SelHelper().CheckNotFlags(SelHelper.HAS_IMG | SelHelper.IS_ONE_TYPE):
      self.setEnabled(False)          

  def GetText(self):
    return self.text_

  def GetRenderMode(self):
    return self.mode_
