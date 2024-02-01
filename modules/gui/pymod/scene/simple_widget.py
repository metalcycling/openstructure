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
from .render_mode_widget import RenderModeWidget

#Simple Render Options
class SimpleWidget(RenderModeWidget):
  def __init__(self, parent=None):
    RenderModeWidget.__init__(self, parent)
    
    #Title
    self.text_ = "Fast Bonds"
        
    #Defaults
    min_line_width = 0.01
    max_line_width = 20
    min_bo_dist = 0.01
    max_bo_dist = 0.20
    
    #Set Render Mode
    self.mode_ = gfx.RenderMode.SIMPLE
        
    #Create Ui elements
    self.aa_rendering_cb_ = QtWidgets.QCheckBox()
  
    self.radius_spinbox_ = QtWidgets.QDoubleSpinBox()
    self.radius_spinbox_.setRange(min_line_width, max_line_width)
    self.radius_spinbox_.setDecimals(2)
    self.radius_spinbox_.setSingleStep(0.1)

    self.bo_rendering_cb_ = QtWidgets.QCheckBox()
    
    self.bo_distance_spinbox_ = QtWidgets.QDoubleSpinBox()
    self.bo_distance_spinbox_.setRange(min_bo_dist, max_bo_dist)
    self.bo_distance_spinbox_.setDecimals(2)
    self.bo_distance_spinbox_.setSingleStep(0.01)

    simple_label = QtWidgets.QLabel("Fast Bonds Simple Settings")
    font = simple_label.font()
    font.setBold(True)
    
    radius_label = QtWidgets.QLabel("Line Width")
    aa_label = QtWidgets.QLabel("AA-Lines")
    bo_label = QtWidgets.QLabel("Show Bond Order")
    bo_distance = QtWidgets.QLabel("Bond Order Distance")

    grid = QtWidgets.QGridLayout()
    grid.addWidget(simple_label,0,0,1,3)
    grid.addWidget(aa_label, 1, 0, 1, 3)
    grid.addWidget(self.aa_rendering_cb_, 1, 2, 1, 1)
    grid.addWidget(radius_label,2,0,1,3)
    grid.addWidget(self.radius_spinbox_,2,2,1,1)
    grid.addWidget(bo_label, 3, 0, 1, 3)
    grid.addWidget(self.bo_rendering_cb_, 3,2,1,1)
    grid.addWidget(bo_distance, 4, 0, 1, 3)
    grid.addWidget(self.bo_distance_spinbox_, 4,2,1,1)
    grid.setRowStretch(5,1)
    self.setLayout(grid)
    
    self.radius_spinbox_.valueChanged.connect(self.UpdateLineWidth)
    self.aa_rendering_cb_.stateChanged.connect(self.UpdateAA)
    self.bo_distance_spinbox_.valueChanged.connect(self.UpdateBODistance)
    self.bo_rendering_cb_.stateChanged.connect(self.UpdateBO)    
    self.setMinimumSize(250,140)
    
    self.setMinimumSize(250,140)
    
  def UpdateAA(self, value):
    self.GetOptions().SetAALines(value)
    self.ApplyOptions()
    
  def UpdateLineWidth(self, value):
    self.GetOptions().SetLineWidth(value)
    self.ApplyOptions()
    
  def UpdateBO(self, value):
    self.GetOptions().SetBondOrderFlag(value)
    self.ApplyOptions()
    
  def UpdateBODistance(self, value):
    self.GetOptions().SetBondOrderDistance(value)
    self.ApplyOptions()

  def UpdateGui(self,options):
    self.aa_rendering_cb_.setChecked(options.GetAALines())    
    self.radius_spinbox_.setValue(options.GetLineWidth())
    self.bo_rendering_cb_.setChecked(options.GetBondOrderFlag())    
    self.bo_distance_spinbox_.setValue(options.GetBondOrderDistance())

  def GetText(self):
    return self.text_

  def GetRenderMode(self):
    return self.mode_
