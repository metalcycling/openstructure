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
class SlineWidget(RenderModeWidget):
  def __init__(self, parent=None):
    RenderModeWidget.__init__(self, parent)
    
    #Title
    self.text_ = "Fast Spline"
    
    #Set Render Mode
    self.mode_ = gfx.RenderMode.SLINE
    
    #Defaults
    min_detail = 1
    max_detail = 20
    
    min_line_width = 0.01
    max_line_width = 20
        
    #Create Ui elements
    self.detail_spinbox_ = QtWidgets.QSpinBox()
    self.detail_spinbox_.setRange(min_detail, max_detail)
    self.detail_spinbox_.setSingleStep(1)
    
    self.aa_rendering_cb_ = QtWidgets.QCheckBox()
    
    self.radius_spinbox_ = QtWidgets.QDoubleSpinBox()
    self.radius_spinbox_.setRange(min_line_width, max_line_width)
    self.radius_spinbox_.setDecimals(2)
    self.radius_spinbox_.setSingleStep(0.1)
    
    sline_label = QtWidgets.QLabel("Spline Settings")
    font = sline_label.font()
    font.setBold(True)
    
    detail_label = QtWidgets.QLabel("Spline Detail")
    aa_label = QtWidgets.QLabel("AA-Lines")
    radius_label = QtWidgets.QLabel("Line Width")
    grid = QtWidgets.QGridLayout()
    grid.addWidget(sline_label, 0, 0, 1, 1)
    grid.addWidget(detail_label, 1, 0, 1, 3)
    grid.addWidget(self.detail_spinbox_, 1, 2, 1, 1)
    grid.addWidget(aa_label, 2, 0, 1, 3)
    grid.addWidget(self.aa_rendering_cb_, 2, 2, 1, 1)
    grid.addWidget(radius_label, 3, 0, 1, 3)
    grid.addWidget(self.radius_spinbox_, 3, 2, 1, 1)
    grid.setRowStretch(4,1)
    self.setLayout(grid)

    self.detail_spinbox_.valueChanged.connect(self.UpdateDetail)
    self.aa_rendering_cb_.stateChanged.connect(self.UpdateAA)
    self.radius_spinbox_.valueChanged.connect(self.UpdateLineWidth)

    self.setMinimumSize(250,120)
  
  def UpdateDetail(self, value):
    self.GetOptions().SetSplineDetail(value)
    self.ApplyOptions()
    
  def UpdateAA(self, value):
    self.GetOptions().SetAALines(value)
    self.ApplyOptions()
    
  def UpdateLineWidth(self, value):
    self.GetOptions().SetLineWidth(value)
    self.ApplyOptions()
    
  def UpdateGui(self,options):
    self.detail_spinbox_.setValue(options.GetSplineDetail())
    self.aa_rendering_cb_.setChecked(options.GetAALines())    
    self.radius_spinbox_.setValue(options.GetLineWidth())

  def GetText(self):
    return self.text_

  def GetRenderMode(self):
    return self.mode_
