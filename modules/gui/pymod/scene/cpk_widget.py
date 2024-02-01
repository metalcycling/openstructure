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

#CPK Render Options
class CPKWidget(RenderModeWidget):
  def __init__(self, parent=None):
    RenderModeWidget.__init__(self, parent)
    
    #Title
    self.text_ = "Spheres"
    
    #Defaults
    min_sphere_detail = 1
    max_sphere_detail = 20
    
    #Set Render Mode
    self.mode_ = gfx.RenderMode.CPK
        
    #Create Ui elements
    
    self.sphere_spinbox_ = QtWidgets.QSpinBox()
    self.sphere_spinbox_.setRange(min_sphere_detail, max_sphere_detail)
    cpk_mode_label = QtWidgets.QLabel("Rendering Mode")
    self.cpk_mode_ = QtWidgets.QComboBox()
    self.cpk_mode_.addItem("Triangles")
    self.cpk_mode_.addItem("3D Sprites")
    cpk_label = QtWidgets.QLabel(self.text_)
    font = cpk_label.font()
    font.setBold(True)
    
    sphere_label = QtWidgets.QLabel("Sphere Detail")
    grid = QtWidgets.QGridLayout()
    grid.addWidget(cpk_label,0,0,1,1)
    grid.addWidget(sphere_label, 1, 0, 1, 3)
    grid.addWidget(self.sphere_spinbox_, 1, 2, 1, 1)
    grid.setRowStretch(2,1)
    self.setLayout(grid)
  
    self.sphere_spinbox_.valueChanged.connect(self.UpdateSphereDetail)
    self.cpk_mode_.currentIndexChanged.connect(self.UpdateSphereMode)
    
    self.setMinimumSize(250,60)
    
  def UpdateSphereDetail(self, value):
    self.GetOptions().SetSphereDetail(value)
    self.ApplyOptions()

  def UpdateSphereMode(self, value):
    self.GetOptions().SetSphereMode(value)
    self.ApplyOptions()
    
  def UpdateGui(self,options):
    self.sphere_spinbox_.setValue(options.GetSphereDetail())
    self.cpk_mode_.setCurrentIndex(options.GetSphereMode())

  def GetText(self):
    return self.text_

  def GetRenderMode(self):
    return self.mode_
