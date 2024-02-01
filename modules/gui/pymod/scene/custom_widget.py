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

#Custom Render Options
class CustomWidget(RenderModeWidget):
  def __init__(self, parent=None):
    RenderModeWidget.__init__(self, parent)
    
    #Title
    self.text_ = "Ball & Stick"
    
    #Defaults
    min_sphere_detail = 1
    max_sphere_detail = 12
    
    min_arc_detail = 1
    max_arc_detail = 12
    
    min_rad = 0.05
    max_bond_rad = 5.0
    max_sphere_rad = 5.0
    
    
    
    #Set Render Mode
    self.mode_ = gfx.RenderMode.CUSTOM
        
    #Create Ui elements
    
    self.sphere_spinbox_ = QtWidgets.QSpinBox()
    self.sphere_spinbox_.setRange(min_sphere_detail, max_sphere_detail)
    
    self.arc_spinbox_ = QtWidgets.QSpinBox()
    self.arc_spinbox_.setRange(min_arc_detail, max_arc_detail)
  
    #Bond Radius
    radius_bond_label = QtWidgets.QLabel("Bond Radius")
    
    self.radius_bond_spinbox_ = QtWidgets.QDoubleSpinBox()
    self.radius_bond_spinbox_.setRange(min_rad, max_bond_rad)
    self.radius_bond_spinbox_.setDecimals(2)
    self.radius_bond_spinbox_.setSingleStep(0.05)
    
    self.radius_bond_slider_ = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
    self.radius_bond_slider_.setRange(int(min_rad*100.0), int(max_bond_rad*100.0))
    self.radius_bond_slider_.setTickPosition(QtWidgets.QSlider.NoTicks)
    self.radius_bond_slider_.setTickInterval(5)
    
    
    #Sphere Radius
    radius_sphere_label = QtWidgets.QLabel("Sphere Radius")
    
    self.radius_sphere_spinbox_ = QtWidgets.QDoubleSpinBox()
    self.radius_sphere_spinbox_.setRange(min_rad, max_sphere_rad)
    self.radius_sphere_spinbox_.setDecimals(2)
    self.radius_sphere_spinbox_.setSingleStep(0.05)
    
    self.radius_sphere_slider_ = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
    self.radius_sphere_slider_.setRange(int(min_rad*100.0), int(max_sphere_rad*100.0))
    self.radius_sphere_slider_.setTickPosition(QtWidgets.QSlider.NoTicks)
    self.radius_sphere_slider_.setTickInterval(5)
    
    custom_label = QtWidgets.QLabel(self.text_)
    font = custom_label.font()
    font.setBold(True)
    
    sphere_label = QtWidgets.QLabel("Sphere Detail")
    arc_label = QtWidgets.QLabel("Arc Detail")
    grid = QtWidgets.QGridLayout()
    grid.addWidget(custom_label,0,0,1,1)
    grid.addWidget(sphere_label, 1, 0, 1, 3)
    grid.addWidget(self.sphere_spinbox_, 1, 2, 1, 1)
    grid.addWidget(arc_label,2,0,1,3)
    grid.addWidget(self.arc_spinbox_,2,2,1,1)
    
    grid.addWidget(radius_bond_label, 3, 0, 1, 1)
    grid.addWidget(self.radius_bond_slider_, 3, 1, 1, 3)
    grid.addWidget(self.radius_bond_spinbox_, 3, 4, 1, 1)
    
    grid.addWidget(radius_sphere_label, 4, 0, 1, 1)
    grid.addWidget(self.radius_sphere_slider_, 4, 1, 1, 3)
    grid.addWidget(self.radius_sphere_spinbox_, 4, 4, 1, 1)
    
    grid.setRowStretch(3,1)
    self.setLayout(grid)
  
    self.sphere_spinbox_.valueChanged.connect(self.UpdateSphereDetail)
    self.arc_spinbox_.valueChanged.connect(self.UpdateArcDetail)
    
    self.radius_bond_spinbox_.valueChanged.connect(self.UpdateBondRadius)
    self.radius_bond_slider_.valueChanged.connect(self.UpdateSliderBondRadius)
    
    self.radius_sphere_spinbox_.valueChanged.connect(self.UpdateSphereRadius)
    self.radius_sphere_slider_.valueChanged.connect(self.UpdateSliderSphereRadius)
    
    self.setMinimumSize(250,150)
    
  def UpdateSphereDetail(self, value):
    self.GetOptions().SetSphereDetail(value)
    self.ApplyOptions()

  def UpdateArcDetail(self, value):
    self.GetOptions().SetArcDetail(value)
    self.ApplyOptions()

  def UpdateBondRadius(self, value):
    self.GetOptions().SetBondRad(value)
    if(self.GetOptions().GetSphereRad()<self.GetOptions().GetBondRad()):
      self.GetOptions().SetSphereRad(value)
    self.ApplyOptions()
    
  def UpdateSliderBondRadius(self, value):
    self.GetOptions().SetBondRad(value/100.0)
    self.ApplyOptions()
    
  def UpdateSphereRadius(self, value):
    self.GetOptions().SetSphereRad(value)
    if(self.GetOptions().GetSphereRad()<self.GetOptions().GetBondRad()):
      self.GetOptions().SetBondRad(value)
    self.ApplyOptions()
    
  def UpdateSliderSphereRadius(self, value):
    self.GetOptions().SetSphereRad(value/100.0)
    self.ApplyOptions()

  def UpdateSphereRadiusGui(self,value):
    value = round(value, 2)
    if(abs(value*100.0 - self.radius_sphere_slider_.value())>=self.radius_sphere_slider_.singleStep()):
      self.radius_sphere_slider_.setValue(int(value*100.0))
    if (abs(value - self.radius_sphere_spinbox_.value())>=self.radius_sphere_spinbox_.singleStep()):
      self.radius_sphere_spinbox_.setValue(value)

  def UpdateBondRadiusGui(self,value):
    value = round(value, 2)
    if(abs(value*100.0 - self.radius_bond_slider_.value())>=self.radius_bond_slider_.singleStep()):
      self.radius_bond_slider_.setValue(int(value*100.0))
    if (abs(value - self.radius_bond_spinbox_.value())>=self.radius_bond_spinbox_.singleStep()):
      self.radius_bond_spinbox_.setValue(value)

  def UpdateGui(self,options):
    self.sphere_spinbox_.setValue(options.GetSphereDetail())    
    self.arc_spinbox_.setValue(options.GetArcDetail())
    
    self.UpdateBondRadiusGui(options.GetBondRad())
    self.UpdateSphereRadiusGui(options.GetSphereRad())

  def GetText(self):
    return self.text_

  def GetRenderMode(self):
    return self.mode_
