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

#Tube Render Options
class HSCWidget(RenderModeWidget):
  def __init__(self, parent=None):
    RenderModeWidget.__init__(self, parent)
    
    #Title
    self.text_ = "Helix & Strand Cartoon"
     
    #Set Render Mode
    self.mode_ = gfx.RenderMode.HSC
    
    #Defaults
    min_spline_detail = 1
    max_spline_detail = 20
    
    min_arc_detail = 1
    max_arc_detail = 20
    
    min_width = 0.1
    max_width = 5.0
    max_tube_width= 5.0

    min_tube_ratio = 0.2
    max_tube_ratio = 4
    min_ratio = 0.1
    max_ratio = 2.5
    
    min_ecc = 0.1
    max_ecc = 5

    min_profile=0
    max_profile=4
    
    #########UI##########
    
    #Poly Mode
    poly_mode_label = QtWidgets.QLabel("Poly Mode")
    
    self.poly_mode_cb_ = QtWidgets.QComboBox()
    self.poly_mode_cb_.addItem("Points")
    self.poly_mode_cb_.addItem("Wireframe")
    self.poly_mode_cb_.addItem("Surface")
    
    #Sphere Label
    spline_label = QtWidgets.QLabel("Spline Detail")
    
    self.spline_spinbox_ = QtWidgets.QSpinBox()
    self.spline_spinbox_.setRange(min_spline_detail, max_spline_detail)
    
    #Arc Label
    arc_label = QtWidgets.QLabel("Arc Detail")
    
    self.arc_spinbox_ = QtWidgets.QSpinBox()
    self.arc_spinbox_.setRange(min_arc_detail, max_arc_detail)
    
    #Tube
    tube_label = QtWidgets.QLabel("Tube")
    font = tube_label.font()
    font.setBold(True)
    
    #Tube Radius
    radius_tube_label = QtWidgets.QLabel("Radius")
    
    self.width_tube_spinbox_ = QtWidgets.QDoubleSpinBox()
    self.width_tube_spinbox_.setRange(min_width, max_tube_width)
    self.width_tube_spinbox_.setDecimals(1)
    self.width_tube_spinbox_.setSingleStep(0.1)
    
    self.width_tube_slider_ = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
    self.width_tube_slider_.setRange(int(min_width*10.0), int(max_tube_width*10.0))
    self.width_tube_slider_.setTickPosition(QtWidgets.QSlider.NoTicks)
    self.width_tube_slider_.setTickInterval(1)
    
    #Tube Ratio
    ratio_tube_label = QtWidgets.QLabel("Ratio")
    
    self.thickness_tube_spinbox_ = QtWidgets.QDoubleSpinBox()
    self.thickness_tube_spinbox_.setRange(min_tube_ratio, max_tube_ratio)
    self.thickness_tube_spinbox_.setDecimals(1)
    self.thickness_tube_spinbox_.setSingleStep(0.1)
    
    self.thickness_tube_slider_ = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
    self.thickness_tube_slider_.setRange(int(min_tube_ratio*10.0), int(max_tube_ratio*10))
    self.thickness_tube_slider_.setTickPosition(QtWidgets.QSlider.NoTicks)
    self.thickness_tube_slider_.setTickInterval(1)

    # Tube Profile Type
    tube_profile_label = QtWidgets.QLabel("Tube Profile Type")
    self.tube_profile_spinbox_ = QtWidgets.QSpinBox()
    self.tube_profile_spinbox_.setRange(min_profile, max_profile)
    
    # Helix
    helix_label = QtWidgets.QLabel("Helix")
    font = helix_label.font()
    font.setBold(True)
    
    #Helix Radius
    radius_helix_label = QtWidgets.QLabel("Width")
    
    self.width_helix_spinbox_ = QtWidgets.QDoubleSpinBox()
    self.width_helix_spinbox_.setRange(min_width, max_width)
    self.width_helix_spinbox_.setDecimals(1)
    self.width_helix_spinbox_.setSingleStep(0.1)
    
    self.width_helix_slider_ = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
    self.width_helix_slider_.setRange(int(min_width*10.0), int(max_width*10.0))
    self.width_helix_slider_.setTickPosition(QtWidgets.QSlider.NoTicks)
    self.width_helix_slider_.setTickInterval(1)
    
    #Helix Ratio
    ratio_helix_label = QtWidgets.QLabel("Thickness")
    
    self.thickness_helix_spinbox_ = QtWidgets.QDoubleSpinBox()
    self.thickness_helix_spinbox_.setRange(min_ratio,max_ratio)
    self.thickness_helix_spinbox_.setDecimals(1)
    self.thickness_helix_spinbox_.setSingleStep(0.1)
    
    self.thickness_helix_slider_ = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
    self.thickness_helix_slider_.setRange(int(min_ratio*10.0), int(max_ratio*10.0))
    self.thickness_helix_slider_.setTickPosition(QtWidgets.QSlider.NoTicks)
    self.thickness_helix_slider_.setTickInterval(1)
    
    #Helix ECC
    ecc_helix_label = QtWidgets.QLabel("ECC")
    
    self.ecc_helix_spinbox_ = QtWidgets.QDoubleSpinBox()
    self.ecc_helix_spinbox_.setRange(min_ecc,max_ecc)
    self.ecc_helix_spinbox_.setDecimals(1)
    self.ecc_helix_spinbox_.setSingleStep(0.1)
    
    self.ecc_helix_slider_ = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
    self.ecc_helix_slider_.setRange(int(min_ecc*10), int(max_ecc*10))
    self.ecc_helix_slider_.setTickPosition(QtWidgets.QSlider.NoTicks)
    self.ecc_helix_slider_.setTickInterval(1)
    
    # Helix Profile Type
    helix_profile_label = QtWidgets.QLabel("Helix Profile Type")
    self.helix_profile_spinbox_ = QtWidgets.QSpinBox()
    self.helix_profile_spinbox_.setRange(min_profile, max_profile)
    
    # Helix Rendering Mode
    helix_mode_label = QtWidgets.QLabel("Helix Render Mode")
    self.helix_mode_spinbox_ = QtWidgets.QSpinBox()
    self.helix_mode_spinbox_.setRange(0, 1)
    
    #Strand
    strand_label = QtWidgets.QLabel("Strand")
    font = strand_label.font()
    font.setBold(1)
    
    #Strand Radius
    radius_strand_label = QtWidgets.QLabel("Width")
    
    self.width_strand_spinbox_ = QtWidgets.QDoubleSpinBox()
    self.width_strand_spinbox_.setRange(min_width, max_width)
    self.width_strand_spinbox_.setDecimals(1)
    self.width_strand_spinbox_.setSingleStep(0.1)
    
    self.width_strand_slider_ = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
    self.width_strand_slider_.setRange(int(min_width*10), int(max_width*10))
    self.width_strand_slider_.setTickPosition(QtWidgets.QSlider.NoTicks)
    self.width_strand_slider_.setTickInterval(1)
    
    #Strand Ratio
    ratio_strand_label = QtWidgets.QLabel("Thickness")
    
    self.thickness_strand_spinbox_ = QtWidgets.QDoubleSpinBox()
    self.thickness_strand_spinbox_.setRange(min_ratio,max_ratio)
    self.thickness_strand_spinbox_.setDecimals(1)
    self.thickness_strand_spinbox_.setSingleStep(0.1)
    
    self.thickness_strand_slider_ = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
    self.thickness_strand_slider_.setRange(int(min_ratio*10),int(max_ratio*10))
    self.thickness_strand_slider_.setTickPosition(QtWidgets.QSlider.NoTicks)
    self.thickness_strand_slider_.setTickInterval(1)
    
    #Strand ECC
    ecc_strand_label = QtWidgets.QLabel("ECC")
    
    self.ecc_strand_spinbox_ = QtWidgets.QDoubleSpinBox()
    self.ecc_strand_spinbox_.setRange(min_ecc,max_ecc)
    self.ecc_strand_spinbox_.setDecimals(1)
    self.ecc_strand_spinbox_.setSingleStep(0.1)
    
    self.ecc_strand_slider_ = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
    self.ecc_strand_slider_.setRange(int(min_ecc*10), int(max_ecc*10))
    self.ecc_strand_slider_.setTickPosition(QtWidgets.QSlider.NoTicks)
    self.ecc_strand_slider_.setTickInterval(1)  
    
    # Strand Profile Type
    strand_profile_label = QtWidgets.QLabel("Strand Profile Type")
    self.strand_profile_spinbox_ = QtWidgets.QSpinBox()
    self.strand_profile_spinbox_.setRange(min_profile, max_profile)

    # Strand Rendering Mode
    strand_mode_label = QtWidgets.QLabel("Strand Render Mode")
    self.strand_mode_spinbox_ = QtWidgets.QSpinBox()
    self.strand_mode_spinbox_.setRange(0, 1)

    #
    row=1
    grid = QtWidgets.QGridLayout()
    grid.addWidget(poly_mode_label,row,0,1,1)
    grid.addWidget(self.poly_mode_cb_,row,3,1,2)
    row+=1
    grid.addWidget(spline_label, row, 0, 1, 3)
    grid.addWidget(self.spline_spinbox_, row, 4, 1, 1)
    row+=1
    grid.addWidget(arc_label,row,0,1,3)
    grid.addWidget(self.arc_spinbox_,row,4,1,1)
    row+=1
    grid.addWidget(tube_label, row, 0, 1, 3)
    row+=1
    grid.addWidget(radius_tube_label, row, 0, 1, 1)
    grid.addWidget(self.width_tube_slider_, row, 1, 1, 3)
    grid.addWidget(self.width_tube_spinbox_, row, 4, 1, 1)
    row+=1
    grid.addWidget(ratio_tube_label, row, 0, 1, 1)
    grid.addWidget(self.thickness_tube_slider_, row, 1, 1, 3)
    grid.addWidget(self.thickness_tube_spinbox_, row, 4, 1, 1)
    row+=1
    grid.addWidget(tube_profile_label, row, 0, 1, 3)
    grid.addWidget(self.tube_profile_spinbox_, row, 4, 1, 1)
    row+=1
    grid.addWidget(helix_label, row, 0, 1, 3)
    row+=1
    grid.addWidget(radius_helix_label, row, 0, 1, 1)
    grid.addWidget(self.width_helix_slider_, row, 1, 1, 3)
    grid.addWidget(self.width_helix_spinbox_, row, 4, 1, 1)
    row+=1
    grid.addWidget(ratio_helix_label, row, 0, 1, 1)
    grid.addWidget(self.thickness_helix_slider_, row, 1, 1, 3)
    grid.addWidget(self.thickness_helix_spinbox_, row, 4, 1, 1)
    row+=1
    grid.addWidget(ecc_helix_label, row, 0, 1, 1)
    grid.addWidget(self.ecc_helix_slider_, row, 1, 1, 3)
    grid.addWidget(self.ecc_helix_spinbox_, row, 4, 1, 1)
    row+=1
    grid.addWidget(helix_profile_label, row, 0, 1, 3)
    grid.addWidget(self.helix_profile_spinbox_, row, 4, 1, 1)
    row+=1
    grid.addWidget(helix_mode_label, row, 0, 1, 3)
    grid.addWidget(self.helix_mode_spinbox_, row, 4, 1, 1)
    row+=1
    grid.addWidget(strand_label, row, 0, 1, 3)
    row+=1
    grid.addWidget(radius_strand_label, row, 0, 1, 1)
    grid.addWidget(self.width_strand_slider_, row, 1, 1, 3)
    grid.addWidget(self.width_strand_spinbox_, row, 4, 1, 1)
    row+=1
    grid.addWidget(ratio_strand_label, row, 0, 1, 1)
    grid.addWidget(self.thickness_strand_slider_, row, 1, 1, 3)
    grid.addWidget(self.thickness_strand_spinbox_, row, 4, 1, 1)
    row+=1
    grid.addWidget(ecc_strand_label, row, 0, 1, 1)
    grid.addWidget(self.ecc_strand_slider_, row, 1, 1, 3)
    grid.addWidget(self.ecc_strand_spinbox_, row, 4, 1, 1)
    row+=1
    grid.addWidget(strand_profile_label, row, 0, 1, 3)
    grid.addWidget(self.strand_profile_spinbox_, row, 4, 1, 1)
    row+=1
    grid.addWidget(strand_mode_label, row, 0, 1, 3)
    grid.addWidget(self.strand_mode_spinbox_, row, 4, 1, 1)
    
    grid.setRowStretch(row+1,1)
    self.setLayout(grid)
    
    self.spline_spinbox_.valueChanged.connect(self.UpdateSplineDetail)
    self.arc_spinbox_.valueChanged.connect(self.UpdateArcDetail)
    self.poly_mode_cb_.currentIndexChanged.connect(self.UpdatePolyMode)
    
    self.width_tube_spinbox_.valueChanged.connect(self.UpdateTubeRadius)
    self.width_tube_slider_.valueChanged.connect(self.UpdateSliderTubeRadius)
    self.thickness_tube_spinbox_.valueChanged.connect(self.UpdateTubeRatio)
    self.thickness_tube_slider_.valueChanged.connect(self.UpdateSliderTubeRatio)
    self.tube_profile_spinbox_.valueChanged.connect(self.UpdateTubeProfileType)
    
    self.width_helix_spinbox_.valueChanged.connect(self.UpdateHelixWidth)
    self.width_helix_slider_.valueChanged.connect(self.UpdateSliderHelixWidth)
    self.thickness_helix_spinbox_.valueChanged.connect(self.UpdateHelixThickness)
    self.thickness_helix_slider_.valueChanged.connect(self.UpdateSliderHelixThickness)
    self.ecc_helix_spinbox_.valueChanged.connect(self.UpdateHelixEcc)
    self.ecc_helix_slider_.valueChanged.connect(self.UpdateSliderHelixEcc)
    self.helix_profile_spinbox_.valueChanged.connect(self.UpdateHelixProfileType)
    self.helix_mode_spinbox_.valueChanged.connect(self.UpdateHelixMode)

    self.width_strand_spinbox_.valueChanged.connect(self.UpdateStrandWidth)
    self.width_strand_slider_.valueChanged.connect(self.UpdateSliderStrandWidth)
    self.thickness_strand_spinbox_.valueChanged.connect(self.UpdateStrandThickness)
    self.thickness_strand_slider_.valueChanged.connect(self.UpdateSliderStrandThickness)
    self.ecc_strand_spinbox_.valueChanged.connect(self.UpdateStrandEcc)
    self.ecc_strand_slider_.valueChanged.connect(self.UpdateSliderStrandEcc)    
    self.strand_profile_spinbox_.valueChanged.connect(self.UpdateStrandProfileType)
    self.strand_mode_spinbox_.valueChanged.connect(self.UpdateStrandMode)
    
    self.setMinimumSize(250,550)
    ########/UI########
  
  def UpdateGui(self,options):
    self.poly_mode_cb_.setCurrentIndex(options.GetPolyMode())
    self.spline_spinbox_.setValue(options.GetSplineDetail())
    self.arc_spinbox_.setValue(options.GetArcDetail())
    
    self.UpdateTubeRadiusGui(options.GetTubeRadius())
    self.UpdateTubeRatioGui(options.GetTubeRatio())
    self.tube_profile_spinbox_.setValue(options.GetTubeProfileType())
    self.UpdateHelixWidthGui(options.GetHelixWidth())
    self.UpdateHelixThicknessGui(options.GetHelixThickness())
    self.UpdateHelixEccGui(options.GetHelixEcc())
    self.helix_profile_spinbox_.setValue(options.GetHelixProfileType())
    self.helix_mode_spinbox_.setValue(options.GetHelixMode())
    self.UpdateStrandWidthGui(options.GetStrandWidth())
    self.UpdateStrandThicknessGui(options.GetStrandThickness())
    self.UpdateStrandEccGui(options.GetStrandEcc())
    self.strand_profile_spinbox_.setValue(options.GetStrandProfileType())
    self.strand_mode_spinbox_.setValue(options.GetStrandMode())
    
  def UpdatePolyMode(self, value):
    self.GetOptions().SetPolyMode(value)
    self.ApplyOptions()
    
  def UpdateSplineDetail(self, value):
    self.GetOptions().SetSplineDetail(value)
    self.ApplyOptions()
    
  def UpdateArcDetail(self, value):
    self.GetOptions().SetArcDetail(value)
    self.ApplyOptions()
    
  def UpdateTubeRadius(self, value):
    self.GetOptions().SetTubeRadius(value)
    self.ApplyOptions()
    
  def UpdateSliderTubeRadius(self, value):
    self.GetOptions().SetTubeRadius(value/10.0)
    self.ApplyOptions()
    
  def UpdateTubeRatio(self, value):
    self.GetOptions().SetTubeRatio(value)
    self.ApplyOptions()
    
  def UpdateSliderTubeRatio(self, value):
    self.GetOptions().SetTubeRatio(value/10.0)
    self.ApplyOptions()
    
  def UpdateTubeProfileType(self, value):
    self.GetOptions().SetTubeProfileType(value)
    self.ApplyOptions()
    
  def UpdateHelixWidth(self, value):
    self.GetOptions().SetHelixWidth(value) 
    self.ApplyOptions()
    
  def UpdateSliderHelixWidth(self, value):
    self.GetOptions().SetHelixWidth(value/10.0)
    self.ApplyOptions()
    
  def UpdateHelixThickness(self, value):
    self.GetOptions().SetHelixThickness(value)
    self.ApplyOptions()
    
  def UpdateSliderHelixThickness(self, value):
    self.GetOptions().SetHelixThickness(value/10.0)
    self.ApplyOptions()
    
  def UpdateHelixEcc(self, value):
    self.GetOptions().SetHelixEcc(value)
    self.ApplyOptions()
    
  def UpdateHelixProfileType(self, value):
    self.GetOptions().SetHelixProfileType(value)
    self.ApplyOptions()

  def UpdateHelixMode(self, value):
    self.GetOptions().SetHelixMode(value)
    self.ApplyOptions()
    
  def UpdateSliderHelixEcc(self, value):
    self.GetOptions().SetHelixEcc(value/10.0)
    self.ApplyOptions()
    
  def UpdateStrandWidth(self, value):
    self.GetOptions().SetStrandWidth(value)
    self.ApplyOptions()
    
  def UpdateSliderStrandWidth(self, value):
    self.GetOptions().SetStrandWidth(value/10.0)
    self.ApplyOptions()
    
  def UpdateStrandThickness(self, value):
    self.GetOptions().SetStrandThickness(value)
    self.ApplyOptions()
    
  def UpdateSliderStrandThickness(self, value):
    self.GetOptions().SetStrandThickness(value/10.0)
    self.ApplyOptions()
    
  def UpdateStrandEcc(self, value):
    self.GetOptions().SetStrandEcc(value)
    self.ApplyOptions()
    
  def UpdateSliderStrandEcc(self, value):
    self.GetOptions().SetStrandEcc(value/10.0)
    self.ApplyOptions()
    
  def UpdateStrandProfileType(self, value):
    self.GetOptions().SetStrandProfileType(value)
    self.ApplyOptions()
    
  def UpdateStrandMode(self, value):
    self.GetOptions().SetStrandMode(value)
    self.ApplyOptions()
    
  def UpdateTubeRadiusGui(self,value):
    value = round(value, 2)
    if(abs(value*10.0 - self.width_tube_slider_.value())>=self.width_tube_slider_.singleStep()):
      self.width_tube_slider_.setValue(int(value*10.0))
    if(abs(value - self.width_tube_spinbox_.value())>=self.width_tube_spinbox_.singleStep()):
      self.width_tube_spinbox_.setValue(value)

  def UpdateTubeRatioGui(self,value):
    value = round(value, 2)
    if(abs(value*10.0 - self.thickness_tube_slider_.value())>=self.thickness_tube_slider_.singleStep()):
      self.thickness_tube_slider_.setValue(int(value*10.0))
    if(abs(value - self.thickness_tube_spinbox_.value())>=self.thickness_tube_spinbox_.singleStep()):
      self.thickness_tube_spinbox_.setValue(value)

  def UpdateHelixWidthGui(self, value):
    value = round(value, 2)
    if(abs(value*10.0 - self.width_helix_slider_.value())>=self.width_helix_slider_.singleStep()):
      self.width_helix_slider_.setValue(int(value*10.0))
    if(abs(value - self.width_helix_spinbox_.value())>=self.width_helix_spinbox_.singleStep()):
      self.width_helix_spinbox_.setValue(value)
    
  def UpdateHelixThicknessGui(self, value):
    value = round(value, 2)
    if(abs(value*10.0 - self.thickness_helix_slider_.value())>=self.thickness_helix_slider_.singleStep()):
      self.thickness_helix_slider_.setValue(int(value*10.0))
    if(abs(value - self.thickness_helix_spinbox_.value())>=self.thickness_helix_spinbox_.singleStep()):
      self.thickness_helix_spinbox_.setValue(value)
    
  def UpdateHelixEccGui(self, value):
    value = round(value, 2)
    if(abs(value*10.0 - self.ecc_helix_slider_.value())>=self.ecc_helix_slider_.singleStep()):
      self.ecc_helix_slider_.setValue(int(value*10.0))
    if(abs(value - self.ecc_helix_spinbox_.value())>=self.ecc_helix_spinbox_.singleStep()):
      self.ecc_helix_spinbox_.setValue(value)
    
  def UpdateStrandWidthGui(self, value):
    value = round(value, 2)
    if(abs(value*10.0 - self.width_strand_slider_.value())>=self.width_strand_slider_.singleStep()):
      self.width_strand_slider_.setValue(int(value*10.0))
    if(abs(value - self.width_strand_spinbox_.value())>=self.width_strand_spinbox_.singleStep()):
      self.width_strand_spinbox_.setValue(value)

  def UpdateStrandThicknessGui(self, value):
    value = round(value, 2)
    if(abs(value*10.0 - self.thickness_strand_slider_.value())>=self.thickness_strand_slider_.singleStep()):
      self.thickness_strand_slider_.setValue(int(value*10.0))
    if(abs(value - self.thickness_strand_spinbox_.value())>=self.thickness_strand_spinbox_.singleStep()):
      self.thickness_strand_spinbox_.setValue(value)
    
  def UpdateStrandEccGui(self, value):
    value = round(value, 2)
    if(abs(value*10.0 - self.ecc_strand_slider_.value())>=self.ecc_strand_slider_.singleStep()):
      self.ecc_strand_slider_.setValue(int(value*10.0))
    if(abs(value - self.ecc_strand_spinbox_.value())>=self.ecc_strand_spinbox_.singleStep()):
      self.ecc_strand_spinbox_.setValue(value)

  def GetText(self):
    return self.text_

  def GetRenderMode(self):
    return self.mode_
