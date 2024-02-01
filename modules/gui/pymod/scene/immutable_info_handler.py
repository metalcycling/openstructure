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
from ost import info
from PyQt5 import QtCore, QtWidgets

#Gradient Editor
class ImmutableInfoHandler:
  def __init__(self, fileName):
    
    self.FILE_NAME = fileName    
    self.GRADIENTS_GROUP_NAME = "Gradients"
    self.GRADIENT_GROUP_NAME = "Gradient"
    self.NAME_ATTRIBUTE_NAME = "Name"
    
    self.handle_ = info.LoadOrCreateInfo(self.FILE_NAME)
        
    self.root_ = self.handle_.Root()
    
    self.gradients_ = self.root_.RetrieveGroup(self.GRADIENTS_GROUP_NAME)

    
  def GetGfxGradient(self, name):
    group = self.gradients_.GetGroup(name)
    gfxGradient = gfx.Gradient.GradientFromInfo(group)
    return gfxGradient
  
  def GetQGradient(self, name):
    gfxGradient = self.GetGfxGradient(name)
    return self.ConvertToQGradient(gfxGradient)
  
  def ConvertToQGradient(self, gfxGradient):
    stops = gfxGradient.GetStops()
    qgrad = QtWidgets.QLinearGradient()
    for stop in stops:
      color = stop.GetColor()
      qcolor = QtWidgets.QColor(color.Red()*255,
                                color.Green()*255,
                                color.Blue()*255,
                                color.Alpha()*255)
      qgrad.setColorAt(stop.GetRel(), qcolor) 
    return qgrad 
  
  def ConvertToGfxGradient(self, gradient):  
    gfxgradient = gfx.Gradient()    
    for s in gradient.stops():
      rel=s[0]
      color=s[1]
      gfxgradient.SetColorAt(s[0], gfx.RGB(s[1].redF(), s[1].greenF(), s[1].blueF()));
    return gfxgradient;
  
  def GetQGradients(self):
    groupList = self.gradients_.GetGroups(self.GRADIENT_GROUP_NAME)
    qGradients = dict()
    
    for group in groupList:
      if group.HasAttribute(self.NAME_ATTRIBUTE_NAME):
        gfxGradient = gfx.Gradient.GradientFromInfo(group)
        qGradients[group.GetAttribute(self.NAME_ATTRIBUTE_NAME)] = (self.ConvertToQGradient(gfxGradient))
      
    return qGradients
