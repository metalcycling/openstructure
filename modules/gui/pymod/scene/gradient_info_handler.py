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
from PyQt5 import QtCore
from .immutable_gradient_info_handler import ImmutableGradientInfoHandler

#Gradient Info Handler
class GradientInfoHandler(ImmutableGradientInfoHandler):
  def __init__(self, fileName):    
    ImmutableGradientInfoHandler.__init__(self, fileName)
  
  def StoreQGradient(self,gradient,name):
    group = self.gradients_.CreateGroup(self.GRADIENT_GROUP_NAME)
    group.SetAttribute(self.NAME_ATTRIBUTE_NAME, name)
    gfx_gradient = self.ConvertToGfxGradient(gradient)
    
    gfx_gradient.GradientToInfo(group)
    self.handle_.Export(self.FILE_NAME)
    
  def RemoveGradient(self, name):
    group_list = self.gradients_.GetGroups(self.GRADIENT_GROUP_NAME)
    for group in group_list:
      if group.HasAttribute(self.NAME_ATTRIBUTE_NAME):
        groupname = group.GetAttribute(self.NAME_ATTRIBUTE_NAME)
        if name == groupname:
          groupToDel = group
          break
    
    if groupToDel != None:
      self.gradients_.Remove(groupToDel)      
      self.handle_.Export(self.FILE_NAME)
      
  def RenameGradient(self, old, new):
    group_list = self.gradients_.GetGroups(self.GRADIENT_GROUP_NAME)
    for group in group_list:
      if group.HasAttribute(self.NAME_ATTRIBUTE_NAME):
        groupname = group.GetAttribute(self.NAME_ATTRIBUTE_NAME)
        if old == groupname:
          group.SetAttribute(self.NAME_ATTRIBUTE_NAME, new)
          self.handle_.Export(self.FILE_NAME)
