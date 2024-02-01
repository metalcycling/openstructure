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
from .immutable_preset_info_handler import ImmutablePresetInfoHandler

#Preset Info Handler
class PresetInfoHandler(ImmutablePresetInfoHandler):
  def __init__(self, file_name):    
    ImmutablePresetInfoHandler.__init__(self, file_name)
  
  def StorePreset(self,preset):
    group = self.presets_.CreateGroup(ImmutablePresetInfoHandler.PRESET_GROUP_NAME)    
    preset.ToInfo(group)
    self.handle_.Export(self.file_name_)
    
  def RemovePreset(self, name):
    group_list = self.presets_.GetGroups(ImmutablePresetInfoHandler.PRESET_GROUP_NAME)
    group_to_del = None
    for group in group_list:
      if group.HasAttribute(ImmutablePresetInfoHandler.NAME_ATTRIBUTE_NAME):
        groupname = group.GetAttribute(ImmutablePresetInfoHandler.NAME_ATTRIBUTE_NAME)
        if name == groupname:
          group_to_del = group
          break
    
    if group_to_del != None:
      self.presets_.Remove(group_to_del)      
      self.handle_.Export(self.file_name_)
      
  def RenamePreset(self, old, new):
    group_list = self.presets_.GetGroups(ImmutablePresetInfoHandler.PRESET_GROUP_NAME)
    for group in group_list:
      if group.HasAttribute(ImmutablePresetInfoHandler.NAME_ATTRIBUTE_NAME):
        groupname = group.GetAttribute(ImmutablePresetInfoHandler.NAME_ATTRIBUTE_NAME)
        if old == groupname:
          group.SetAttribute(ImmutablePresetInfoHandler.NAME_ATTRIBUTE_NAME, new)
          self.handle_.Export(self.file_name_)
