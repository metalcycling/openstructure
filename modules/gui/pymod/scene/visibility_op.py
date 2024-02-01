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

from ost import info
from ost import gfx

class VisibilityOp:
  VISIBLE_ATTRIBUTE_NAME = "Visible"
  FLAGS_ATTRIBUTE_NAME = "Flags"
    
  def __init__(self, selection, flags, visible=False):    
    self.selection_ = selection
    self.visible_ = visible
    self.flags_ = flags
    
  def GetName(self):
    return "Visible: %s"%str(self.IsVisible())
  
  def SetSelection(self, selection):
    self.selection_ = selection
    
  def GetSelection(self):
    return self.selection_

  def SetVisible(self, visible):
    self.visible_ = visible

  def SetSelectionFlags(self, flags):
    self.flags_ = flags
    
  def GetSelectionFlags(self):
    return self.flags_

  def IsVisible(self):
    return self.visible_

  def ApplyOn(self, entity):
    if (entity is not None) and isinstance(entity, gfx.Entity):
      entity.SetVisible(entity.view.Select(self.GetSelection(),self.GetSelectionFlags()),self.IsVisible())
  
  def ToInfo(self,group):
      group.SetAttribute(VisibilityOp.VISIBLE_ATTRIBUTE_NAME, str(int(self.IsVisible())))
      group.SetAttribute(VisibilityOp.FLAGS_ATTRIBUTE_NAME, str(self.GetSelectionFlags()))
      group.SetTextData(str(self.GetSelection()))
    
  @staticmethod
  def FromInfo(group):
    visible_op = None
    if group.HasAttribute(VisibilityOp.VISIBLE_ATTRIBUTE_NAME):
      visible = bool(int(group.GetAttribute(VisibilityOp.VISIBLE_ATTRIBUTE_NAME)))
      flags = 0
      if group.HasAttribute(VisibilityOp.FLAGS_ATTRIBUTE_NAME):
        flags = int(group.GetAttribute(VisibilityOp.FLAGS_ATTRIBUTE_NAME))
      selection = group.GetTextData()
      visible_op = VisibilityOp(selection,flags,visible)
    return visible_op
