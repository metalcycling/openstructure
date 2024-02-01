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

class RenderOp:
  RENDERMODE_ATTRIBUTE_NAME = "RenderMode"
  KEEP_ATTRIBUTE_NAME = "Keep"
  FLAGS_ATTRIBUTE_NAME = "Flags"
    
  def __init__(self, render_mode, selection, flags, keep=False):
    self.render_mode_ = render_mode
    self.selection_ = selection
    self.keep_ = keep
    self.flags_ = flags

  def GetName(self):
    return "Render mode: %s"%str(self.GetRenderMode())

  def SetRenderMode(self, render_mode):
    self.render_mode_ = render_mode
    
  def GetRenderMode(self):
    return self.render_mode_
 
  def SetSelection(self, selection):
    self.selection_ = selection
    
  def GetSelection(self):
    return self.selection_

  def SetSelectionFlags(self, flags):
    self.flags_ = flags
    
  def GetSelectionFlags(self):
    return self.flags_

  def SetKeep(self, keep):
    self.keep_ = keep

  def IsKept(self):
    return self.keep_

  def ApplyOn(self, entity):
    if (entity is not None) and isinstance(entity, gfx.Entity):
      entity.SetRenderMode(self.GetRenderMode(),entity.view.Select(self.GetSelection(),self.GetSelectionFlags()),self.IsKept())
  
  def ToInfo(self,group):
    group.SetAttribute(RenderOp.RENDERMODE_ATTRIBUTE_NAME, str(self.GetRenderMode().name))
    group.SetAttribute(RenderOp.KEEP_ATTRIBUTE_NAME, str(int(self.IsKept())))
    group.SetAttribute(RenderOp.FLAGS_ATTRIBUTE_NAME, str(self.GetSelectionFlags()))
    group.SetTextData(str(self.GetSelection()))
    
  @staticmethod
  def FromInfo(group):
    render_op = None
    if (group.HasAttribute(RenderOp.RENDERMODE_ATTRIBUTE_NAME) 
    and  group.HasAttribute(RenderOp.KEEP_ATTRIBUTE_NAME)):
      render_mode = getattr(gfx.RenderMode, group.GetAttribute(RenderOp.RENDERMODE_ATTRIBUTE_NAME))
      keep = bool(int(group.GetAttribute(RenderOp.KEEP_ATTRIBUTE_NAME)))
      flags=0
      if group.HasAttribute(RenderOp.FLAGS_ATTRIBUTE_NAME):
        flags = int(group.GetAttribute(RenderOp.FLAGS_ATTRIBUTE_NAME))
      selection = group.GetTextData()
      render_op = RenderOp(render_mode,selection,flags,keep)
    return render_op
