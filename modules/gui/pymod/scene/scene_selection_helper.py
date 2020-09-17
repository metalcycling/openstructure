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

import sys
from ost import gui
from ost import gfx
import ost
import os
try: 
  from ost import img
  _img_present=True
except ImportError:
  _img_present=False
  pass

class SelHelper:
  __shared_state = {}
  
  NO_SELECTION = 0
  HAS_ENTITY = 1
  HAS_VIEW = 2
  HAS_IMG = 4
  HAS_SURFACE = 8
  IS_ONE_TYPE = 16
  IS_MULTI_TYPE = 32
  SINGLE_SELECTION = 64
  MULTI_SELECTION = 128
    
  def __init__(self):
    self.__dict__ = self.__shared_state
    if not '_ready' in dir(self):
      self.scene_sel_ = gui.SceneSelection.Instance()
      self.current_flags_ = 0
      self._ready = True

  def Update(self):
    self.current_flags_ = 0
    if self.scene_sel_.GetActiveNodeCount() == 0 and self.scene_sel_.GetActiveViewCount() == 0:
      return
    
    for i in range(0,self.scene_sel_.GetActiveNodeCount()):
      node = self.scene_sel_.GetActiveNode(i)
      if isinstance(node, gfx.Entity):
        self.current_flags_ = self.current_flags_ | SelHelper.HAS_ENTITY 
      if isinstance(node, gfx.Surface):
        self.current_flags_ = self.current_flags_ | SelHelper.HAS_SURFACE
      if (_img_present) and isinstance(node, gfx.MapIso):
        self.current_flags_ = self.current_flags_ | SelHelper.HAS_IMG
    
    if self.scene_sel_.GetActiveViewCount() > 0:
      self.current_flags_ = self.current_flags_ | SelHelper.HAS_VIEW
    
    cnt = 0
    if self.current_flags_ & SelHelper.HAS_ENTITY:
      cnt += 1
    if self.current_flags_ & SelHelper.HAS_SURFACE:
      cnt += 1
    if self.current_flags_ & SelHelper.HAS_IMG:
      cnt += 1
    if self.current_flags_ & SelHelper.HAS_VIEW:
      cnt += 1
    
    if cnt == 1:
      self.current_flags_ = self.current_flags_ | SelHelper.IS_ONE_TYPE
    elif cnt > 1:
      self.current_flags_ = self.current_flags_ | SelHelper.IS_MULTI_TYPE
    
    if self.scene_sel_.GetActiveNodeCount()==1:
      self.current_flags_ = self.current_flags_ | SelHelper.SINGLE_SELECTION
    elif self.scene_sel_.GetActiveNodeCount()>1:
      self.current_flags_ = self.current_flags_ | SelHelper.MULTI_SELECTION
    
  def CheckAllFlags(self, flags):
    if(flags == self.current_flags_ & flags) and (flags == self.current_flags_ | flags):
      return True
    return False
  
  def CheckNotFlags(self, flags):
    return not self.CheckFlags(flags)
    
  def CheckFlags(self, flags):
    if(flags == self.current_flags_ & flags):
      return True
    return False
  
  def CheckMinOneFlag(self, flags):
    if((self.current_flags_ - (self.current_flags_ & flags)) < self.current_flags_):
      return True
    return False
