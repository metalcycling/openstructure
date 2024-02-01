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

import sys
from ost import gui
from ost import gfx

class SceneObserverImpl(gfx.SceneObserver):
  def __init__(self):
    gfx.SceneObserver.__init__(self)
    self.obs_ = list()                              
  
  def AttachObserver(self, ob):
    self.obs_.append(ob)
  
  def DettachObserver(self, ob):
    self.obs_.remove(ob)
  
  def NodeAdded(self, node):
    for ob in self.obs_ :
      if hasattr(ob, "NodeAdded"):
        ob.NodeAdded(node)

  def NodeRemoved(self, node):
    for ob in self.obs_ :
      if hasattr(ob, "NodeRemoved"):
        ob.NodeRemoved(node)
  
  def NodeChanged(self, node):
    for ob in self.obs_ :
      if hasattr(ob, "NodeChanged"):
        ob.NodeChanged(node)
    
  def SelectionModeChanged(self, sel_mode):
    for ob in self.obs_ :
      if hasattr(ob, "SelectionModeChanged"):
        ob.SelectionModeChanged(node)    
  
  def SelectionChanged(self, obj, sel):
    for ob in self.obs_ :
      if hasattr(ob, "SelectionChanged"):
        ob.SelectionChanged(node)    

  def NodeTransformed(self, node):
      for ob in self.obs_ :
        if hasattr(ob, "NodeTransformed"):
          ob.NodeTransformed(node) 
    
  def NodeRemovedDefault(self, node):
    for ob in self.obs_ :
      if hasattr(ob, "NodeRemovedDefault"):
        ob.NodeRemovedDefault(node)
    
  def RenderModeChanged(self, node):
    for ob in self.obs_ :
      if hasattr(ob, "RenderModeChanged"):
        ob.RenderModeChanged(node)
