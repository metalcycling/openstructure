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


#Tube Render Options
class RenderModeWidget(QtWidgets.QWidget):
  def __init__(self, parent=None):
    QtWidgets.QWidget.__init__(self, parent)
    self.options_ = None
    self.entities_ = set()

  def GetText(self):
    raise NotImplementedError("Subclasses must define GetText()")
  
  def GetRenderMode(self):
    raise NotImplementedError("Subclasses must define GetRenderMode()")
    
  def UpdateGui(self, options):
    pass
  
  def Update(self):
    self.entities_.clear()
    
    scene_selection = gui.SceneSelection.Instance()
    if scene_selection.GetActiveNodeCount() == 0 and scene_selection.GetActiveViewCount() == 0:
      self.setEnabled(False)
      return
    
    if scene_selection.GetActiveNodeCount() > 0 :
      for i in range(0,scene_selection.GetActiveNodeCount()):
        entity = scene_selection.GetActiveNode(i)
        if isinstance(entity, gfx.Entity):
          self.entities_.add(entity)
        else:
          self.setEnabled(False)
          return

    if scene_selection.GetActiveViewCount() > 0 :
      entity = scene_selection.GetViewEntity()
      self.entities_.add(entity)
    
    if len(self.entities_)>0:
      entity = self.entities_.pop()
      self.options_=entity.GetOptions(self.GetRenderMode())
      self.UpdateGui(self.options_)
      QtWidgets.QWidget.setEnabled(self,True)
    else:
      QtWidgets.QWidget.setEnabled(self,False)
    
  def ApplyOptions(self):
    for entity in self.entities_:
      entity.ApplyOptions(self.GetRenderMode(), self.GetOptions())

  def GetOptions(self):
    return self.options_
