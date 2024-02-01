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
try: 
  from ost import img
  from .wireframe_widget import WireframeWidget
  _img_present=True
except ImportError:
  _img_present=False
  pass
from PyQt5 import QtCore, QtWidgets
from .scene_selection_helper import SelHelper
from .combo_options_widget import ComboOptionsWidget
from .custom_widget import CustomWidget
from .cpk_widget import CPKWidget
from .tube_widget import TubeWidget
from .simple_widget import SimpleWidget
from .sline_widget import SlineWidget
from .hsc_widget import HSCWidget
from .trace_widget import TraceWidget
from .line_trace_widget import LineTraceWidget

class RenderOptionsWidget(ComboOptionsWidget):
  def __init__(self, parent=None):
    ComboOptionsWidget.__init__(self, parent)
    
    #Title
    self.text_ = "Render Options"
    
    self.keep_action_ = QtWidgets.QAction("K",self)
    self.keep_action_.setCheckable(True);
    self.keep_action_.setChecked(False)
    self.keep_action_.setToolTip("Keep rendermodes and add current")
    self.keep_button_ = QtWidgets.QToolButton(self)
    self.keep_button_.setDefaultAction(self.keep_action_)
    
    self.grid_layout_.addWidget(self.keep_button_, 0, 1, 1, 1)
    self.grid_layout_.addWidget(self.stacked_widget_, 1, 0, 1, 2)
    
    #Add options to menu    
    self.entity_widgets_ = list()
    self.entity_widgets_.append(["", EmptyMode()])
    self.entity_widgets_.append([gfx.RenderMode.SIMPLE, SimpleWidget()])
    self.entity_widgets_.append([gfx.RenderMode.CUSTOM, CustomWidget()])
    self.entity_widgets_.append([gfx.RenderMode.CPK, CPKWidget()])
    self.entity_widgets_.append([gfx.RenderMode.LINE_TRACE, LineTraceWidget()])
    self.entity_widgets_.append([gfx.RenderMode.TRACE, TraceWidget()])
    self.entity_widgets_.append([gfx.RenderMode.SLINE, SlineWidget()])
    self.entity_widgets_.append([gfx.RenderMode.TUBE, TubeWidget()])
    self.entity_widgets_.append([gfx.RenderMode.HSC, HSCWidget()])

    self.img_widgets_ = list()
    if _img_present:
      self.img_widgets_.append(["", EmptyMode()])
      self.img_widgets_.append([gfx.RenderMode.SIMPLE, WireframeWidget()])
      self.img_widgets_.append([gfx.RenderMode.FILL, EmptyMode("Fill",gfx.RenderMode.FILL)])

    self._in_view_method = False
    self.setMinimumSize(250,200)
  def OnActivate(self, item):
    self.DoResize()
  def OnComboChange(self, item):
    scene_selection = gui.SceneSelection.Instance()
    for i in range(0,scene_selection.GetActiveNodeCount()):
      node = scene_selection.GetActiveNode(i)
      if isinstance(node, gfx.Entity) or (_img_present and isinstance(node, gfx.MapIso)):
        render_mode = item.GetRenderMode()
        if render_mode is not None:
          node.SetRenderMode(render_mode)
    if(scene_selection.GetActiveViewCount() > 0):
      entity = scene_selection.GetViewEntity()
      view = scene_selection.GetViewUnion()        
      render_mode = item.GetRenderMode()
      if render_mode is not None:
        entity.SetRenderMode(item.GetRenderMode(),view,self.keep_action_.isChecked())
        
    item.Update()
    self.DoResize()
  
  def Update(self):
    if hasattr(self, "keep_button_"):
      self.keep_button_.setEnabled(True)
    
    ComboOptionsWidget.setEnabled(self,True)
    
    cur_widget = self.GetCurrentWidget()
    new_render_mode = None
    if cur_widget is not None:
      new_render_mode = cur_widget.GetRenderMode()
    
    if SelHelper().CheckAllFlags(SelHelper.NO_SELECTION):
      ComboOptionsWidget.setEnabled(self,False)
      return
    
    if not self._in_view_method:
      for w in self.entity_widgets_:
        self.RemoveWidget(w[0])
      for w in self.img_widgets_:
        self.RemoveWidget(w[0])

    scene_selection = gui.SceneSelection.Instance()        
    if scene_selection.GetActiveNodeCount() > 0 :
      if hasattr(self, "keep_button_"):
        self.keep_button_.setEnabled(False)
      render_mode_valid = True
      render_mode = None
      for i in range(0,scene_selection.GetActiveNodeCount()):
        node = scene_selection.GetActiveNode(i)
        if isinstance(scene_selection.GetActiveNode(i), gfx.GfxObj):
          if render_mode is None:
            render_mode = node.GetRenderMode()
          elif render_mode != node.GetRenderMode():
            render_mode_valid = False
 
    if SelHelper().CheckFlags(SelHelper.HAS_IMG | SelHelper.IS_ONE_TYPE):
      for w in self.img_widgets_:
        self.AddWidget(w[0], w[1])
    elif SelHelper().CheckMinOneFlag(SelHelper.HAS_ENTITY| SelHelper.HAS_VIEW) and SelHelper().CheckNotFlags(SelHelper.HAS_IMG):
      if not self._in_view_method:
        for w in self.entity_widgets_:
          self.AddWidget(w[0], w[1])
    else:
      ComboOptionsWidget.setEnabled(self,False)
      return
    
    if SelHelper().CheckMinOneFlag(SelHelper.HAS_ENTITY| SelHelper.HAS_IMG) and SelHelper().CheckNotFlags(SelHelper.HAS_VIEW):
      if(render_mode_valid):
        ComboOptionsWidget.ChangeSelectedItem(self,render_mode)
      else:
        ComboOptionsWidget.ChangeSelectedItem(self,"")
    else:
      if not self._in_view_method:
        self._in_view_method = True
      ComboOptionsWidget.ChangeSelectedItem(self,new_render_mode)
    self.GetCurrentWidget().Update()
      
  def GetText(self):
    return self.text_
        
class EmptyMode(QtWidgets.QWidget):
  def __init__(self, text="", render_mode=None, parent=None):
    QtWidgets.QLabel.__init__(self, parent)
    self.setMinimumSize(250,30)
    self.text = text
    self.render_mode = render_mode
    if(render_mode):
      text_label = QtWidgets.QLabel(text)
      font = text_label.font()
      font.setBold(True)
      grid = QtWidgets.QGridLayout()
      grid.addWidget(text_label,0,0,1,1)
      grid.addWidget(QtWidgets.QLabel("No Settings available"), 1, 0, 1, 3)
      grid.setRowStretch(2,1)
      self.setLayout(grid)
      self.setMinimumSize(250,60)
      
  def Update(self):
   True #Do Nothing

  def GetText(self):
    return self.text

  def GetRenderMode(self):
    return self.render_mode
