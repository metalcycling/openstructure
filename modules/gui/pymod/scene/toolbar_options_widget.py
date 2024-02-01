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
from PyQt5 import QtCore, QtWidgets, QtGui

class ToolBarOptionsWidget(QtWidgets.QWidget):
  """QWidget with a ToolBar and a show area.
    
   This abstract QWidget has a toolbar and a show area. Whenever a button of the tool bar is pressed, the Widget corresponding to the button's value is shown in the show area.  
  """
  def __init__(self, parent=None):
    QtWidgets.QWidget.__init__(self, parent)
    
    #Setup ui_
    self.parent_=parent
    self.resize(400, 300)
    self.setMinimumSize(QtCore.QSize(250, 200))
    self.gridLayout = QtWidgets.QGridLayout(self)
    self.gridLayout.setHorizontalSpacing(0)
    self.gridLayout.setVerticalSpacing(0)
    self.gridLayout.setContentsMargins(0,0,0,0)
    self.gridLayout.setSpacing(0)    
    self.tool_bar_ = QtWidgets.QToolBar(self)
    self.tool_bar_.setIconSize(QtCore.QSize(16, 16))
    self.gridLayout.addWidget(self.tool_bar_, 0, 0, 1, 1)
    self.stackedWidget = QtWidgets.QStackedWidget(self)
    self.gridLayout.addWidget(self.stackedWidget, 1, 0, 1, 1)

    self.current_action_ = None
    self.actions_ = list()
    
    self.tool_bar_.actionTriggered.connect(self.ChangeSelectedItem)
    
    self.setEnabled(False)
           
    self.Update()
  
  def Update(self):
    """Updates the active widget of the show area.
    
     This method calls the Update method of the active widget.
    """
    self.setEnabled(True)
    widget = self.__GetCurrentWidget()
    if hasattr(widget, "Update"):
      widget.Update()
  
    
  def AddWidget(self, ident, widget, text=None):
    """Adds a Widget to this Options Widget.
    
     The Widget must have a identifier. If another Widget has the same identifier,
     the old widget will be removed and the new widget gets the identifier.
     Returns True, if widget is added. Otherwise it returns False
    """
    if isinstance(widget, QtWidgets.QWidget) and ident is not None:
      if text is not None:
        string = text
      elif hasattr(widget, "GetText"):
        string = widget.GetText()
      else:
        string = ident
      
      self.stackedWidget.addWidget(widget)
      action = self.tool_bar_.addAction(ident)
      action.setIcon(QtGui.QIcon(ident))
      action.setToolTip(string)
      pair = ident, widget
      action.setData(QtCore.QVariant(pair))
      action.setCheckable(True);
      if(len(self.actions_) == 0):
        self.ChangeSelectedItem(action)
      self.actions_.append(action)
      return True
    return False
  
  def OnComboChange(self, item):
    """This abstract method is called whenever the View is updated.
    
     This abstract method must be implemented by all subclasses. 
     It can be used to do something ;-) whenever the combobox changes its value.
    """
    raise NotImplementedError("Subclasses must define OnComboChange()")
  
  def DoResize(self):
    item = self.__GetCurrentWidget()
    width = 0
    height = 0
    if(hasattr(item,"minimumHeight")):
      height=item.minimumHeight()
    if(hasattr(item,"minimumWidth")):
      width=item.minimumWidth()
    self.setMinimumSize(width,self.tool_bar_.height()+height)
    if(hasattr(self.parent_,"DoResize")):
      self.parent_.DoResize()
  
  def ChangeSelectedItem(self, action):
    """Change Current Selected Item.
    
     Shows the widget which corresponds to the action in the show area.
    """
    if(self.current_action_ != None):
        self.current_action_.setChecked(False)
    else:
       self.current_action_ = action
    widget = action.data()[1]
    self.stackedWidget.setCurrentWidget(widget)
    if hasattr(widget, "Update"): 
      widget.Update()
    if(self.current_action_ == action):
      self.current_action_.setChecked(True)
    else:
      self.current_action_=action
    self.OnComboChange(widget)
  #Private Methods
  def __GetCurrentWidget(self):
      return self.stackedWidget.currentWidget();

  #Overwritten Methods
  def setEnabled(self, bool):
    QtWidgets.QWidget.setEnabled(self, bool)

