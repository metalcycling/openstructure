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
from PyQt5 import QtCore, QtWidgets

class ComboOptionsWidget(QtWidgets.QWidget):
  """QWidget with a Combobox and a show area.
    
   This abstract QWidget has a Combobox and a show area. Whenever the value of 
   the Combobox changes, the Widget corresponding to the Combobox's value is 
   shown in the show area.  
  """
  def __init__(self, parent=None):
    QtWidgets.QWidget.__init__(self, parent)
    #Setup ui_
    self.parent_ = parent
    self.grid_layout_ = QtWidgets.QGridLayout(self)
    self.grid_layout_.setHorizontalSpacing(0)
    self.grid_layout_.setVerticalSpacing(0)
    self.grid_layout_.setContentsMargins(0,0,0,0)
    self.combo_box_ = QtWidgets.QComboBox(self)
    self.grid_layout_.addWidget(self.combo_box_, 0, 0, 1, 1)
    self.stacked_widget_ = QtWidgets.QStackedWidget(self)
    self.grid_layout_.addWidget(self.stacked_widget_, 1, 0, 1, 1)

    self.__UpdateView(self.combo_box_.currentIndex())
    
    self.combo_box_.activated.connect(self.__UpdateView)
    
    self.setEnabled(False)
       
    self.Update()
    
  def Update(self):
    """Updates the ComboOptionsWidget and the active widget of the show area.
    
     This method calls the Update method of the active widget. 
    """
    self.__GetCurrentPair()[1].Update()
    
  def AddWidget(self, ident, widget):
    """Adds a Widget to this Options Widget.
    
     The Widget must have a identifier. If another Widget has the same identifier,
     the old widget will be removed and the new widget gets the identifier.
     Returns True, if widget is added. Otherwise it returns False
    """
    if isinstance(widget, QtWidgets.QWidget) and ident is not None:
      if hasattr(widget, "GetText"):
        string = widget.GetText()
      else:
        string = ident
      
      self.RemoveWidget(ident)
      self.stacked_widget_.addWidget(widget)
      qpair = ident, widget
      self.combo_box_.addItem(string, QtCore.QVariant(qpair))
      return True
    return False
  
  def RemoveWidget(self,ident):
    index = self.__GetIndex(ident)
    if(index >= 0):
      self.stacked_widget_.removeWidget(self.combo_box_.itemData(index)[1])
      self.combo_box_.removeItem(index)
  
  def OnComboChange(self, item):
    """This abstract method is called whenever the View is updated.
    
     This abstract method must be implemented by all subclasses. 
     It can be used to do something ;-) whenever the combobox changes its value.
    """
    raise NotImplementedError("Subclasses must define OnComboChange()")

  def OnActivate(self, item):
    return self.OnComboChange(self, item)
    
  def ChangeSelectedItem(self, ident):
    """
    Change Current Selected Item.
    
    Shows the widget which corresponds to the ident in the show area. If ident 
    is not valid, nothing happens.
    """
    i = self.__GetIndex(ident)
    if(i>=0) and self.combo_box_.currentIndex() != i:
      self.combo_box_.setCurrentIndex(i)
      if (self.combo_box_.count() > 0):
        pair = self.__GetCurrentPair()
        self.stacked_widget_.setCurrentWidget(pair[1])
        self.OnActivate(pair[1])
  
  def GetCurrentWidget(self):
    if(self.combo_box_.currentIndex() >= 0):
      return self.__GetCurrentPair()[1]
    return None
  
  def DoResize(self):
    item = self.GetCurrentWidget()
    width = 0
    height = 0
    if(hasattr(item,"minimumHeight")):
      height=item.minimumHeight()
    if(hasattr(item,"minimumWidth")):
      width=item.minimumWidth()
    self.setMinimumSize(width,40+height)
    if(hasattr(self.parent_,"DoResize")):
      self.parent_.DoResize()
  
  #Private Methods
  def __UpdateView(self, item):
    if (self.combo_box_.count() > 0):
      pair = self.__GetCurrentPair()
      self.stacked_widget_.setCurrentWidget(pair[1])
      self.OnComboChange(pair[1])

  def __GetIndex(self, ident):
    for i in range(self.combo_box_.count()):
      pair = self.combo_box_.itemData(i)
      if ident == pair[0]:
        return i
    return -1
  
  def __GetCurrentPair(self):
      current_index = self.combo_box_.currentIndex()
      return self.combo_box_.itemData(current_index)

  #Overwritten Methods
  def setEnabled(self, bool):
    QtWidgets.QWidget.setEnabled(self, bool)
    for i in range(self.combo_box_.count()):
        pair = self.combo_box_.itemData(i)
        pair[1].setEnabled(bool)
