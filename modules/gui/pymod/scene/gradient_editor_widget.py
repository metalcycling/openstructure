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
from ost import mol
from PyQt5 import QtCore, QtGui, QtWidgets
from .color_select_widget import ColorSelectWidget
from .gradient_preset_widget import GradientPresetWidget

#Gradient Editor
class GradientEditor(QtWidgets.QWidget):
  def __init__(self, parent=None):
    QtWidgets.QWidget.__init__(self, parent)
    self.parent_ = parent
    
    #Create Ui elements
    gradient_label = QtWidgets.QLabel("Gradient Editor")
    font = gradient_label.font()
    font.setBold(True)
    
    self.prop_combo_box_ = QtWidgets.QComboBox()
    self.gradient_preview_ = GradientPreview()
    self.gradient_edit_ = GradientEdit(self.gradient_preview_,self)
    self.gradient_preset_ = GradientPresetWidget(self.gradient_edit_)
    
    self.prop_combo_box_.addItem("atom B-factor")
    self.prop_combo_box_.addItem("average residue B-factor")
    self.prop_combo_box_.addItem("X-Coordinate")
    self.prop_combo_box_.addItem("Y-Coordinate")
    self.prop_combo_box_.addItem("Z-Coordinate")
    self.prop_combo_box_.addItem("Residue Number")
    self.prop_combo_box_.addItem("Atom Charge")
    
    self.props=["abfac","rbfac","x","y","z","rnum","acharge"]
    
    grid = QtWidgets.QGridLayout()
    grid.setContentsMargins(0,5,0,0)
    grid.addWidget(gradient_label, 0, 0, 1, 1)
    grid.addWidget(self.prop_combo_box_, 1, 0, 1, 1)
    grid.addWidget(self.gradient_preview_, 2, 0, 1, 1)
    grid.addWidget(self.gradient_edit_, 3, 0, 1, 1)
    grid.setRowStretch(4, 1)
    grid.addWidget(self.gradient_preset_, 5, 0, 1, 1)
    self.setLayout(grid)
    
    self.prop_combo_box_.currentIndexChanged.connect(self.Update)
    self.gradient_edit_.gradientUpdated.connect(self.Update)
    self.gradient_preset_.gradientSelected.connect(self.gradient_edit_.LoadGradient)
    
    self.setMinimumSize(250,300)
        
  def Update(self):
    scene_selection = gui.SceneSelection.Instance()
    for i in range(0,scene_selection.GetActiveNodeCount()):
      node = scene_selection.GetActiveNode(i)
      self.ChangeColor(node)
      
    if(scene_selection.GetActiveViewCount() > 0):
      entity = scene_selection.GetViewEntity()
      view = scene_selection.GetViewUnion()
      self.ChangeViewColor(entity,view)
  
  def ChangeColor(self,node):
    if isinstance(node, gfx.Entity) or isinstance(node, gfx.Surface):
      node.CleanColorOps()
      if self.parent_.GetCarbonsOnly():
        node.ColorBy(self.props[self.prop_combo_box_.currentIndex()],
                     self.gradient_edit_.GetGfxGradient(), "ele=C")
      else:
        node.ColorBy(self.props[self.prop_combo_box_.currentIndex()],
                     self.gradient_edit_.GetGfxGradient())

  def ChangeViewColor(self, entity, view):
    if isinstance(entity, gfx.Entity) and isinstance(view, mol.EntityView):
      if self.parent_.GetCarbonsOnly():
        glco=gfx.GradientLevelColorOp(mol.QueryViewWrapper(mol.Query("ele=C"),view),self.props[self.prop_combo_box_.currentIndex()],self.gradient_edit_.GetGfxGradient(),mol.Prop.Level.UNSPECIFIED)
      else:
        glco=gfx.GradientLevelColorOp(mol.QueryViewWrapper(view),self.props[self.prop_combo_box_.currentIndex()],self.gradient_edit_.GetGfxGradient(),mol.Prop.Level.UNSPECIFIED)
      entity.Apply(glco)
      
#Gradient Preview
class GradientPreview(QtWidgets.QWidget):
  def __init__(self, parent=None):
    QtWidgets.QWidget.__init__(self, parent)
    
    #Defaults
    self.border_offset_ = 3
    self.preview_height_ = 25
    QtWidgets.QWidget.__init__(self, parent)
    
    self.gradient_ = QtGui.QLinearGradient()
    
    #Ui
    self.setMinimumSize(0, self.preview_height_ + 4)

  def SetGradient(self, gradient):
    self.gradient_ = gradient
    self.update()

  def paintEvent(self, event):   
    if self.isEnabled() and self.gradient_ is not None:
      size = self.size()
      paint = QtGui.QPainter()
      if paint.begin(self):
        brush = QtGui.QBrush(self.gradient_)
        paint.setBrush(brush)
        paint.drawRect(self.border_offset_,
                       self.border_offset_,
                       size.width() - 2 * self.border_offset_,
                       self.preview_height_)
        paint.end()
    
#Gradient Edit
class GradientEdit(QtWidgets.QWidget):
  
  gradientUpdated = QtCore.pyqtSignal()

  def __init__(self, gradient_preview, parent=None):
    QtWidgets.QWidget.__init__(self, parent)

    #Defaults
    self.gradient_preview_ = gradient_preview
    self.border_offset_ = self.gradient_preview_.border_offset_
    self.edit_height_ = 20
        
    #Stops
    self.stops = list()
    
    
    #ContextMenu    
    self.add_ = QtWidgets.QAction("Add", self)
    self.add_.triggered.connect(self.Add)
    self.addAction(self.add_)
    self.setContextMenuPolicy(QtCore.Qt.ActionsContextMenu)
    
    #Ui
    pal = QtGui.QPalette(self.palette())
    pal.setColor(QtGui.QPalette.Window, QtGui.QColor(255, 255, 255, 128))
    self.setPalette(pal)
    self.setAutoFillBackground(True) 
    self.gradient_ = QtGui.QLinearGradient()
    
    self.setMinimumSize(0, self.edit_height_)
    
    #Add Points
    self.AddStop(self.border_offset_+3, QtGui.QColor("Red"))
    self.AddStop(self.width()-(self.border_offset_+3), QtGui.QColor("Green"))

    #Width Problem
    self.width_ = self.width()

  def RemoveStop(self, stop):
    self.stops.remove(stop)
    stop.hide()
    del(stop)

    self.UpdateGradient()
    self.gradientUpdated.emit()
      
  def RemoveStopGui(self, stop):
    if(len(self.stops)>2):
      self.RemoveStop(stop)
    else:
      QtWidgets.QMessageBox.question(self, "Information", "Please keep in mind, at least two stops are needed for a gradient!")

  def AddStop(self, pos, color=None):
    stop = MyGradientStop(pos, self.border_offset_, self, color)
    stop.gradientChanged.connect(self.UpdateGradient)
    stop.colorChanged.connect(self.UpdateGradient)
    stop.colorChanged.connect(self.parent().Update)
    stop.gradientUpdated.connect(self.parent().Update)
    stop.removeStop.connect(self.RemoveStopGui)
    stop.MoveToNewPos()
    
    self.stops.append(stop)
    self.UpdateGradient() 
    self.gradientUpdated.emit()
  
  def Add(self):
    self.AddStop(QtWidgets.QCursor.pos().x() - self.mapToGlobal(QtCore.QPoint(0, 0)).x())
  
  def UpdateGradient(self):
    self.gradient_preview_.SetGradient(self.GetGradient())      

  def GetGfxGradient(self):  
    gradient = gfx.Gradient()
    for s in self.stops:
      c=s.GetColor();
      gradient.SetColorAt(s.GetRel(), gfx.RGB(c.redF(), c.greenF(), c.blueF()));
    return gradient;

  def GetGradient(self):
    gradient = QtGui.QLinearGradient(self.border_offset_, 
                                     0,
                                     self.size().width() - 2 * self.border_offset_, 
                                     0)
    for s in self.stops:
      gradient.setColorAt(s.GetRel(), s.GetColor())
    return gradient

  def LoadGradient(self, gradient):
    stopcpy = self.stops[:]
    for stop in stopcpy:
      self.RemoveStop(stop)
    
    gradStops = gradient.stops()    
    for stop in gradStops:
      self.AddStop((self.width()-(2*self.border_offset_))*stop[0], stop[1])
  
  #Overwritten methods
  def mouseDoubleClickEvent(self, event):
    self.AddStop(event.x())

  def resizeEvent(self, event):
    factor = 1/(float(self.width_)/event.size().width())
    if(factor >=0):
      for s in self.stops:
        s.SetPos(s.GetPos() * factor)
      self.UpdateGradient()
      self.width_ = event.size().width()
    
#Gradient Stop  
class MyGradientStop(ColorSelectWidget):

  removeStop = QtCore.pyqtSignal()
  gradientChanged = QtCore.pyqtSignal()
  gradientUpdated = QtCore.pyqtSignal()

  def __init__(self, pos, offset, parent, color=None):
    #Defaults
    self.length_ = 20
    self.halflength_ = self.length_ / 2
    
    self.pos_ = pos
    self.offset_ = offset
    
    ColorSelectWidget.__init__(self, self.length_,self.length_, color, parent)
    
    if(color is None):
      self.color_ = QtGui.QColor("White")
      self.ChangeColor()
    else:
      self.color_ = color
    
    #ContextMenu    
    self.remove_ = QtWidgets.QAction("Remove", self)
    self.remove_.triggered.connect(self.Remove)
    self.addAction(self.remove_)
    
  def Remove(self):
    self.removeStop.emit()

  def GetPos(self):
    return self.pos_
  
  def SetPos(self, pos):
    self.pos_ = pos
    self.MoveToNewPos()
    
  def GetRel(self):
    rel = self.pos_ / float(self.parent().width() - 2 * self.offset_)
    if rel > 1:
      return 1
    elif rel < 0:
      return 0
    return rel
  
  def MoveToNewPos(self):
    self.move(int(self.pos_ - self.halflength_), 0)
    self.update()
    self.gradientChanged.emit()
 
  def mouseMoveEvent(self, event):
    self.pos_ += event.pos().x() - self.halflength_
    if self.pos_ <= self.offset_:
      self.pos_ = self.offset_ + 1
    elif self.pos_ >= self.parent().width() - 2*self.offset_:
      self.pos_ = self.parent().width() - (2*self.offset_ + 1)
    self.MoveToNewPos()
    
  def mouseDoubleClickEvent(self, event):    
    self.ChangeColor()
    
  def mouseReleaseEvent(self, entity):
    if entity.button() == QtCore.Qt.LeftButton:
      self.gradientUpdated.emit()
