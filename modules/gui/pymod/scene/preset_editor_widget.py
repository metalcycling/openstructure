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
from datetime import datetime
from datetime import datetime
from PyQt5 import QtCore, QtWidgets, QtGui
from .color_select_widget import ColorSelectWidget
from .gradient_preset_widget import GradientPresetWidget
from .gradient_editor_widget import GradientPreview
from .gradient_editor_widget import GradientEdit
from .preset_editor_list_model import PresetEditorListModel
from .immutable_gradient_info_handler import ImmutableGradientInfoHandler
from .query_editor import QueryEditorWidget
from ost.mol import Prop
from ost.gfx import ByElementColorOp
from ost.gfx import ByChainColorOp
from ost.gfx import GradientLevelColorOp
from ost.gfx import UniformColorOp
from .preset import Preset
from .render_op import RenderOp
from .visibility_op import VisibilityOp

#Preset Editor
class PresetEditor(QtWidgets.QDialog):
  def __init__(self, parent=None):
    QtWidgets.QDialog.__init__(self, parent)
        
    self.setWindowTitle("Preset Editor")
    
    #Create Ui Elements
    self.list_view_ = QtWidgets.QListView()
    
    self.combo_box_ = QtWidgets.QComboBox()
    
    self.ufcow_=UniformColorOpWidget(self)
    self.glcow_=GradientLevelColorOpWidget(self)
    self.beow_=ByElementColorOpWidget(self)
    self.bcow_=ByChainColorOpWidget(self)
    self.row_=RenderOpWidget(self)
    self.vow_=VisibilityOpWidget(self)
    self.combo_box_.addItem("Uniform Color Operation", QtCore.QVariant(self.ufcow_))
    self.combo_box_.addItem("Gradient Operation", QtCore.QVariant(self.glcow_))
    self.combo_box_.addItem("By Element Operation", QtCore.QVariant(self.beow_))
    self.combo_box_.addItem("By Chain Operation", QtCore.QVariant(self.bcow_))
    self.combo_box_.addItem("RenderMode Operation", QtCore.QVariant(self.row_))
    self.combo_box_.addItem("Visibility Operation", QtCore.QVariant(self.vow_))
    
    self.add_button_ = QtWidgets.QPushButton("Add")
    
    #Create Model
    self.list_view_.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
    
    self.list_view_.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
    self.list_view_.customContextMenuRequested.connect(self.contextMenuEvent)
    
    self.hbox_ = QtWidgets.QHBoxLayout()
    self.ok_button_ = QtWidgets.QPushButton("OK")
    self.cancel_button_ = QtWidgets.QPushButton("Cancel")
    self.hbox_.addWidget(self.ok_button_)
    self.hbox_.addStretch()
    self.hbox_.addWidget(self.cancel_button_)
    
    grid = QtWidgets.QGridLayout()
    grid.setContentsMargins(0,5,0,0)
    grid.addWidget(self.combo_box_,0,0,1,1)
    grid.addWidget(self.add_button_,0,1,1,1)
    grid.addWidget(self.list_view_,1,0,3,3)
    grid.addLayout(self.hbox_,5,0,1,3)
    grid.setRowStretch(1, 1)
    self.setLayout(grid)
    
    self.add_button_.clicked.connect(self.Add)
    self.ok_button_.clicked.connect(self.Ok)
    self.cancel_button_.clicked.connect(self.Cancel)
  
    self.CreateContextMenu()
      
  def SetPreset(self, preset):
    self.list_model_ = PresetEditorListModel(preset, self)
    self.list_view_.setModel(self.list_model_)
  
  def CreateContextMenu(self):
    self.context_menu_ = QtWidgets.QMenu("Context menu", self)
    self.edit_ = QtWidgets.QAction("Edit", self.list_view_)
    self.remove_ = QtWidgets.QAction("Remove", self.list_view_)
    self.moveup_ = QtWidgets.QAction("Move Up", self.list_view_)
    self.movedown_ = QtWidgets.QAction("Move Down", self.list_view_)
    self.context_menu_.addAction(self.edit_)
    self.context_menu_.addAction(self.remove_)
    self.context_menu_.addAction(self.moveup_)
    self.context_menu_.addAction(self.movedown_)
    #Connect Signals with Slots
    self.edit_.triggered.connect(self.Edit)
    self.remove_.triggered.connect(self.Remove)
    self.moveup_.triggered.connect(self.MoveUp)
    self.movedown_.triggered.connect(self.MoveDown)

  def contextMenuEvent(self, pos):
    #ContextMenu
    index = self.list_view_.indexAt(pos)
    if index.isValid(): 
        self.context_menu_.popup(QtWidgets.QCursor.pos())

  def Add(self):
    dialog = self.combo_box_.itemData(self.combo_box_.currentIndex()).toPyObject()
    if(dialog.exec_()):
      row = self.list_model_.rowCount()
      op = dialog.GetOp()
      self.list_model_.AddItem(op, row)
  
  def Edit(self):
    current_index = self.list_view_.currentIndex()
    op = self.list_model_.GetOp(current_index)
    if isinstance(op, gfx.UniformColorOp):
       self.ufcow_.SetOp(op)
       if self.ufcow_.exec_():
         self.list_model_.SetItem(current_index, self.ufcow_.GetOp())
    elif isinstance(op, gfx.GradientLevelColorOp):
       self.glcow_.SetOp(op)
       if self.glcow_.exec_():
         self.list_model_.SetItem(current_index, self.glcow_.GetOp())
    elif isinstance(op, gfx.ByElementColorOp):
       self.beow_.SetOp(op)
       if self.beow_.exec_():
         self.list_model_.SetItem(current_index, self.beow_.GetOp())
    elif isinstance(op, gfx.ByChainColorOp):
       self.bcow_.SetOp(op)
       if self.bcow_.exec_():
         self.list_model_.SetItem(current_index, self.bcow_.GetOp())
    elif isinstance(op, RenderOp):
       self.row_.SetOp(op)
       if self.row_.exec_():
         self.list_model_.SetItem(current_index, self.row_.GetOp())
    elif isinstance(op, VisibilityOp):
       self.vow_.SetOp(op)
       if self.vow_.exec_():
         self.list_model_.SetItem(current_index, self.vow_.GetOp())
    
  def Remove(self):
    current_index = self.list_view_.currentIndex()
    self.list_model_.RemoveItem(current_index.row())
    
  def MoveDown(self):
    current_index = self.list_view_.currentIndex()
    if self.list_model_.GetLastRow != current_index.row():
      op = self.list_model_.GetOp(current_index)
      self.list_model_.RemoveItem(current_index.row())
      self.list_model_.AddItem(op, current_index.row()+1)
 
  def MoveUp(self):
    current_index = self.list_view_.currentIndex()
    if self.list_model_.GetLastRow != 0:
      op = self.list_model_.GetOp(current_index)
      self.list_model_.RemoveItem(current_index.row())
      self.list_model_.AddItem(op, current_index.row()-1)
           
  def Ok(self):
    self.accept()
    
  def Cancel(self):
    self.reject()

class UniformColorOpWidget(QtWidgets.QDialog):
  def __init__(self, parent=None):
    QtWidgets.QDialog.__init__(self, parent)
    self.query_editor_ = QueryEditorWidget(self)

    detail_label = QtWidgets.QLabel("Parts")
    self.detail_selection_cb_ = QtWidgets.QComboBox()
    self.detail_selection_cb_.addItem("Main and Detail",QtCore.QVariant(3))
    self.detail_selection_cb_.addItem("Main",QtCore.QVariant(2))
    self.detail_selection_cb_.addItem("Detail",QtCore.QVariant(1))
    
    color_label = QtWidgets.QLabel("Color")
    self.color_select_widget_ = ColorSelectWidget(30,30,QtGui.QColor("White"))    
    
    self.hbox_ = QtWidgets.QHBoxLayout()
    self.ok_button_ = QtWidgets.QPushButton("OK")
    self.cancel_button_ = QtWidgets.QPushButton("Cancel")
    self.hbox_.addWidget(self.ok_button_)
    self.hbox_.addStretch()
    self.hbox_.addWidget(self.cancel_button_)

        
    grid = QtWidgets.QGridLayout()
    grid.setContentsMargins(0,5,0,0)
    grid.addWidget(self.query_editor_, 0, 0, 1, 2)
    grid.addWidget(detail_label, 1, 0, 1, 1)
    grid.addWidget(self.detail_selection_cb_, 1, 1, 1, 1)
    grid.addWidget(color_label, 2, 0, 1, 1)
    grid.addWidget(self.color_select_widget_, 2, 1, 1, 1)
    grid.addLayout(self.hbox_,3,0,1,2)
    grid.setRowStretch(2, 1)
    self.setLayout(grid)
    
    self.ok_button_.clicked.connect(self.Ok)
    self.cancel_button_.clicked.connect(self.Cancel)
    
  def GetOp(self):
    qv=mol.QueryViewWrapper(self.query_editor_.GetQuery(),self.query_editor_.GetQueryFlags())
    ufco = UniformColorOp(qv,gfx.RGBA(1,1,1,1))
    
    detail = self.detail_selection_cb_.itemData(self.detail_selection_cb_.currentIndex()).toPyObject()
    ufco.SetMask(detail)
    qcolor = self.color_select_widget_.GetColor()
    color=gfx.RGBAb(qcolor.red(),qcolor.green(),qcolor.blue(),qcolor.alpha())
    ufco.SetColor(color)
    return ufco

  def SetOp(self, ufco):
    self.query_editor_.SetQuery(ufco.GetSelection())
    self.query_editor_.SetQueryFlags(ufco.GetSelectionFlags())
    found=False
    for i in range(0,self.detail_selection_cb_.count()):
      mask = self.detail_selection_cb_.itemData(i).toPyObject()
      if mask == ufco.GetMask():
        self.detail_selection_cb_.setCurrentIndex(i)
        found = True
        break
    if not found:
      self.detail_selection_cb_.setCurrentIndex(0)
    
    color = ufco.GetColor()
    qcolor = QtGui.QColor(color.Red()*255, color.Green()*255, color.Blue()*255, color.Alpha()*255)
    self.color_select_widget_.SetColor(qcolor)

  def Ok(self):
    self.accept()
    
  def Cancel(self):
    self.reject()

class GradientLevelColorOpWidget(QtWidgets.QDialog):
  def __init__(self, parent=None):
    QtWidgets.QDialog.__init__(self, parent)
    
    self.query_editor_ = QueryEditorWidget(self)
    
    detail_label = QtWidgets.QLabel("Parts")
    self.detail_selection_cb_ = QtWidgets.QComboBox()
    self.detail_selection_cb_.addItem("Main and Detail",QtCore.QVariant(3))
    self.detail_selection_cb_.addItem("Main",QtCore.QVariant(2))
    self.detail_selection_cb_.addItem("Detail",QtCore.QVariant(1))
    
    property_label = QtWidgets.QLabel("Property")
    self.property_edit_ = QtWidgets.QLineEdit()

    self.prop_combo_box_ = QtWidgets.QComboBox()
    self.prop_combo_box_.addItem("atom B-factor",QtCore.QVariant("abfac"))
    self.prop_combo_box_.addItem("average residue B-factor",QtCore.QVariant("rbfac"))
    self.prop_combo_box_.addItem("X-Coordinate",QtCore.QVariant("x"))
    self.prop_combo_box_.addItem("Y-Coordinate",QtCore.QVariant("y"))
    self.prop_combo_box_.addItem("Z-Coordinate",QtCore.QVariant("z"))
    self.prop_combo_box_.addItem("Residue Number",QtCore.QVariant("rnum"))
    self.prop_combo_box_.addItem("Atom Charge",QtCore.QVariant("acharge"))
    self.prop_combo_box_.addItem("Custom",QtCore.QVariant("custom"))
    

    level_label = QtWidgets.QLabel("Level")
    self.combo_box_ = QtWidgets.QComboBox(self);
    self.combo_box_.addItem("Atom",QtCore.QVariant(Prop.Level.ATOM))
    self.combo_box_.addItem("Residue",QtCore.QVariant(Prop.Level.RESIDUE))
    self.combo_box_.addItem("Chain",QtCore.QVariant(Prop.Level.CHAIN))
    self.combo_box_.addItem("Unspecified",QtCore.QVariant(Prop.Level.UNSPECIFIED))
        
    gradient_label = QtWidgets.QLabel("Gradient")
    self.gradient_preview_ = GradientPreview()    
    self.gradient_edit_ = GradientEdit(self.gradient_preview_,self)
    
    self.minmax_label_ = QtWidgets.QLabel("Min Max")
    self.auto_calc_ = QtWidgets.QCheckBox("Auto calculate")
    self.auto_calc_.setChecked(True)
    
    self.minv_label_ = QtWidgets.QLabel("Min Value")
    self.maxv_label_ = QtWidgets.QLabel("Max Value")
    
    self.minv_ = QtWidgets.QDoubleSpinBox(self)
    self.minv_.setDecimals(2)
    self.minv_.setMinimum(-9999.99)
    self.minv_.setValue(0)
    self.maxv_ = QtWidgets.QDoubleSpinBox(self)
    self.maxv_.setDecimals(2)
    self.maxv_.setValue(1)
    self.maxv_.setMinimum(-9999.99)
    
    self.hbox_ = QtWidgets.QHBoxLayout()
    self.ok_button_ = QtWidgets.QPushButton("OK")
    self.cancel_button_ = QtWidgets.QPushButton("Cancel")
    self.hbox_.addWidget(self.ok_button_)
    self.hbox_.addStretch()
    self.hbox_.addWidget(self.cancel_button_)

        
    grid = QtWidgets.QGridLayout()
    grid.setContentsMargins(0,5,0,0)
    grid.addWidget(self.query_editor_, 0, 0, 1, 2)
    grid.addWidget(detail_label, 1, 0, 1, 1)
    grid.addWidget(self.detail_selection_cb_, 1, 1, 1, 1)
    grid.addWidget(property_label, 2, 0, 1, 1)
    grid.addWidget(self.prop_combo_box_, 2, 1, 1, 1)
    grid.addWidget(self.property_edit_, 3, 1, 1, 1)
    grid.addWidget(level_label, 4, 0, 1, 1)
    grid.addWidget(self.combo_box_, 4, 1, 1, 1)
    grid.addWidget(gradient_label, 5, 0, 1, 1)
    grid.addWidget(self.gradient_preview_, 5, 1, 1, 1)
    grid.addWidget(self.gradient_edit_, 6, 1, 1, 1)
    grid.addWidget(self.gradient_edit_, 6, 1, 1, 1)
    grid.addWidget(self.minmax_label_,7,0,1,1)
    grid.addWidget(self.auto_calc_,7,1,1,1)
    grid.addWidget(self.minv_label_,8,0,1,1)
    grid.addWidget(self.minv_,8,1,1,1)
    grid.addWidget(self.maxv_label_,9,0,1,1)
    grid.addWidget(self.maxv_,9,1,1,1)
    grid.addLayout(self.hbox_,10,0,1,2)
    grid.setRowStretch(1, 1)
    self.setLayout(grid)
    
    self.prop_combo_box_.currentIndexChanged.connect(self.UpdateGui)
    self.ok_button_.clicked.connect(self.Ok)
    self.cancel_button_.clicked.connect(self.Cancel)
    self.auto_calc_.stateChanged.connect(self.UpdateGui)
    
    self.UpdateGui()
    
  def GetOp(self):
    gradient = self.gradient_edit_.GetGfxGradient()
    
    detail = self.detail_selection_cb_.itemData(self.detail_selection_cb_.currentIndex()).toPyObject()
    qv=mol.QueryViewWrapper(self.query_editor_.GetQuery(),self.query_editor_.GetQueryFlags())
    prop = str()
    level = Prop.Level()
    if(self.property_edit_.isEnabled()):
      prop = str(self.property_edit_.text())
      level = Prop.Level(self.combo_box_.itemData(self.combo_box_.currentIndex()).toPyObject())
    else:
      prop = str(self.prop_combo_box_.itemData(self.prop_combo_box_.currentIndex()).toPyObject())
    
    if(self.auto_calc_.isChecked()):
      return GradientLevelColorOp(qv, detail, prop, gradient, level)
    else:
      minv = self.minv_.value()
      maxv = self.maxv_.value()
      return GradientLevelColorOp(qv, detail, prop, gradient, minv, maxv, level)
      
    
  def SetOp(self, glco):
    self.query_editor_.SetQuery(glco.GetSelection())
    self.query_editor_.SetQueryFlags(glco.GetSelectionFlags())
    found=False
    for i in range(0,self.detail_selection_cb_.count()):
      mask = self.detail_selection_cb_.itemData(i).toPyObject()
      if mask == glco.GetMask():
        self.detail_selection_cb_.setCurrentIndex(i)
        found = True
        break
    if not found:
      self.detail_selection_cb_.setCurrentIndex(0)
    
    found = False
    for i in range(0,self.prop_combo_box_.count()):
      prop = str(self.prop_combo_box_.itemData(i).toPyObject())
      if prop == glco.GetProperty():
        self.prop_combo_box_.setCurrentIndex(i)
        found = True
    if not found:
      self.prop_combo_box_.setCurrentIndex(self.prop_combo_box_.count()-1)
      self.property_edit_.setText(glco.GetProperty())
      
    self.combo_box_.setCurrentIndex(glco.GetLevel())
    self.gradient_edit_.LoadGradient(ImmutableGradientInfoHandler.ConvertToQGradient(glco.GetGradient()))
    self.auto_calc_.setChecked(glco.GetCalculateMinMax());
    if(not glco.GetCalculateMinMax()):
      self.minv_.setValue(glco.GetMinV())
      self.maxv_.setValue(glco.GetMaxV())
    self.UpdateGui()
    
  def UpdateGui(self):
    prop = str(self.prop_combo_box_.itemData(self.prop_combo_box_.currentIndex()))
    if(prop == "custom"):
      self.combo_box_.setEnabled(True)
      self.property_edit_.setEnabled(True)
    else:
      self.combo_box_.setEnabled(False)
      self.property_edit_.setEnabled(False)
    
    if(self.auto_calc_.isChecked()):
      self.minv_.setEnabled(False)
      self.maxv_.setEnabled(False)
    else:
      self.minv_.setEnabled(True)
      self.maxv_.setEnabled(True)
      
  def Ok(self):
    self.accept()
    
  def Cancel(self):
    self.reject()
    
  def Update(self):
    pass #Do Nothing

class ByElementColorOpWidget(QtWidgets.QDialog):
  def __init__(self, parent=None):
    QtWidgets.QDialog.__init__(self, parent)
    self.query_editor_ = QueryEditorWidget(self)

    detail_label = QtWidgets.QLabel("Parts")
    self.detail_selection_cb_ = QtWidgets.QComboBox()
    self.detail_selection_cb_.addItem("Main and Detail",QtCore.QVariant(3))
    self.detail_selection_cb_.addItem("Main",QtCore.QVariant(2))
    self.detail_selection_cb_.addItem("Detail",QtCore.QVariant(1))
  
    self.hbox_ = QtWidgets.QHBoxLayout()
    self.ok_button_ = QtWidgets.QPushButton("OK")
    self.cancel_button_ = QtWidgets.QPushButton("Cancel")
    self.hbox_.addWidget(self.ok_button_)
    self.hbox_.addStretch()
    self.hbox_.addWidget(self.cancel_button_)
    
    grid = QtWidgets.QGridLayout()
    grid.setContentsMargins(0,5,0,0)
    grid.addWidget(self.query_editor_, 0, 0, 1, 2)
    grid.addWidget(detail_label, 1, 0, 1, 1)
    grid.addWidget(self.detail_selection_cb_, 1, 1, 1, 1)
    grid.addLayout(self.hbox_,2,0,1,2)
    grid.setRowStretch(1, 1)
    self.setLayout(grid)
    
    self.ok_button_.clicked.connect(self.Ok)
    self.cancel_button_.clicked.connect(self.Cancel)
    
  def GetOp(self):
    detail = self.detail_selection_cb_.itemData(self.detail_selection_cb_.currentIndex()).toPyObject()
    qv=mol.QueryViewWrapper(self.query_editor_.GetQuery(),self.query_editor_.GetQueryFlags())
    beco = ByElementColorOp(qv, detail)
    return beco

  def SetOp(self, beco):
    self.query_editor_.SetQuery(beco.GetSelection())
    self.query_editor_.SetQueryFlags(beco.GetSelectionFlags())
    found=False
    for i in range(0,self.detail_selection_cb_.count()):
      mask = self.detail_selection_cb_.itemData(i).toPyObject()
      if mask == beco.GetMask():
        self.detail_selection_cb_.setCurrentIndex(i)
        found = True
        break
    if not found:
      self.detail_selection_cb_.setCurrentIndex(0)
  
  def Ok(self):
    self.accept()
    
  def Cancel(self):
    self.reject()
    
    
class ByChainColorOpWidget(QtWidgets.QDialog):
  def __init__(self, parent=None):
    QtWidgets.QDialog.__init__(self, parent)
    self.query_editor_ = QueryEditorWidget(self)

    detail_label = QtWidgets.QLabel("Parts")
    self.detail_selection_cb_ = QtWidgets.QComboBox()
    self.detail_selection_cb_.addItem("Main and Detail",QtCore.QVariant(3))
    self.detail_selection_cb_.addItem("Main",QtCore.QVariant(2))
    self.detail_selection_cb_.addItem("Detail",QtCore.QVariant(1))
  
    self.hbox_ = QtWidgets.QHBoxLayout()
    self.ok_button_ = QtWidgets.QPushButton("OK")
    self.cancel_button_ = QtWidgets.QPushButton("Cancel")
    self.hbox_.addWidget(self.ok_button_)
    self.hbox_.addStretch()
    self.hbox_.addWidget(self.cancel_button_)
    
    grid = QtWidgets.QGridLayout()
    grid.setContentsMargins(0,5,0,0)
    grid.addWidget(self.query_editor_, 0, 0, 1, 2)
    grid.addWidget(detail_label, 1, 0, 1, 1)
    grid.addWidget(self.detail_selection_cb_, 1, 1, 1, 1)
    grid.addLayout(self.hbox_,2,0,1,2)
    grid.setRowStretch(1, 1)
    self.setLayout(grid)
    
    self.ok_button_.clicked.connect(self.Ok)
    self.cancel_button_.clicked.connect(self.Cancel)
    
  def GetOp(self):
    detail = self.detail_selection_cb_.itemData(self.detail_selection_cb_.currentIndex()).toPyObject()
    qv=mol.QueryViewWrapper(self.query_editor_.GetQuery(),self.query_editor_.GetQueryFlags())
    bcco = ByChainColorOp(qv, detail)
    return bcco

  def SetOp(self, bcco):
    self.query_editor_.SetQuery(bcco.GetSelection())
    self.query_editor_.SetQueryFlags(bcco.GetSelectionFlags())
    found=False
    for i in range(0,self.detail_selection_cb_.count()):
      mask = self.detail_selection_cb_.itemData(i).toPyObject()
      if mask == bcco.GetMask():
        self.detail_selection_cb_.setCurrentIndex(i)
        found = True
        break
    if not found:
      self.detail_selection_cb_.setCurrentIndex(0)
  
  def Ok(self):
    self.accept()
    
  def Cancel(self):
    self.reject()


class RenderOpWidget(QtWidgets.QDialog):
  def __init__(self, parent=None):
    QtWidgets.QDialog.__init__(self, parent)
    self.query_editor_ = QueryEditorWidget(self)
  
    self.keep_ = QtWidgets.QCheckBox("Keep")
    self.keep_.setChecked(False)
    
    render_label = QtWidgets.QLabel("Rendermode")
    self.render_modes_ = QtWidgets.QComboBox()
    self.render_modes_.addItem("Fast Bonds")
    self.render_modes_.addItem("Ball & Stick")
    self.render_modes_.addItem("Spheres")
    self.render_modes_.addItem("Fast Trace")
    self.render_modes_.addItem("Trace")
    self.render_modes_.addItem("Fast Spline")
    self.render_modes_.addItem("Smooth Tube")
    self.render_modes_.addItem("Helix & Strand Cartoon")
    
    self.render_modes_list_ = [gfx.RenderMode.SIMPLE,
                               gfx.RenderMode.CUSTOM,
                               gfx.RenderMode.CPK,
                               gfx.RenderMode.LINE_TRACE,
                               gfx.RenderMode.TRACE,
                               gfx.RenderMode.SLINE,
                               gfx.RenderMode.TUBE,
                               gfx.RenderMode.HSC]
    
    self.hbox_ = QtWidgets.QHBoxLayout()
    self.ok_button_ = QtWidgets.QPushButton("OK")
    self.cancel_button_ = QtWidgets.QPushButton("Cancel")
    self.hbox_.addWidget(self.ok_button_)
    self.hbox_.addStretch()
    self.hbox_.addWidget(self.cancel_button_)
    
    grid = QtWidgets.QGridLayout()
    grid.setContentsMargins(0,5,0,0)
    grid.addWidget(self.query_editor_, 0, 0, 1, 2)
    grid.addWidget(self.keep_, 1, 1, 1, 1)
    grid.addWidget(render_label,2,0,1,1)
    grid.addWidget(self.render_modes_,2,1,1,1)
    grid.addLayout(self.hbox_,3,0,1,2)
    grid.setRowStretch(1, 1)
    self.setLayout(grid)
    
    self.ok_button_.clicked.connect(self.Ok)
    self.cancel_button_.clicked.connect(self.Cancel)
    
  def GetOp(self):
    selection = self.query_editor_.GetQueryText()
    flags = self.query_editor_.GetQueryFlags()
    render_mode = self.render_modes_list_[self.render_modes_.currentIndex()]
    ro = RenderOp(render_mode, selection, flags, self.keep_.isChecked())
    return ro

  def SetOp(self, ro):
    self.query_editor_.SetQuery(ro.GetSelection())
    self.query_editor_.SetQueryFlags(ro.GetSelectionFlags())
    found=False
    for i in range(0,self.render_modes_.count()):
      render_mode = self.render_modes_list_[i]
      if render_mode == ro.GetRenderMode():
        self.render_modes_.setCurrentIndex(i)
        found = True
        break
    if not found:
      self.render_modes_.setCurrentIndex(0)
    self.keep_.setChecked(ro.IsKept())
    
  def Ok(self):
    self.accept()
    
  def Cancel(self):
    self.reject()
    
class VisibilityOpWidget(QtWidgets.QDialog):
  def __init__(self, parent=None):
    QtWidgets.QDialog.__init__(self, parent)
    self.query_editor_ = QueryEditorWidget(self)
  
    self.visible_ = QtWidgets.QCheckBox("Visible")
    self.visible_.setChecked(True)
    
    self.hbox_ = QtWidgets.QHBoxLayout()
    self.ok_button_ = QtWidgets.QPushButton("OK")
    self.cancel_button_ = QtWidgets.QPushButton("Cancel")
    self.hbox_.addWidget(self.ok_button_)
    self.hbox_.addStretch()
    self.hbox_.addWidget(self.cancel_button_)
    
    grid = QtWidgets.QGridLayout()
    grid.setContentsMargins(0,5,0,0)
    grid.addWidget(self.query_editor_, 0, 0, 1, 2)
    grid.addWidget(self.visible_, 1, 1, 1, 1)
    grid.addLayout(self.hbox_,2,0,1,2)
    grid.setRowStretch(1, 1)
    self.setLayout(grid)
    
    self.ok_button_.clicked.connect(self.Ok)
    self.cancel_button_.clicked.connect(self.Cancel)
    
  def GetOp(self):
    selection = self.query_editor_.GetQueryText()
    flags = self.query_editor_.GetQueryFlags()
    vo = VisibilityOp(selection, flags, self.visible_.isChecked())
    return vo

  def SetOp(self, vo):
    self.query_editor_.SetQuery(vo.GetSelection())
    self.query_editor_.SetQueryFlags(vo.GetSelectionFlags())
    self.visible_.setChecked(vo.IsVisible())
  
  def Ok(self):
    self.accept()
    
  def Cancel(self):
    self.reject()
