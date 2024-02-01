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

from ost.gfx import ColorOp

#Rendering Preset
class Preset:
  NAME_ATTRIBUTE_NAME = "Name"
  OP_GROUP_NAME = "Op"
  CLASS_NAME_ATTRIBUTE_NAME = "ClassName"
  INDEX_ATTRIBUTE_NAME = "Index"
  
  MODULE_NAMES = ["ost.gfx","ost.gui.scene.visibility_op","ost.gui.scene.render_op"]
  
  def __init__(self, name, parent=None):    
    self.name_ = name
    self.ops_ = list()
    
  def SetName(self, name):
    self.name_ = name
    
  def GetName(self):
    return self.name_
    
  def InsertOp(self, index, op):
    self.ops_.insert(index, op)
    
  def RemoveOp(self, op):
    self.ops_.remove(op)
  
  def RemoveOpAt(self, index):
    del(self.ops_[index])
    
  def GetOp(self, index):
    return self.ops_[index]
  
  def GetOpCount(self):
    return len(self.ops_)
  
  def SetOp(self, index, op):
    self.ops_[index] = op
  
  def AddOp(self, op):
    self.ops_.append(op)
  
  def GetOps(self):
    return self.ops_
  
  def ApplyOn(self, entity):
    if (entity is not None) and isinstance(entity, gfx.Entity):
      for op in self.ops_:
        if isinstance(op,ColorOp):
          entity.Apply(op)
        else:
          op.ApplyOn(entity)
  
  def ToInfo(self,group):
    group.SetAttribute(Preset.NAME_ATTRIBUTE_NAME, self.name_)
    for i in range(0,len(self.ops_)):
      op_group = group.CreateGroup(Preset.OP_GROUP_NAME)
      op_group.SetAttribute(Preset.INDEX_ATTRIBUTE_NAME, str(i))
      op_group.SetAttribute(Preset.CLASS_NAME_ATTRIBUTE_NAME, "%s"%(self.ops_[i].__class__.__name__))
      self.ops_[i].ToInfo(op_group)
    
  @staticmethod
  def FromInfo(group):
    preset = None
    if group.HasAttribute(Preset.NAME_ATTRIBUTE_NAME):
      name = group.GetAttribute(Preset.NAME_ATTRIBUTE_NAME)
      preset = Preset(name)
      group_list = group.GetGroups(Preset.OP_GROUP_NAME)
      
      class_order_dict = dict()
      for op_group in group_list:
        if(op_group.HasAttribute(Preset.CLASS_NAME_ATTRIBUTE_NAME) and op_group.HasAttribute(Preset.INDEX_ATTRIBUTE_NAME)):
          class_name = op_group.GetAttribute(Preset.CLASS_NAME_ATTRIBUTE_NAME)
          index = int(op_group.GetAttribute(Preset.INDEX_ATTRIBUTE_NAME))
          op_class = None
          for module in Preset.MODULE_NAMES:
            try:
              op_class = Preset.__get_op_class("%s.%s"%(module,class_name))
              break
            except AttributeError:
              pass
          if op_class is not None:
            op = op_class.FromInfo(op_group)
            class_order_dict[index]=op
    for i in range(0, len(class_order_dict)):
      if(i in class_order_dict):
        preset.AddOp(class_order_dict[i])
    return preset
    
  @staticmethod
  def __get_op_class( cls ):
    parts = cls.split('.')
    module = ".".join(parts[:-1])
    m = __import__( module )
    for comp in parts[1:]:
      m = getattr(m, comp)
    return m
