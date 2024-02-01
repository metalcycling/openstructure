from ost import gui
from ost import gfx
import os
import ost
from PyQt5 import QtCore
from .preset import Preset

class PresetEditorListModel(QtCore.QAbstractListModel):
  def __init__(self, preset, parent=None, *args): 
    QtCore.QAbstractListModel.__init__(self, parent, *args)
    
    self.preset_ = preset
    
  def AddItem(self, op, row):
    self.insertRow(row, QtCore.QModelIndex(),op)
    model_index = self.createIndex(row,0)
    end_index = self.createIndex(self.rowCount(),0)
    self.emit(QtCore.SIGNAL("dataChanged"),model_index, end_index)
    return True
  
  def RemoveItem(self, row):
    self.removeRow(row, QtCore.QModelIndex())
    model_index = self.createIndex(row,0)
    self.emit(QtCore.SIGNAL("dataChanged"),model_index, model_index)
    return True
    
  def GetOp(self, model_index):
    if model_index.isValid():
      return self.preset_.GetOp(model_index.row())
  
  def SetItem(self, model_index, op):
    self.preset_.SetOp(model_index.row(), op)
  
  def GetLastRow(self):
    return self.rowCount()
    
  #Overwritten Methods  
  def rowCount(self, parent=QtCore.QModelIndex()): 
    return self.preset_.GetOpCount() 
 
  def data(self, index, role):
    if index.isValid() and index.row()< self.rowCount():
      data = self.preset_.GetOp(index.row())
      if role == QtCore.Qt.DisplayRole:
        name=data.GetName()
        selection=str(data.GetSelection())
        if(len(selection)>0):
          return QtCore.QVariant(str("%s (%s)"%(name,selection.replace(os.linesep, " "))))
        else:
          return QtCore.QVariant(str("%s (all)"%name))
    return QtCore.QVariant()

  def setData(self, index, value, role):
    if index.isValid():
      row = index.row()
      data = self.preset_.GetOp(row)
      if role == QtCore.Qt.EditRole:
        return True
      elif role == QtCore.Qt.DisplayRole:
        data
      elif role == QtCore.Qt.DecorationRole:
        data
    return False

  def flags(self, index):
    if index.isValid():
      flags = QtCore.QAbstractItemModel.flags(self,index)
      return flags | QtCore.Qt.ItemIsEditable
    return QtCore.Qt.ItemIsEnabled
  
  def insertRow(self, position, index, op):
    self.beginInsertRows(index, position, position)
    self.preset_.InsertOp(position, op)
    self.endInsertRows()
    return True
  
  def removeRow(self, position, index):
    self.beginRemoveRows(index, position, position)
    self.preset_.RemoveOpAt(position)
    self.endRemoveRows()
    return True
  
