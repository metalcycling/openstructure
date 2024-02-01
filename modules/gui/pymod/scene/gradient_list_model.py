from ost import gui
from ost import gfx
import os
import ost
from PyQt5 import QtCore, QtWidgets, QtGui
from .immutable_gradient_info_handler import ImmutableGradientInfoHandler
from .gradient_info_handler import GradientInfoHandler

class GradientListModel(QtCore.QAbstractListModel):
  
  IMMUTABLE_GRADIENTS_PATH = os.path.join(ost.GetSharedDataPath(),"scene", 
                                          "gradients.xml")
  MUTABLE_GRADIENTS_PATH = "user_gradients.xml"

  dataChanged = QtCore.pyqtSignal(int, int, name="dataChanged")
  
  def __init__(self, parent=None, *args): 
    QtCore.QAbstractListModel.__init__(self, parent, *args)
    
    self.data_ = list()
    
    #Info Handler
    self.immutable_infoh_ = ImmutableGradientInfoHandler(GradientListModel.IMMUTABLE_GRADIENTS_PATH)
    self.infoh_ = GradientInfoHandler(GradientListModel.MUTABLE_GRADIENTS_PATH)
    self.LoadGradientFromInfo()
    
  def AddItem(self, name, gradient, row, editable, save):
    if self.NameIsValid(name):
      icon = self.GetIcon(gradient)
      self.insertRow(row, QtCore.QModelIndex())
      self.data_[row] = [name, gradient, icon, editable]
      model_index = self.createIndex(row,0)
      index = self.index(row)
      end_index = self.createIndex(self.rowCount(),0)
      if save:
        self.AddGradientToInfo(gradient,name)
      self.dataChanged.emit(model_index, end_index)
      return True
    return False
  
  def IsEditable(self, row):
    return self.data_[row][3]
  
  def RemoveItem(self, row):
    if self.IsEditable(row):
      name = self.data_[row][0]
      self.removeRow(row, QtCore.QModelIndex())
      model_index = self.createIndex(row,0)
      self.infoh_.RemoveGradient(name)
      self.dataChanged.emit(model_index, model_index)
      return True
    return False

  def AddGradientToInfo(self,gradient,name):
    self.infoh_.StoreQGradient(gradient,str(name))
    
  def LoadGradientFromInfo(self):
    if self.immutable_infoh_:
      qgradients = self.immutable_infoh_.GetQGradients()
      for k, v in qgradients.items():
        self.AddItem(k, v, self.GetLastRow(), False, False)
    
    qgradients = self.infoh_.GetQGradients()
    for k, v in qgradients.items():
      self.AddItem(k, v, self.GetLastRow(), True, False)

  def GetGradient(self, model_index):
    if model_index.isValid():
      return self.data_[model_index.row()][1]
  
  def GetLastRow(self):
    return self.rowCount()
  
  #Helper
  def GetIcon(self, gradient):
    pixmap = QtGui.QPixmap(64, 64)
    pixmap.fill(QtCore.Qt.transparent)
    painter = QtGui.QPainter()
    if painter.begin(pixmap):
      gradient.setStart(QtCore.QPointF(0, 0))
      gradient.setFinalStop(QtCore.QPointF(64, 0))
      brush = QtGui.QBrush(gradient)
      painter.setBrush(brush)
      painter.drawRect(0, 0, 64, 64)
      painter.end()
    
    return QtGui.QIcon(pixmap)
  
  def NameIsValid(self, string):
    if len(string)==0:
      return False
    for values in self.data_:
      if string == values[0]:
        return False
    return True
  
  #Overwritten Methods  
  def rowCount(self, parent=QtCore.QModelIndex()): 
    return len(self.data_) 
 
  def data(self, index, role):
    if index.isValid() and index.row()< self.rowCount():
      data = self.data_[index.row()]
      if role == QtCore.Qt.DisplayRole:
        return QtCore.QVariant(data[0])
      elif role == QtCore.Qt.DecorationRole:
        return QtCore.QVariant(data[2])
    return QtCore.QVariant()

  def setData(self, index, value, role):
    if index.isValid():
      row = index.row()
      if not self.data_[row]:
        self.data_[row] = list()
      if role == QtCore.Qt.EditRole and self.NameIsValid(value.toString()):
        old_name = str(self.data_[row][0])
        new_name = value.toString()
        self.data_[row][0] = new_name
        self.infoh_.RenameGradient(old_name,str(new_name))
        self.dataChanged.emit(index, index)
        return True
      elif role == QtCore.Qt.DisplayRole:
        self.data_[row][0] = value.toString()
      elif role == QtCore.Qt.DecorationRole:
        self.data_[row][2] = value.toPyObject()
    return False

  def flags(self, index):
    if index.isValid():
      flags = QtCore.QAbstractItemModel.flags(self,index)
      if self.IsEditable(index.row()):
        return flags | QtCore.Qt.ItemIsEditable
      else:
        return flags
    return QtCore.Qt.ItemIsEnabled
  
  def insertRow(self, position, index):
    self.beginInsertRows(index, position, position)
    self.data_.insert(position,list())
    self.endInsertRows()
    return True
  
  def removeRow(self, position, index):
    self.beginRemoveRows(index, position, position)
    del self.data_[position]
    self.endRemoveRows()
    return True

  
  
  
