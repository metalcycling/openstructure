from ost import gui
from ost import gfx
import os
import ost
from PyQt5 import QtCore, QtWidgets, QtGui
from .immutable_preset_info_handler import ImmutablePresetInfoHandler
from .preset_info_handler import PresetInfoHandler
class PresetListModel(QtCore.QAbstractListModel):
  
  IMMUTABLE_PRESET_PATH = os.path.join(ost.GetSharedDataPath(),"scene", 
                                          "presets.xml")
  MUTABLE_PRESET_PATH = "user_presets.xml"

  dataChanged = QtCore.pyqtSignal(int, int, name="dataChanged")
  
  def __init__(self, parent=None, *args): 
    QtCore.QAbstractListModel.__init__(self, parent, *args)
    
    self.data_ = list()
    
    #Info Handler

    self.immutable_infoh_=ImmutablePresetInfoHandler(PresetListModel.IMMUTABLE_PRESET_PATH)
    data_loc=str(QtCore.QStandardPaths.writableLocation(QtCore.QStandardPaths.DataLocation))

    mutable_path=os.path.join(str(data_loc), 'config',
                              PresetListModel.MUTABLE_PRESET_PATH)
    qdir=QtCore.QDir(data_loc)
    if not qdir.exists():
      qdir.mkpath("config")
    self.infoh_=PresetInfoHandler(mutable_path)
    self.LoadPresetsFromInfo()
    
  def AddItem(self, preset, row, editable, save):
    if self.NameIsValid(preset.GetName()):
      icon = self.GetIcon(preset)
      self.insertRow(row, QtCore.QModelIndex())
      self.data_[row] = [preset, icon, editable]
      model_index = self.createIndex(row,0)
      index = self.index(row)
      end_index = self.createIndex(self.rowCount(),0)
      if save:
        self.AddPresetToInfo(preset)
      self.dataChanged.emit(model_index, end_index)
      return True
    return False
  
  def IsEditable(self, row):
    return self.data_[row][2]
  
  def RemoveItem(self, row):
    if self.IsEditable(row):
      name = self.data_[row][0].GetName()
      self.removeRow(row, QtCore.QModelIndex())
      model_index = self.createIndex(row,0)
      self.infoh_.RemovePreset(name)
      self.dataChanged.emit(model_index, model_index)
      return True
    return False

  def AddPresetToInfo(self,preset):
    self.infoh_.StorePreset(preset)
    
  def LoadPresetsFromInfo(self):
    if self.immutable_infoh_:
      presets = self.immutable_infoh_.GetPresets()
      for k, v in presets.items():
        self.AddItem(v, self.GetLastRow(), False, False)
    
    presets = self.infoh_.GetPresets()
    for k, v in presets.items():
      self.AddItem(v, self.GetLastRow(), True, False)

  def GetPreset(self, model_index):
    if model_index.isValid():
      return self.data_[model_index.row()][0]
  
  def GetLastRow(self):
    return self.rowCount()
  
  #Helper
  def GetIcon(self, preset):
    pixmap = QtGui.QPixmap(64, 64)
    pixmap.fill(QtCore.Qt.transparent)
    return QtGui.QIcon(pixmap)
  
  def NameIsValid(self, string):
    if len(string)==0:
      return False
    for values in self.data_:
      if string == values[0].GetName():
        return False
    return True
  
  #Overwritten Methods  
  def rowCount(self, parent=QtCore.QModelIndex()): 
    return len(self.data_) 
 
  def data(self, index, role):
    if index.isValid() and index.row()< self.rowCount():
      data = self.data_[index.row()]
      if role == QtCore.Qt.DisplayRole:
        return QtCore.QVariant(data[0].GetName())
      elif role == QtCore.Qt.DecorationRole:
        return QtCore.QVariant(data[1])
    return QtCore.QVariant()

  def setData(self, index, value, role):
    if index.isValid():
      row = index.row()
      if not self.data_[row]:
        self.data_[row] = list()
      if role == QtCore.Qt.EditRole and self.NameIsValid(value.toString()):
        old_name = str(self.data_[row][0].GetName())
        new_name = value.toString()
        self.data_[row][0].SetName(str(new_name))
        self.infoh_.RenamePreset(old_name,str(new_name))
        self.dataChanged.emit(index, index)
        return True
      elif role == QtCore.Qt.DisplayRole:
        self.data_[row][0].SetName(value.toString())
      elif role == QtCore.Qt.DecorationRole:
        self.dat_[row][1] = value.toPyObject()
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

  
  
  
