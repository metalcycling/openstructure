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

from ost import gui
from ost import info
import ost
import re
from PyQt5 import QtCore, QtGui, QtNetwork
from ost.gui import FileLoader

"""
class BaseRemoteLoader(gui.RemoteSiteLoader):
  def __init__(self):
    gui.RemoteSiteLoader.__init__(self)
    self.networkmanager_=QtNetwork.QNetworkAccessManager()
    self.downloads_=dict()
    QtCore.QObject.connect(self.networkmanager_, QtCore.SIGNAL("finished(QNetworkReply*)"), self.DownloadFinished)

  def LoadById(self, id, selection=""):
    self.ById(id, selection)
  
  def ById(self, ids, selection=""):
    for id in re.findall('\w+', ids):
      file_name=self.GetFileName(id)
      file = QtCore.QFile(file_name)
      if(file.size()==0):
          url = QtCore.QUrl(self.GetUrl(id))
          request = QtNetwork.QNetworkRequest(url)
          reply = self.networkmanager_.get(request)
          self.downloads_[reply]=[id,selection]
      else:
        gui.FileLoader.LoadObject(str(file_name),str(selection))
    return None
    
  def IsImg(self):
    return False
    
  def IsDefault(self):
    return False

  def GetRemoteSiteName(self):
    pass  
  
  def GetUrl(self,id):
    pass
      
  def HandleError(self,message):
    pass
    
  def DownloadFinished(self, reply):
    file_name=self.GetFileName(self.downloads_[reply][0])
    file = QtCore.QFile(file_name)
    if(reply.error()!=QtNetwork.QNetworkReply.NoError or reply.bytesAvailable()==0):
      self.HandleError(reply.errorString());
      file.remove()
    else:
      if(self.downloads_.has_key(reply)):
        if(file.open(QtCore.QIODevice.WriteOnly)):
           file.write(reply.readAll());
        file.close();
        selection=self.downloads_[reply][1]
        del(self.downloads_[reply])
        gui.FileLoader.LoadObject(str(file_name),str(selection))

class GenericLoader(BaseRemoteLoader):
  EXT_NAME_ATTRIBUTE_NAME = "ExtName"
  URL_ATTRIBUTE_NAME = "Url"
  UP_CASE_ATTRIBUTE_NAME = "UpCase"
  FILE_TYPE_ATTRIBUTE_NAME = "FileType"
  TMP_URL_ATTRIBUTE_NAME= "TmpUrl"
  IMG_ATTRIBUTE_NAME= "Img"
  DEFAULT_ATTRIBUTE_NAME= "Default"
  def __init__(self, name, url, up_case, file_type, tmp_url=None, img=False, default=False):
    BaseRemoteLoader.__init__(self)
    self.name_ = QtCore.QString(name)
    self.url_ = QtCore.QString(url)
    self.up_case_ = up_case
    self.img_ = img
    self.default_ = default
    if file_type is not None:
      self.file_type_ = QtCore.QString(file_type)
    else:
      self.file_type_ = None
    
    if tmp_url is not None:
      self.tmp_url_ = QtCore.QString(tmp_url)
    else:
      self.tmp_url_ = tmp_url
    
  def IsImg(self):
    return self.img_    

  def IsDefault(self):
    return self.default_
  
  def GetRemoteSiteName(self):
    return self.name_
  
  def GetUrl(self, id):
    formatted_id = id
    if self.up_case_:
      formatted_id = QtCore.QString(id).toUpper()
    url = QtCore.QString(self.url_)
    url.replace("${ID}",formatted_id)
    return str(url)
    
  def GetFileName(self,id):
    formatted_id = id
    if self.up_case_:
      formatted_id = QtCore.QString(id).toUpper()
    
    url = None
    if self.tmp_url_ is not None:
      url = self.tmp_url_ + QtCore.QDir.separator()
    else:
      url = QtCore.QDir.tempPath() + QtCore.QDir.separator()
    
    file_type = None 
    if self.file_type_ is None:
      remote_url = QtCore.QString(self.GetUrl(id))
      index = remote_url.lastIndexOf(".")
      if(index >=0):
        file_type = remote_url.right(remote_url.size()-index-1)
    if file_type is None:
      file_type = self.file_type_
    if len(self.name_)>0:
      return url+formatted_id +"_"+self.name_+"."+file_type
    else:
      return url+formatted_id+"."+file_type
  
  def HandleError(self, message):
    ost.LogError(str("Could not download file\n"+message+""))
   
  def ToInfo(self,group):
    group.SetAttribute(GenericLoader.EXT_NAME_ATTRIBUTE_NAME, str(self.name_))
    group.SetAttribute(GenericLoader.URL_ATTRIBUTE_NAME, str(self.url_))
    if self.up_case_:
      group.SetAttribute(GenericLoader.UP_CASE_ATTRIBUTE_NAME, str(self.up_case_))
    if self.file_type_ is not None:
      group.SetAttribute(GenericLoader.FILE_TYPE_ATTRIBUTE_NAME, str(self.file_type_))
    if self.tmp_url_ is not None:
      group.SetAttribute(GenericLoader.TMP_URL_ATTRIBUTE_NAME, str(self.tmp_url_))
    if self.img_ is not None:
      group.SetAttribute(GenericLoader.IMG_ATTRIBUTE_NAME, str(int(self.img_)))
    if self.default_ is not None:
      group.SetAttribute(GenericLoader.DEFAULT_ATTRIBUTE_NAME, str(int(self.default_)))
      
  @staticmethod
  def FromInfo(group):
    name = None
    url = None
    up_case = False
    file_type = None
    tmp_url = None
    img = False
    default = False
    
    if group.HasAttribute(GenericLoader.EXT_NAME_ATTRIBUTE_NAME):
      name = QtCore.QString(group.GetAttribute(GenericLoader.EXT_NAME_ATTRIBUTE_NAME))
    
    if group.HasAttribute(GenericLoader.URL_ATTRIBUTE_NAME):
      url = QtCore.QString(group.GetAttribute(GenericLoader.URL_ATTRIBUTE_NAME))
      
    if group.HasAttribute(GenericLoader.UP_CASE_ATTRIBUTE_NAME):
      up_case = bool(group.GetAttribute(GenericLoader.UP_CASE_ATTRIBUTE_NAME))  
    
    if group.HasAttribute(GenericLoader.FILE_TYPE_ATTRIBUTE_NAME):
      file_type = QtCore.QString(group.GetAttribute(GenericLoader.FILE_TYPE_ATTRIBUTE_NAME))
      
    if group.HasAttribute(GenericLoader.TMP_URL_ATTRIBUTE_NAME):
      tmp_url = QtCore.QString(group.GetAttribute(GenericLoader.TMP_URL_ATTRIBUTE_NAME))
    
    if group.HasAttribute(GenericLoader.IMG_ATTRIBUTE_NAME):
      img = bool(int(group.GetAttribute(GenericLoader.IMG_ATTRIBUTE_NAME)))
      
    if group.HasAttribute(GenericLoader.DEFAULT_ATTRIBUTE_NAME):
      default = bool(int(group.GetAttribute(GenericLoader.DEFAULT_ATTRIBUTE_NAME)))

    return GenericLoader(name, url, up_case, file_type, tmp_url, img, default)
"""
