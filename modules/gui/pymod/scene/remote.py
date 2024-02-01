from PyQt5 import QtCore, QtWidgets
import threading
from ost import LogError
from ost import gui, gfx
from ost.io.remote import RemoteLoad, REMOTE_REPOSITORIES
import re

class RemoteLoader(QtWidgets.QWidget):
  def __init__(self):
    QtWidgets.QWidget.__init__(self)
    self._line = QtWidgets.QLineEdit(self)
    self._load = QtWidgets.QToolButton(self)
    self._line.returnPressed.connect(self._DoLoad)
    self._load.clicked.connect(self._DoLoad)
    hbox = QtWidgets.QHBoxLayout(self)
    hbox.addWidget(self._line,1)
    hbox.addWidget(self._load,0)
    self.setLayout(hbox)
    hbox.setContentsMargins(3,3,3,3)
    hbox.setSpacing(3)
    self._menu = self._RemoteMenu()
    self._load.setMenu(self._menu)
    self._current_repo = 'pdb'
    self._UpdateLoadButton(self._current_repo)
    self.setFixedHeight(self._load.height())
    self.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)

  def _RemoteMenu(self):
    menu = QtWidgets.QMenu()
    action_group = QtWidgets.QActionGroup(menu)
    for k,v in REMOTE_REPOSITORIES.items():
      action = menu.addAction(v.name)
      action.setCheckable(True)
      if k == 'pdb':
        action.setChecked(True)
      action.setData(QtCore.QVariant(k))
      action_group.addAction(action)
    menu.triggered.connect(self._ToggleRepo)
    return menu
  def _UpdateLoadButton(self, current):
    name = REMOTE_REPOSITORIES[current].name
    self._load.setText('load from %s  ' % name)

  def _ToggleRepo(self, action):
    self._current_repo = str(action.data())
    self._UpdateLoadButton(self._current_repo)


  def _DoLoad(self):
    split_ids = []
    for p in str(self._line.text()).split(','):
      for k in p.split():
        if len(k.strip()):
          split_ids.append(k.strip())
    self._line.setText('')
    for split_id in split_ids:
      try:
        ent = RemoteLoad(split_id, from_repo=self._current_repo)
      except Exception as e:
        LogError(str(e))
        continue
      g = gfx.Entity(split_id, ent)
      try:
        gfx.Scene().Add(g)
      except Exception as e:
        LogError(str(e))

remote_loader=RemoteLoader()
remote_loader_for_panel=gui.Widget(remote_loader)
remote_loader_for_panel.qobject.setSizePolicy(QtWidgets.QSizePolicy.Expanding, 
                                              QtWidgets.QSizePolicy.Fixed)
panels=gui.GostyApp.Instance().perspective.panels
panels.AddWidgetToPool("Remote Loader",remote_loader_for_panel)
