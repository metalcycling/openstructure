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
import ost.gui.trajectory_viewer
from ._ost_gui import *

## \brief Opens a DataViewer 
# \sa \example fft_li.py "View Fourier Transform Example" \sa \ref modulate_image.py "Modulate Image Example"
def _close_event_override_(event):
  print("close event")
def _set_data_override_(data):
  print("set data")

def CreateDataViewer(ih,flag=False):
      viewer=GostyApp.Instance().CreateDataViewer(ih)
      if flag:
        viewer.image_=ih
        sip_viewer=viewer.qobject
        sip_viewer.closeEvent=_close_event_override_
        sip_viewer.setData=_set_data_override_
      return viewer

def ClearMessageWidget():
      gosty=GostyApp.Instance()
      gosty.message_widget.Clear()


from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from ost import gfx

def PickColor(default=gfx.WHITE):
  """
  Pops up a color chooser that lets' the user pick a color and returns the 
  selected color. If the user cancels the color chooser, None is returned.
  
  :rtype: :class:`~ost.gfx.Color`
  """
  dialog=QColorDialog()
  qt_color=QColor(int(min(255,default.Red()*256)),
                  int(min(255,default.Green()*256)),
                  int(min(255,default.Blue()*256)))
  qt_color=dialog.getColor(qt_color)
  if not qt_color.isValid():
    return None
  return gfx.RGBb(qt_color.red(), qt_color.green(),qt_color.blue())
                   
def GetMenu(menu_name, create=False):
  persp=GostyApp.Instance().perspective
  if isinstance(menu_name[0], QAction):
    return menu_name[0]
  else:
    node=persp.GetMenu(menu_name[0])
    for name in menu_name[1:]:
      found=False
      for action in node.actions():
        if str(action.text())==str(name):
          found=True
          node=action
          break
      if not found:
        if create:
          node=node.addMenu(name)  
        else:
          raise ValueError("Menu '%s' doesn't exist" % ', '.join(menu_name))
    return node



def AddMenuAction(*args, **kwargs):
  """
  Add menu action to main menu.

  This function lets you conveniently add actions to the main menu. To add a new 
  new action "Background Color" to the "Scene" menu, simply use
  
  .. code-block:: python
  
    def SetBackgroundColor():
      scene.bg=gfx.PickColor(scene.bg)

    AddMenuAction('Scene', "Background Color", SetBackgroundColor)
    
  This will add the menu "Scene" if it does not exist yet, register the action 
  "Background Color" and execute the function SetBackgroundColor whenever the 
  action is triggered.
  
  To assign a keyboard shortcut to the action, you can use the shortcut 
  argument:
  
  .. code-block:: python
  
    AddMenuAction('Scene', 'Background Color', SetBackgroundColor,
                  shortcut='Ctrl+B')

  Whenever you press Ctrl+B (Cmd+B on macOS), the action will be executed.
  
  Very often menu actions are coupled to the current selected objects in the 
  scene menu. These menu actions are either enabled or disabled depending on the 
  type of the selected objects. To easily support this scenario, the "enable" 
  state of the menu action can be tightly coupled to the scene menu by providing 
  a callable to the enabled argument. As an example, the following menu action 
  is only enabled when exactly one gfx.Entity is selected.
  
  .. code-block:: python
  
    AddMenuAction('Scene', 'PrintAtomCount', PrintAtomCount,
                  enabled=OneOf(gfx.Entity))
  
  OneOf is a simple function object that takes any number of classes and returns 
  true when the selected object is an instance. Alterantively, you might want to 
  use ManyOf or supply a custom function that evaluates the state of the menu 
  action to suit your needs.
  """
  class MenuActionEnabler(QObject):
    def __init__(self, predicate, action):
      QObject.__init__(self, action)
      self.predicate=predicate
      self.action=action
      app=GostyApp.Instance()
      app.scene_win.qobject.ActiveNodesChanged.connect(self.TestEnable)
      self.TestEnable()

    def TestEnable(self):
      self.action.setEnabled(self.predicate())
  persp=GostyApp.Instance().perspective
  menu_name=args[:-1]
  function=args[-1]
  if isinstance(menu_name[0], QMenu):
    node=menu_name[0]
  else:
    node=persp.GetMenu(args[0])
    for name in menu_name[1:-1]:
      found=False
      for action in node.actions():
        if str(action.text())==str(name):
          node=action
          break
      if not found:
        node=node.addMenu(name)
  action=node.addAction(str(menu_name[-1]))
  if 'shortcut' in kwargs:
    action.setShortcut(QKeySequence(kwargs['shortcut']))
  if 'checkable' in kwargs:
    action.setCheckable(kwargs['checkable'])
  if 'checked' in kwargs:
    action.setChecked(kwargs['checked'])
  if 'enabled' in kwargs:
    if callable(kwargs['enabled']):
      enabler=MenuActionEnabler(kwargs['enabled'], action)
    else:
      action.setEnabled(kwargs['enabled'])
  action.triggered.connect(function)
  return action


class OneOf:
  def __init__(self, *classes):
    self.classes=classes
  def __call__(self):
    sel=SceneSelection.Instance()
    if sel.GetActiveNodeCount()!=1:
      return False
    node=sel.GetActiveNode(0)
    for cl in self.classes:
      if isinstance(node, cl):
        return True
    return False

class TwoOf:
  def __init__(self, *classes):
    self.classes=classes
  def __call__(self):
    sel=SceneSelection.Instance()
    act_count=sel.GetActiveNodeCount()
    if act_count<2:
      return False
    found=0
    for i in range(0, act_count):
      node=sel.GetActiveNode(i)
      for cl in self.classes:
        if isinstance(node, cl):
          found += 1
      if found > 2:
        return False
    if found == 2:
      return True
    return False

class ManyOf:
  def __init__(self, *classes):
    self.classes=classes
  def __call__(self):
    sel=SceneSelection.Instance()
    if sel.GetActiveNodeCount()==0:
      return False
    for  i in range(sel.GetActiveNodeCount()):
      node=sel.GetActiveNode(i)
      found=False
      for cl in self.classes:
        if isinstance(node, cl):
          found=True
          break
      if not found:
        return False
    return True

from ost import PushVerbosityLevel as _PushVerbosityLevel
from ost import PopVerbosityLevel as _PopVerbosityLevel
from ost import GetVerbosityLevel as _GetVerbosityLevel


def PushVerbosityLevel(value):
  GostyApp.Instance().perspective.ChangeVerbositySlider(value)

def PopVerbosityLevel():
  _PopVerbosityLevel()
  GostyApp.Instance().perspective.ChangeVerbositySlider(_GetVerbosityLevel())
  _PopVerbosityLevel() # the slider change pushes a new level :-(
