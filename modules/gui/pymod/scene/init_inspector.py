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
from ost import gfx
from PyQt5 import QtCore
from ost.gui.scene.inspector_widget import InspectorDialog
  
def _InitInspector(app):
  mywidget = InspectorDialog(app.gl_win.qobject)
  mywidget.setWindowFlags(QtCore.Qt.Dialog | QtCore.Qt.Tool)
  mywidget.show()
  app.AddWidgetToApp("InspectorDialog", mywidget)  
  menu_bar=app.perspective.GetMenuBar()

