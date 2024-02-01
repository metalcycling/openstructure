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

from PyQt5 import QtCore, QtGui, QtWidgets

#Get Panels (Class which manages widgets)
panels=gui.GostyApp.Instance().perspective.panels

#Create Widget
my_awesome_widget = QtWidgets.QWidget()
layout = QtWidgets.QVBoxLayout()
layout.addWidget(QtWidgets.QPushButton('I dont do anything'))
layout.addWidget(QtWidgets.QPushButton('I do even less'))
my_awesome_widget.setLayout(layout)

#Wrap widget to Qt Widget
wid=gui.Widget(my_awesome_widget)

#Add Widget to widget pool
panels.AddWidgetToPool("Widget that doesnt do a thing", wid)

#Add Widget to right panel
panels.AddWidget(gui.PanelPosition.RIGHT_PANEL,wid,False)

