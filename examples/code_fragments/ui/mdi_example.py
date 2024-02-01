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
from PyQt5 import QtGui, QtCore, QtWidgets

from ost import gui
import sys, random

class Points(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowTitle('Some Points')

    def paintEvent(self, event):
        paint = QtGui.QPainter()
        paint.begin(self)
        size = self.size()
        paint.setPen(QtCore.Qt.red)
        for i in range(5000):
            x = random.randint(1, size.width()-1)
            y = random.randint(1, size.height()-1)
            paint.drawPoint(x, y)
        paint.setPen(QtCore.Qt.green)
        for i in range(5000):
            x = random.randint(1, size.width()-1)
            y = random.randint(1, size.height()-1)
            paint.drawPoint(x, y)
        paint.setPen(QtCore.Qt.blue)
        for i in range(5000):
            x = random.randint(1, size.width()-1)
            y = random.randint(1, size.height()-1)
            paint.drawPoint(x, y)
        paint.end()


app=gui.GostyApp.Instance()

#Get main area widget
main_area=app.perspective.main_area

pts=Points(main_area.qobject)

#Add Widget
main_area.AddWidget("Some Points", pts)


