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
import traceback
from ._ost_gfx import *

class PyGfxObj(GfxObj):
  def __init__(self,name):
    """
    requires a unique name not yet utilized in the Scene;
    do not place OpenGL calls in the ctor, use InitGL for
    that purpose
    """
    GfxObj.__init__(self,name)
    self._valid_flag=False

  def _InitGL(self):
    try:
      self.InitGL()
      self._valid_flag=True
    except:
      traceback.print_exc()
      
  def InitGL(self):
    """
    called once for each OpenGL context (usually one time),
    allows one-time OpenGL initialization to be implemented,
    such as vbo allocation
    """
    pass


  def _CustomPreRenderGL(self,rebuild):
    if not self._valid_flag:
      return
    try:
      self.CustomPreRenderGL(rebuild)
    except:
      self._valid_flag=False
      traceback.print_exc()

  def CustomPreRenderGL(self,rebuild):
    """
    called just before CustomRenderGL is called; the flag
    indicates that a rebuild is required or was requested
    """
    pass


  def _CustomRenderGL(self,render_pass):
    if not self._valid_flag:
      return
    try:
      self.CustomRenderGL(render_pass)
    except:
      self._valid_flag=False
      traceback.print_exc()

  def CustomRenderGL(self,render_pass):
    """
    called for each scene refresh
    """
    pass

