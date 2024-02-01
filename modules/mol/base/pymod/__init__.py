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
from ._ost_mol import *
import ost.geom as _geom
from ost.mol import alg


def Transform(tf=None):
  from ost import LogWarning
  if Transform.mol_transform_warning_flag:
    LogWarning("mol.Transform is deprecated, please use geom.Transform instead")
    Transform.mol_transform_warning_flag=False
  if tf:
    return _geom.Transform(tf)
  else:
    return _geom.Transform()
Transform.mol_transform_warning_flag=True

def MergeCoordGroups(*coord_groups):
  """
  Merge several separate coord groups into one. The coord groups must have the 
  same number of atoms. In case no coord group is supplied, None will be 
  returned.
  """
  if len(coord_groups)==0:
    return None
  cg=CreateCoordGroup(coord_groups[0].atoms)
  for coord_group in coord_groups:
    cg.AddFrames(coord_group)
  return cg
