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
from ._ost_conop import *

# The 20 standard amino acids in no particular order
STANDARD_AMINOACIDS=(
 'ALA', 'ARG', 'ASN',
 'ASP', 'GLN', 'GLU',
 'LYS', 'SER', 'CYS',
 'MET', 'TRP', 'TYR',
 'THR', 'VAL', 'ILE',
 'LEU', 'GLY', 'PRO',
 'HIS', 'PHE',
)

def SetDefaultLib(compound_lib):
  '''
  Set the default compound library. The compound library is used by various
  functions of the framework that requires knowledge of naming and 
  connectivity of residues.
  '''
  conop_inst=Conopology.Instance()
  conop_inst.SetDefaultLib(compound_lib)

def GetDefaultLib():
  conop_inst=Conopology.Instance()
  return conop_inst.GetDefaultLib()

