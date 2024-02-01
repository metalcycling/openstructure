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
from ._ost_seq import *
from ost.seq import alg

_CreateAlignment=CreateAlignment
def CreateAlignment(*seqs):
  aln=_CreateAlignment()
  for s in seqs:
    aln.AddSequence(s)
  return aln
  
_CreateSequenceList=CreateSequenceList
def CreateSequenceList(*seqs):
  sl=_CreateSequenceList()
  for s in seqs:
    sl.AddSequence(s)
  return sl

