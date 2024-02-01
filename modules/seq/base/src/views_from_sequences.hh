//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 3.0 of the License, or (at your option)
// any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//------------------------------------------------------------------------------
#ifndef OST_SEQ_VIEWS_FROM_SEQUENCES_HH
#define OST_SEQ_VIEWS_FROM_SEQUENCES_HH

#include <ost/seq/module_config.hh>
#include <ost/seq/alignment_handle.hh>

namespace ost { namespace seq {

std::pair<mol::EntityView, mol::EntityView> DLLEXPORT_OST_SEQ 
ViewsFromSequences(const ConstSequenceHandle& seq1, 
                   const ConstSequenceHandle& seq2);

std::pair<mol::EntityView, mol::EntityView> DLLEXPORT_OST_SEQ
ViewsFromAlignment(const AlignmentHandle& aln, int index1=0, int index2=1);
  
}}

#endif
