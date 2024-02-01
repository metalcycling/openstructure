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
#ifndef OST_SEQ_SEQUENCE_IDENTITY_HH
#define OST_SEQ_SEQUENCE_IDENTITY_HH

/*  
  Author: Marco Biasini
 */
#include <ost/seq/alg/module_config.hh>
#include <ost/seq/alignment_handle.hh>

namespace ost { namespace seq { namespace alg {

struct RefMode {
  enum Type {
    ALIGNMENT=0,
    LONGER_SEQUENCE
  };
};
  
/// \brief calculate sequence identity for two sequences in an alignment
/// \param ref_mode influences the way the sequence identity is calculated. When 
///     set to LONGER_SEQUENCE, the sequence identity is calculated as the 
///     number of matches  divided by the length of the longer sequence. If
///     set to ALIGNMENT, the sequence identity is calculated as the number of
///     matches divided by the number of aligned residues.
/// \param seq_a is the index of the first sequence
/// \param seq_b is the index of the second sequence
/// \param aln is the sequence alignment
/// \return sequence identity in the range 0 to 100
Real DLLEXPORT_OST_SEQ_ALG
SequenceIdentity(const AlignmentHandle& aln, 
                 RefMode::Type ref_mode=RefMode::LONGER_SEQUENCE, 
                 int seq_a=0, int seq_b=1);

}}}

#endif
