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
#ifndef OST_SEQ_SEQUENCE_SIMILARITY_HH
#define OST_SEQ_SEQUENCE_SIMILARITY_HH

#include <ost/seq/alg/module_config.hh>
#include <ost/seq/alignment_handle.hh>
#include <ost/seq/alg/subst_weight_matrix.hh>

namespace ost { namespace seq { namespace alg {
  
/// \brief calculate sequence similarity for two sequences in an alignment
/// \param seq_a is the index of the first sequence
/// \param seq_b is the index of the second sequence
/// \param aln is the sequence alignment
/// \return sequence similarity
Real DLLEXPORT_OST_SEQ_ALG
SequenceSimilarity(const AlignmentHandle& aln, SubstWeightMatrixPtr subst,
                   bool normalize=false, int seq_a=0, int seq_b=1);

}}}

#endif
