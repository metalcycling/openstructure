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
#ifndef OST_SEQ_ALG_MERGE_PAIRWISE_ALIGNMENTS_HH
#define OST_SEQ_ALG_MERGE_PAIRWISE_ALIGNMENTS_HH

#include <ost/seq/alignment_handle.hh>
#include <ost/seq/alg/module_config.hh>
/*
  Author: Marco Biasini
 */
namespace ost { namespace seq { namespace alg {


/// \brief merge a list of pairwise alignments into one multiple sequence 
///     alignment
/// 
/// All sequences in the pairwise sequence alignments are a realigned to the 
/// reference sequence. This is useful to merge the results of a BLAST or HMM 
/// database search into one multiple sequence alignment.
/// 
/// The method does not produce the optimal multiple sequence alignemnt for all 
/// the sequences.
/// 
/// \param pairwise_alns is a list of AlignmentHandles, each containing
///      two sequences
/// \param ref_seq is the reference sequence. The reference sequence must not
///      contain any gaps.
AlignmentHandle DLLEXPORT_OST_SEQ_ALG 
MergePairwiseAlignments(const AlignmentList& pairwise_alns,
                        const ConstSequenceHandle& ref_seq);
}}}

#endif
