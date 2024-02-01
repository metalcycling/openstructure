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
#ifndef OST_SEQ_ALG_INS_DEL_HH
#define OST_SEQ_ALG_INS_DEL_HH

/*  
  Author: Marco Biasini
 */
#include <ost/seq/aligned_region.hh>
#include <ost/seq/alignment_handle.hh>
#include <ost/seq/alg/module_config.hh>

namespace ost { namespace seq { namespace alg {

/// \brief get insertions and deletions of an alignment
///
/// InsDel collects regions for which one of the two sequences contains
/// a gap. Whether such a region is an insertion or a deletion depends on 
/// the perspective. For example, given
/// \verbatim
/// abcXYZdef
/// abc---def
/// \endverbatim
/// the half-closed interval 3:6 is an insertion for the first sequence, but a 
/// deletion for the second sequence. Thus calling GetInsertions() and 
/// GetDeletions() requires you to pass in a sequence of reference.
class DLLEXPORT_OST_SEQ_ALG InsDel {
public:
  /// \brief construct new ins del instance
  InsDel(const AlignmentHandle& aln, int seq_a=0, int seq_b=1);
  /// \brief get regions that are insertions of the reference sequence
  AlignedRegionList GetInsertions(int ref) const;
  /// \brief get regiosn that are deletions for the reference sequence
  AlignedRegionList GetDeletions(int ref) const;
private:
  AlignmentHandle    aln_;
  int                seq_a_;
  int                seq_b_;
  AlignedRegionList  insdels_;
};

}}}

#endif
