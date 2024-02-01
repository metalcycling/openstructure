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

#include <ost/seq/alignment_handle.hh>
#include <ost/seq/alg/module_config.hh>

/*
  Authors: Marco Biasini
 */
namespace ost { namespace seq { namespace alg {

/// \brief Calculates conservation scores for each column in the alignment. 
/// 
/// The conservation score is a value between 0 and 1. The bigger the number 
/// the more conserved the aligned residues are. 
/// 
/// \p assign If true, the conservation scores are assigned to attached 
///     residues. The name of the property can be changed with the prop_name 
///     parameter. Useful when coloring entities based on sequence conservation.
/// \p prop_name The property name for assigning the conservation to 
///     attached residues. Defaults to 'cons'.
/// \p ignore_gap If true, the dissimilarity between two gaps is increased to
///     6.0 instead of 0.5 as defined in the original version. Without this, a
///     stretch where in the alignment there is only one sequence which is
///     aligned to only gaps, is considered highly conserved (depending on the
///     number of gap sequences).
/// 
std::vector<Real> DLLEXPORT_OST_SEQ_ALG Conservation(const AlignmentHandle& aln,
                                             bool assign=true,
                                             const String& prop_name="cons",
                                             bool ignore_gap=false);
}}}

