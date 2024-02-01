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
#ifndef OST_SEQ_ALG_ALIGNMENT_OPTS_HH
#define OST_SEQ_ALG_ALIGNMENT_OPTS_HH

/*
  Author: Marco Biasini
 */
#include <ost/seq/alg/module_config.hh>
#include <ost/seq/alg/subst_weight_matrix.hh>
namespace ost { namespace seq { namespace alg {
  
/// \brief options for local and global sequence alignment algorithms
///
/// \sa LocalSequenceAlign, GlobalSequenceAlign
struct DLLEXPORT_OST_SEQ_ALG AlignmentOpts {
  AlignmentOpts(int gap_open, int gap_ext, 
                const SubstWeightMatrixPtr& weights):
    gap_open_penalty(gap_open), gap_extension_penalty(gap_ext),
    subst_weights(weights)
  {
  }
  int gap_open_penalty;
  int gap_extension_penalty;
  SubstWeightMatrixPtr subst_weights;
};
}}}
#endif
