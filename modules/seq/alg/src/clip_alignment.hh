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

/*
  Author: Gerardo Tauriello, Juergen Haas
 */

#ifndef OST_SEQ_ALG_CLIP_ALIGNMENT_HH
#define OST_SEQ_ALG_CLIP_ALIGNMENT_HH

#include <ost/seq/alignment_handle.hh>
#include <ost/seq/alg/module_config.hh>

namespace ost { namespace seq { namespace alg {

/// \brief Clips alignment so that first and last column have at least the
///        desired number of structures.
/// 
/// \param aln Multiple sequence alignment. Will be cut.
/// \param n_seq_thresh Minimal number of sequences desired.
/// \param set_offset Shall we update offsets for attached views?
/// \param remove_empty Shall we remove sequences with only gaps in cut aln?
/// \returns Starting column (0-indexed), where cut region starts (w.r.t.
///          original aln). -1, if there is no region in the alignment with
///          at least the desired number of structures.
int DLLEXPORT_OST_SEQ_ALG
ClipAlignment(AlignmentHandle& aln, uint n_seq_thresh=2, bool set_offset=true,
              bool remove_empty=true);
}}}

#endif
