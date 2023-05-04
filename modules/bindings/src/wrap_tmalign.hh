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
#ifndef OST_BINDINGS_TMALIGN_H
#define OST_BINDINGS_TMALIGN_H

#include <ost/geom/mat4.hh>
#include <ost/geom/vec3.hh>
#include <ost/seq/alignment_handle.hh>
#include <ost/mol/chain_view.hh>
#include <ost/mol/chain_handle.hh>

namespace ost { namespace bindings {

struct TMAlignResult {

  TMAlignResult() { }

  TMAlignResult(Real rm, Real tm, int aln_l, const geom::Mat4& t, 
                const ost::seq::AlignmentHandle& aln): rmsd(rm),
                                                       tm_score(tm),
                                                       aligned_length(aln_l),
                                                       transform(t),
                                                       alignment(aln) { } 


  Real rmsd;
  Real tm_score;
  int aligned_length;
  geom::Mat4 transform;
  ost::seq::AlignmentHandle alignment;

  Real GetTMScore() { return tm_score; }
  Real GetRMSD() { return rmsd; }
  int GetAlignedLength() { return aligned_length; }
  const geom::Mat4& GetTransform() { return transform; }
  const ost::seq::AlignmentHandle& GetAlignment() { return alignment; }
};

TMAlignResult WrappedTMAlign(const geom::Vec3List& pos_one, 
                             const geom::Vec3List& pos_two, 
                             const ost::seq::SequenceHandle& seq1,
                             const ost::seq::SequenceHandle& seq2,
                             bool fast = false,
                             bool rna = false);

TMAlignResult WrappedTMAlign(const ost::mol::ChainView& ent1,
                             const ost::mol::ChainView& ent2,
                             bool fast = false);
}} //ns

#endif
