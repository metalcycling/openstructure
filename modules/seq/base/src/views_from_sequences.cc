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
#include <ost/integrity_error.hh>
#include <ost/seq/views_from_sequences.hh>

namespace ost { namespace seq {

namespace {

bool skip_gaps(const ConstSequenceHandle& seq1, 
               const ConstSequenceHandle& seq2, int& pos, 
               int& index_a, int& index_b)
{
  while (seq1.GetLength()>pos) {
    if (seq1[pos]=='-') {
      if (seq2[pos]!='-') {
        index_b+=1;
      }
    } else if (seq2[pos]=='-') {
      index_a+=1;
    } else {
      return true;
    }
    pos+=1;
  }
  return false;
}

}

std::pair<mol::EntityView, mol::EntityView> 
ViewsFromSequences(const ConstSequenceHandle& seq1, 
                   const ConstSequenceHandle& seq2)
{
  if (seq1.GetLength()!=seq2.GetLength()) {
    throw IntegrityError("Sequence lengths do not match");
  }
  if (!seq1.HasAttachedView()) {
    throw IntegrityError("First sequence does not have an attached view");
  }
  if (!seq2.HasAttachedView()) {
    throw IntegrityError("Second sequence does not have an attached view");
  }
  mol::EntityView src_a=seq1.GetAttachedView();
  mol::EntityView src_b=seq2.GetAttachedView();
  
  mol::EntityView dst_a=src_a.CreateEmptyView();
  mol::EntityView dst_b=src_b.CreateEmptyView();
  mol::ResidueViewList res_a=src_a.GetResidueList();
  mol::ResidueViewList res_b=src_b.GetResidueList();  
  int pos=0, index_a=0, index_b=0;
  while (skip_gaps(seq1, seq2, pos, index_a, index_b)) {
    dst_a.AddResidue(res_a.at(seq1.GetOffset()+index_a), 
                     mol::ViewAddFlag::INCLUDE_ATOMS);
    dst_b.AddResidue(res_b.at(seq2.GetOffset()+index_b), 
                     mol::ViewAddFlag::INCLUDE_ATOMS);
    pos+=1;
    index_a+=1;
    index_b+=1;
  }
  return std::make_pair(dst_a, dst_b);
}


std::pair<mol::EntityView, mol::EntityView> 
ViewsFromAlignment(const AlignmentHandle& aln, int index1, int index2)
{
  return ViewsFromSequences(aln.GetSequence(index1), aln.GetSequence(index2));
}

}}
