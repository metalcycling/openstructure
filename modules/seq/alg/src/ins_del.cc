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
  Author: Marco Biasini
 */
#include <ost/seq/aligned_column_iterator.hh>
#include <ost/seq/alg/ins_del.hh>

namespace ost { namespace seq { namespace alg {

InsDel::InsDel(const AlignmentHandle& aln, int seq_a, int seq_b):
  aln_(aln), seq_a_(seq_a), seq_b_(seq_b)
{
  AlignedColumnIterator i=aln.begin(), e=aln.end();
  while (i!=e) {
    AlignedColumn col=*i;
    int start=col.GetIndex();
    while (i!=e && (*i)[seq_a]=='-' && (*i)[seq_b]!='-') {
      ++i;
    }
    if (start!=i->GetIndex()) {
      insdels_.push_back(AlignedRegion(aln, start, i->GetIndex()));
      continue;
    }
    while (i!=e && (*i)[seq_b]=='-' && (*i)[seq_a]!='-') {
      ++i;
    }
    if (start!=i->GetIndex()) {
      insdels_.push_back(AlignedRegion(aln, start, i->GetIndex()));
      continue;
    }
    ++i;      
  }
}


AlignedRegionList InsDel::GetInsertions(int ref) const
{
  AlignedRegionList ins;
  for (AlignedRegionList::const_iterator i=insdels_.begin(), 
    e=insdels_.end(); i!=e; ++i) {
    AlignedRegion rgn=*i;
    if (rgn[0][ref]!='-') {
      rgn.SetMaster(ref);       
      ins.push_back(rgn);
    }    
  }
  return ins;
}

AlignedRegionList InsDel::GetDeletions(int ref) const
{
  AlignedRegionList del;
  for (AlignedRegionList::const_iterator i=insdels_.begin(), 
    e=insdels_.end(); i!=e; ++i) {
    AlignedRegion rgn=*i;
    if (rgn[0][ref]=='-') {
      rgn.SetMaster(ref);
      del.push_back(rgn);
    }    
  }
  return del;
}

}}}
