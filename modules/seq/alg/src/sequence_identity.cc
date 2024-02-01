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
#include <ost/seq/alg/sequence_identity.hh>

namespace ost { namespace seq { namespace alg {

Real SequenceIdentity(const AlignmentHandle& aln, RefMode::Type ref_mode, 
                        int seq_a, int seq_b)
{
  int non_gap_count=0;
  int identical_count=0;
  const String& sa=aln.GetSequence(seq_a).GetString();
  const String& sb=aln.GetSequence(seq_b).GetString();  
  for (String::const_iterator 
       i=sa.begin(), e=sa.end(), j=sb.begin(); i!=e; ++i, ++j) {
    if (!((*i)== '-' || (*j)=='-')) {
      non_gap_count++;
    }
    if ((*i)!='-' && (*i)==(*j)) {
      identical_count++;
    }
  }
  Real seq_id=0.0;
  // divide by length of longer sequence:
  if (ref_mode==RefMode::LONGER_SEQUENCE) {
    int wo_gaps_a=aln.GetSequence(seq_a).GetGaplessString().length();
    int wo_gaps_b=aln.GetSequence(seq_b).GetGaplessString().length();
    if(wo_gaps_a>=wo_gaps_b) {
      if (wo_gaps_a==0) {
        seq_id=0;
      } else {
        seq_id=100*static_cast<Real>(identical_count)/wo_gaps_a;        
      }

    } else {
      if (wo_gaps_b==0) {
        seq_id=0;
      } else {
        seq_id=100*static_cast<Real>(identical_count)/wo_gaps_b;        
      }
    }
  } else {
    if (non_gap_count==0) {
      seq_id=0;
    } else {
      // divide by alignment length:
      seq_id=100*static_cast<Real>(identical_count)/non_gap_count;      
    }
  }
  return seq_id;
}

}}}
