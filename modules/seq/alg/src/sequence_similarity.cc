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

#include <ost/seq/aligned_column_iterator.hh> 
#include <ost/seq/alg/sequence_similarity.hh>


namespace ost { namespace seq { namespace alg {

Real SequenceSimilarity(const AlignmentHandle& aln, SubstWeightMatrixPtr subst, 
                        bool normalize, int seq_a, int seq_b){
  int non_gap_count=0;
  int similarity_score = 0;

  const String& sa=aln.GetSequence(seq_a).GetString();
  const String& sb=aln.GetSequence(seq_b).GetString();  

  for (String::const_iterator 
       i=sa.begin(), e=sa.end(), j=sb.begin(); i!=e; ++i, ++j) {
    if (!((*i)== '-' || (*j)=='-')) {
      ++non_gap_count;
      similarity_score += subst->GetWeight(*i,*j);
    }
  }

  Real seq_sim = subst->GetMinWeight();

  if(non_gap_count > 0){
    seq_sim = static_cast<Real>(similarity_score)/non_gap_count;
  }

  if(normalize){
    seq_sim -= subst->GetMinWeight();
    Real weight_range = subst->GetMaxWeight() - subst->GetMinWeight();
    if(weight_range > 0) seq_sim /= weight_range;
  }

  return seq_sim;
}

}}}
