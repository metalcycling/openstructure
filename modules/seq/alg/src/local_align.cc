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
#include <ost/log.hh>
#include "impl/align_impl.hh"
#include "local_align.hh"

namespace ost { namespace seq { namespace alg {

bool CollectPatches(impl::AlnMat& mat, int min_i, int max_i, int min_j, int max_j,
                    const ConstSequenceHandle& s1, 
                    const ConstSequenceHandle& s2,
                    AlignmentList& alignments)
{
  int max_val=0, mmax_i=0, mmax_j=0;
  for (int i=min_i; i<max_i; ++i) {
    for (int j=min_j; j<max_j; ++j) {
      if (mat(i, j).score>=max_val) {
        mmax_i=i;
        mmax_j=j;
        max_val=mat(i, j).score;
      }
    }
  }
  int i=mmax_i; int j=mmax_j;
  String aln_str1;
  String aln_str2;
  while (i>min_i && j>min_j && mat(i, j).score>0) {
    impl::SetRoute(mat, i, j, s1, s2, aln_str1, aln_str2);
  }
  if (aln_str1.size()<2) {
    return false;
  }
  impl::StoreStrAsAln(aln_str1, aln_str2, s1, s2, i, j, alignments);
  for (int x=i; x<=mmax_i; ++x) {
    for (int y=0; y<mat.GetHeight(); ++y) {
      mat(x, y).score=0;
    }
  }
  for (int x=0; x<mat.GetWidth(); ++x) {
    for (int y=j; y<=mmax_j; ++y) {
      mat(x, y).score=0;
    }
  }  
  return true;
}

AlignmentList LocalAlign(const ConstSequenceHandle& s1, 
                         const ConstSequenceHandle& s2,
                         alg::SubstWeightMatrixPtr& subst,
                         int gap_open, int gap_ext)
{
  impl::AlnMat mat(s1.GetLength()+1, s2.GetLength()+1);
  for (int i=0; i<mat.GetWidth()-1; ++i) {
    for (int j=0; j<mat.GetHeight()-1; ++j) {
      char c1=s1[i];
      char c2=s2[j];
      short weight=subst->GetWeight(c1, c2);
      int diag=weight+mat(i, j).score;
      int ins1=mat(i+1, j).score
              +(mat(i+1, j).from==impl::INS1 ? gap_ext : gap_open);
      ins1=std::max(0, ins1);
      int ins2=mat(i, j+1).score
              +(mat(i, j+1).from==impl::INS2 ? gap_ext : gap_open);
      ins2=std::max(0, ins2);
      if (diag>=ins1) {
        if (diag>=ins2) {
          mat(i+1, j+1).score=diag;
          mat(i+1, j+1).from=impl::DIAG;
        } else {
          mat(i+1, j+1).score=ins2;
          mat(i+1, j+1).from=impl::INS2;
        }
      } else if (ins1>ins2) {
        mat(i+1, j+1).score=ins1;
        mat(i+1, j+1).from=impl::INS1;
      } else {
        mat(i+1, j+1).score=ins2;
        mat(i+1, j+1).from=impl::INS2;
      }
    }
  }
#if !defined(NDEBUG)
  DbgWriteAlnMatrix(mat,s1, s2);
#endif
  AlignmentList alignments;
  while (CollectPatches(mat, 0, mat.GetWidth(), 0, mat.GetHeight(), 
         s1, s2, alignments)) {
    LOG_DEBUG(alignments.back().ToString(80));
  }
  return alignments;
}

}}}
