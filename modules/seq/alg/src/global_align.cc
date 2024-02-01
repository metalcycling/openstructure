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
#include "global_align.hh"

namespace ost { namespace seq { namespace alg {

void Traceback(impl::AlnMat& mat, int max_i, int max_j,
               const ConstSequenceHandle& s1, 
               const ConstSequenceHandle& s2,
               AlignmentList& alignments)
{
  int i = (max_i > 0 ? max_i-1 : 0);
  int j = (max_j > 0 ? max_j-1 : 0);
  String aln_str1;
  String aln_str2;
  while (i > 0 || j > 0) {
    impl::SetRoute(mat, i, j, s1, s2, aln_str1, aln_str2);
  }
  impl::StoreStrAsAln(aln_str1, aln_str2, s1, s2, i, j, alignments);
}

AlignmentList GlobalAlign(const ConstSequenceHandle& s1, 
                          const ConstSequenceHandle& s2,
                          alg::SubstWeightMatrixPtr& subst,
                          int gap_open,int gap_ext)
{
  impl::AlnMat mat(s1.GetLength()+1, s2.GetLength()+1);
  // init first column of matrix
  if (mat.GetHeight() > 1) {
    mat(0, 1).score = mat(0, 0).score + gap_open;
    mat(0, 1).from = impl::INS1;
  }
  for (int j = 2; j < mat.GetHeight(); ++j) {
    mat(0, j).score = mat(0, j-1).score + gap_ext;
    mat(0, j).from = impl::INS1;
  }
  // fill the alignment matrix
  for (int i = 0; i < mat.GetWidth()-1; ++i) {
    mat(i+1, 0).score = mat(i, 0).score
                      + (mat(i, 0).from==impl::INS2 ? gap_ext : gap_open);
    mat(i+1,0).from = impl::INS2;
    for (int j = 0; j < mat.GetHeight()-1; ++j) {
      char c1=s1[i];
      char c2=s2[j];
      short weight=subst->GetWeight(c1, c2);
      int diag=weight+mat(i, j).score;
      int ins1=mat(i+1, j).score
              +(mat(i+1, j).from==impl::INS1 ? gap_ext : gap_open);
      int ins2=mat(i, j+1).score
              +(mat(i, j+1).from==impl::INS2 ? gap_ext : gap_open);
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
  // write traceback matrix in debug mode
#if !defined(NDEBUG)
  DbgWriteAlnMatrix(mat,s1, s2);
#endif
  AlignmentList alignments;
  Traceback(mat, mat.GetWidth(), mat.GetHeight(), s1, s2, alignments);
  LOG_DEBUG(alignments.back().ToString(80));

  return alignments;
}

}}}
