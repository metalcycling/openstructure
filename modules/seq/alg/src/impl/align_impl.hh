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
#ifndef OST_ALIGN_IMPL_HH
#define OST_ALIGN_IMPL_HH

#include <ost/seq/alignment_handle.hh>
#include <ost/seq/alg/subst_weight_matrix.hh>
//#include "module_config.hh"

namespace ost { namespace seq { namespace alg { namespace impl {

typedef enum {
  DIAG,
  INS1,
  INS2,
  UNKN
} Path;

struct DLLEXPORT AlnPos {
  AlnPos(): score(0), from(impl::UNKN) { }
  int score;
  impl::Path  from;
};

struct DLLEXPORT AlnMat {
  AlnMat(int width, int height): 
    mat_(width*height), width_(width), height_(height)
  { }
  
  impl::AlnPos& operator()(int x, int y) { return mat_[x*height_+y]; }
  
  const impl::AlnPos& operator()(int x, int y) const {
    return mat_[x*height_+y];
  } 
  
  int GetWidth() const { return width_; }
  int GetHeight() const { return height_; }
  
 private:
  std::vector<impl::AlnPos>  mat_;
  int               width_;
  int               height_;  
};
        
inline void DLLEXPORT SetRoute(impl::AlnMat& mat, int& i, int& j,
                     const ConstSequenceHandle& s1, 
                     const ConstSequenceHandle& s2,
                     String& aln_str1, String& aln_str2)
{
  switch (mat(i, j).from) {
    case impl::DIAG:
      --i;
      --j;
      aln_str1.push_back(s1[i]);
      aln_str2.push_back(s2[j]);
      break;
    case impl::INS1:
      --j;
      aln_str1.push_back('-');
      aln_str2.push_back(s2[j]);
      break;
    case impl::INS2:
      --i;
      aln_str1.push_back(s1[i]);
      aln_str2.push_back('-');
      break;
    default:
      assert(0 && "should never get here");
  }
};

inline void DLLEXPORT StoreStrAsAln(String& aln_str1, String& aln_str2,
                          const ConstSequenceHandle& s1, 
                          const ConstSequenceHandle& s2,
                          const int i, const int j,
                          AlignmentList& alignments)
{
  for (size_t x=0; x<aln_str1.size()/2; ++x) {
    std::swap(aln_str1[x], aln_str1[aln_str1.size()-x-1]);
    std::swap(aln_str2[x], aln_str2[aln_str1.size()-x-1]);    
  }
  AlignmentHandle aln=CreateAlignment();
  aln.AddSequence(CreateSequence(s1.GetName(), aln_str1));
  aln.AddSequence(CreateSequence(s2.GetName(), aln_str2));  

  if (s1.GetAttachedView()) {
    aln.AttachView(0, s1.GetAttachedView());
  }

  if (s2.GetAttachedView()) {
    aln.AttachView(1, s2.GetAttachedView());
  }

  aln.SetSequenceOffset(0, i);
  aln.SetSequenceOffset(1, j);
  alignments.push_back(aln);
};

#if !defined(NDEBUG)
inline void DLLEXPORT DbgWriteAlnMatrix(impl::AlnMat& mat,
                                        const ConstSequenceHandle& s1, 
                                        const ConstSequenceHandle& s2)
{
  std::stringstream ss;
  ss << "   ";
  for (int j = 1; j < mat.GetHeight(); ++j) {
    ss << s2[j-1];
  }
  ss << std::endl;
  for (int i = 0; i < mat.GetWidth(); ++i) {
    for (int j = 0; j < mat.GetHeight(); ++j) {
      if (j==0 && i>0) {
        ss << s1[i-1] << " ";
      }
      if (i==0 && j==0) {
        ss << "  ";
      }
      if (mat(i, j).from==impl::DIAG) {
        ss << "\\";
      } else if (mat(i, j).from==impl::INS1) {
        ss << "-";
      } else if (mat(i, j).from==impl::INS2) {
        ss << "|";
      } else if (mat(i, j).from==impl::UNKN) {
        ss << "$";
      }
    }
    ss << std::endl;
  }
  LOG_DEBUG(ss.str());
};
#endif

}}}}

#endif
