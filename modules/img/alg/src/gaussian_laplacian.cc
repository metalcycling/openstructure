//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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
  Authors: Giani Signorell, Ansgar Philippsen, Valerio Mariani
*/

#include <cmath>
#include <ost/message.hh>
#include "gaussian_laplacian.hh"

namespace ost {
namespace img {
  namespace alg {

GaussianLaplacianFilterBase::GaussianLaplacianFilterBase(Real sigma)
{
  calcBVals( calcQ(sigma) );
}

template <typename T, class D>
void GaussianLaplacianFilterBase::VisitState(ImageStateImpl<T,D>& s)
{

  ImageStateImpl<T,D> a = s;
  ImageStateImpl<T,D> b = s;
  ImageStateImpl<T,D> c = s;

  for(LineIterator li(a.GetExtent(),Axis::X);!li.AtEnd();++li) {
    ExtentIterator bi = li.GetLine();
    lineTransformFW(a,bi);
    lineTransformBW(a,bi);
  }

  for(LineIterator li(b.GetExtent(),Axis::Y);!li.AtEnd();++li) {
    ExtentIterator bi = li.GetLine();
    lineTransformFW(b,bi);
    lineTransformBW(b,bi);
  }

  for(LineIterator li(c.GetExtent(),Axis::Z);!li.AtEnd();++li) {
    ExtentIterator bi = li.GetLine();
    lineTransformFW(c,bi);
    lineTransformBW(c,bi);
  }

  for(ExtentIterator ex(s.GetExtent());!ex.AtEnd();++ex) {
     Real val1 = Val2Val<T,Real>(a.Value(ex));
     Real val2 = Val2Val<T,Real>(b.Value(ex));
     Real val3 = Val2Val<T,Real>(c.Value(ex));
     if (s.GetExtent().GetHeight()==1 ) val2=0;
     if (s.GetExtent().GetDepth()==1 ) val3=0;
     s.Value(ex)=val1 + val2 + val3; 
  }

}

template <typename T, class D>
void GaussianLaplacianFilterBase::lineTransformFW(ImageStateImpl<T,D>& s, ExtentIterator& it) 
{
  T value_prev_it=s.Value(it);
  T wn_1 = s.Value(it);
  T wn_2 = wn_1;
  T wn_3 = wn_1;
  T wn;
  for(;!it.AtEnd(); ++it){
    wn = bb_ * (s.Value(it)-value_prev_it) + (b1_*wn_1 + b2_*wn_2 + b3_*wn_3)/b0_;
    value_prev_it=s.Value(it);
    s.Value(it)=wn;
    wn_3 = wn_2;
    wn_2 = wn_1;
    wn_1 = wn;
  }
}

template <typename T, class D>
void GaussianLaplacianFilterBase::lineTransformBW(ImageStateImpl<T,D>& s, ExtentIterator& it)
{
  --it.ToEnd();
  T value_next_it=s.Value(it);
  T out_1 = s.Value(it);
  T out_2 = out_1;
  T out_3 = out_1;
  T out_n(0);
  
  for(;!it.AtStart(); --it){
    out_n = bb_ *(value_next_it-s.Value(it)) + (b1_*out_1 + b2_*out_2 + b3_*out_3)/b0_;
    value_next_it=s.Value(it);
    s.Value(it)=out_n;
    out_3 = out_2;
    out_2 = out_1;
    out_1 = out_n;
  }
}

// protected and private methods

Real GaussianLaplacianFilterBase::calcQ( Real sigma )
{
  if(sigma >= 2.5) {
    return 0.98711 * sigma -0.9633;
  } else if(sigma >= 0.5) {
    return 3.97156 - 4.14554 * std::sqrt( 1 - 0.26891 * sigma);
  } else {
    //not defined in the paper
    throw(Error("gaussian filter: sigma must be larger than 0.5"));
  }
}

void GaussianLaplacianFilterBase::calcBVals( Real q )
{
  Real q2=q*q;
  Real q3=q2*q;
  b0_ = 1.57825 + (2.44413 * q) + (1.4281 * q2) + (0.422205 * q3);
  b1_ = (2.44413 * q) + (2.85619 * q2) + (1.26661 * q3);
  b2_ = - ( (1.4281 * q2) + (1.26661 * q3) );
  b3_ = 0.422205 * q3;
  bb_ = 1 - ((b1_ + b2_ + b3_) / b0_);
}

} // ns

template class TEMPLATE_DEF_EXPORT image_state::ImageStateModIPAlgorithm<alg::GaussianLaplacianFilterBase>;
}}
