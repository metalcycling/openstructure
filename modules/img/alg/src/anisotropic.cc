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
  authors: Giani Signorell, Ansgar Philippsen
*/

#include <cmath>

#include "anisotropic.hh"

namespace ost { namespace img { namespace alg {

AnisotropicFilter::AnisotropicFilter(Real sigma):
  GaussianFilter(sigma)
{
  b1_b0_ = b1_/b0_;
  b2_b0_ = b2_/b0_;
  b3_b0_ = b3_/b0_;
  n_=0;
}

template <typename T, class D>
void AnisotropicFilter::VisitState(ImageStateImpl<T,D>& s)
{
  Extent ext = s->GetExtent();
  // maximum extent plus border of 3 on each side
  int n = std::max(ext.GetWidth(),std::max(ext.GetHeight(),ext.GetDepth()))+6;

  std::vector<T> tmp_avg(n);
  std::vector<T> tmp_min(n);
  std::vector<T> tmp_max(n);
  
  for(LineIterator li(ext);!li.AtEnd();++li) {
    ExtentIterator bi=li.GetLine();
    lineTransformFW(s, bi, tmp_avg, tmp_min, tmp_max );
    lineTransformBW(s, bi, tmp_avg, tmp_min ,tmp_max );
  }
}

// specialization
template <>
bool AnisotropicFilter::is_larger<Complex>(const Complex& v1, const Complex& v2)
{
  return std::abs(v1)>std::abs(v2);
}

// generic case
template <typename T>
bool AnisotropicFilter::is_larger(const T& v1, const T& v2)
{
  return v1>v2;
}


template <typename T, class D>
void AnisotropicFilter::lineTransformFW(ImageStateImpl<T,D>& s,
					ExtentIterator& it, 
					std::vector<T>& tmp_avg, 
					std::vector<T>& tmp_min, 
					std::vector<T>& tmp_max)
{
  tmp_avg[0] = s->Value(it);
  tmp_avg[1] = tmp_avg[0];
  tmp_avg[2] = tmp_avg[0];
  tmp_min[0] = tmp_avg[0];
  tmp_min[1] = tmp_avg[0];
  tmp_min[2] = tmp_avg[0];
  tmp_max[0] = tmp_avg[0];
  tmp_max[1] = tmp_avg[0];
  tmp_max[2] = tmp_avg[0];
  n_=3; 
  
  for(;!it.AtEnd(); ++it){
    T xn = s->Value(it);
    
    tmp_avg[n_] = xn + b1_b0_ * (tmp_avg[n_-1] - xn) +
      b2_b0_ * (tmp_avg[n_-2] - xn) + 
      b3_b0_ * (tmp_avg[n_-3] - xn);
    
    T xmin = xn;
    if( is_larger(xn,tmp_min[n_-1])) {
      xmin += b1_b0_ * (tmp_min[n_-1] - xn);
    }
    if( is_larger(xn,tmp_min[n_-2])) {
      xmin += b2_b0_ * (tmp_min[n_-2] - xn);
    }
    if( is_larger(xn,tmp_min[n_-3])) {
      xmin += b3_b0_ * (tmp_min[n_-3] - xn);
    }
    tmp_min[n_] = xmin;
    
    T xmax = xn;
    if( is_larger(tmp_max[n_-1],xn)) {
      xmax += b1_b0_ * (tmp_max[n_-1] - xn);
    }
    if( is_larger(tmp_max[n_-2],xn)) {
      xmax += b2_b0_ * (tmp_max[n_-2] - xn);
    }
    if( is_larger(tmp_max[n_-3],xn)) {
      xmax += b3_b0_ * (tmp_max[n_-3] - xn);
    }
    tmp_max[n_] = xmax;
    
    ++n_;
  }
}

template <typename T, class D>
void AnisotropicFilter::lineTransformBW(ImageStateImpl<T,D>& s,
                                        ExtentIterator& it,
                                        std::vector<T>& tmp_avg, 
                                        std::vector<T>& tmp_min, 
std::vector<T>& tmp_max)
{
  --it.ToEnd();
  --n_;
  tmp_avg[n_+1] = s->Value(it);
  tmp_avg[n_+2] = tmp_avg[n_+1];
  tmp_avg[n_+3] = tmp_avg[n_+1];
  tmp_min[n_+1] = tmp_avg[n_+1];
  tmp_min[n_+2] = tmp_avg[n_+1];
  tmp_min[n_+3] = tmp_avg[n_+1];
  tmp_max[n_+1] = tmp_avg[n_+1];
  tmp_max[n_+2] = tmp_avg[n_+1];
  tmp_max[n_+3] = tmp_avg[n_+1];
  
  for(;!it.AtStart(); --it){
    tmp_avg[n_] = tmp_avg[n_] +
      b1_b0_*(tmp_avg[n_+1] - tmp_avg[n_]) + 
      b2_b0_*(tmp_avg[n_+2] - tmp_avg[n_]) +
      b3_b0_*(tmp_avg[n_+3] - tmp_avg[n_]);
    
    T xmin = tmp_min[n_];
    if( is_larger(tmp_min[n_],tmp_min[n_+1])) {
      xmin += b1_b0_ * (tmp_min[n_+1] - tmp_min[n_]);
    }
    if( is_larger(tmp_min[n_],tmp_min[n_+2])) {
      xmin += b2_b0_ * (tmp_min[n_+2] - tmp_min[n_]);
    }
    if( is_larger(tmp_min[n_],tmp_min[n_+3])) {
      xmin += b3_b0_ * (tmp_min[n_+3] - tmp_min[n_]);
    }
    
    tmp_min[n_] = xmin;
    
    T xmax = tmp_max[n_];
    
    if( is_larger(tmp_max[n_], tmp_max[n_+1] )) {
      xmax += b1_b0_ * (tmp_max[n_+1] - tmp_max[n_]);
    }
    if( is_larger(tmp_max[n_], tmp_max[n_+2] )) {
      xmax += b2_b0_ * (tmp_max[n_+2] - tmp_max[n_]);
    }
    if( is_larger(tmp_max[n_], tmp_max[n_+3] )) {
      xmax += b3_b0_ * (tmp_max[n_+3] - tmp_max[n_]);
    }
    tmp_max[n_] = xmax;
    
    if( is_larger(tmp_avg[n_],s->Value(it))) {
      s->Value(it)=xmin;
    } else {
      s->Value(it)=xmax;
    }
    n_--;
  }
}

}}}
