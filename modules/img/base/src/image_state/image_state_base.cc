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
  Author: Ansgar Philippsen
*/

#include "image_state_base.hh"
#include "dispatch.hh"

namespace ost { namespace img { namespace image_state {

// define functors to be dispatched
namespace {
template<typename T1, class D1, typename T2, class D2>
struct fnc_add_ip {
  void operator()(ImageStateImpl<T1,D1>* lhs, const ImageStateImpl<T2,D2>* rhs)
  {
    Extent ov=Overlap(lhs->GetExtent(),rhs->GetExtent());
    for(ExtentIterator it(ov);!it.AtEnd();++it) {
      lhs->Value(it)+=Val2Val<T2,T1>(rhs->Value(it));
    }
  }
};

// sub
template<typename T1, class D1, typename T2, class D2>
struct fnc_sub_ip {
  void operator()(ImageStateImpl<T1,D1>* lhs, const ImageStateImpl<T2,D2>* rhs)
  {
    Extent ov=Overlap(lhs->GetExtent(),rhs->GetExtent());
    for(ExtentIterator it(ov);!it.AtEnd();++it) {
      lhs->Value(it)-=Val2Val<T2,T1>(rhs->Value(it));
    }
  }
};

// mul
template<typename T1, class D1, typename T2, class D2>
struct fnc_mul_ip {
  void operator()(ImageStateImpl<T1,D1>* lhs, const ImageStateImpl<T2,D2>* rhs)
  {
    Extent ov=Overlap(lhs->GetExtent(),rhs->GetExtent());
    for(ExtentIterator it(ov);!it.AtEnd();++it) {
      lhs->Value(it)*=Val2Val<T2,T1>(rhs->Value(it));
    }
  }
};

// div
template<typename T1, class D1, typename T2, class D2>
struct fnc_div_ip {
  void operator()(ImageStateImpl<T1,D1>* lhs, const ImageStateImpl<T2,D2>* rhs)
  {
    Extent ov=Overlap(lhs->GetExtent(),rhs->GetExtent());
    for(ExtentIterator it(ov);!it.AtEnd();++it) {
      lhs->Value(it)/=Val2Val<T2,T1>(rhs->Value(it));
    }
  }
};
} // anon ns

ImageStateBase& ImageStateBase::operator+=(const ImageStateBase& b)
{
  static dispatch::binary_dispatch_ip<fnc_add_ip> f_add_ip;
  f_add_ip(this,&b);
  return *this;
}

ImageStateBase& ImageStateBase::operator-=(const ImageStateBase& b)
{
  static dispatch::binary_dispatch_ip<fnc_sub_ip> f_sub_ip;
  f_sub_ip(this,&b);
  return *this;
}

ImageStateBase& ImageStateBase::operator*=(const ImageStateBase& b)
{
  static dispatch::binary_dispatch_ip<fnc_mul_ip> f_mul_ip;
  f_mul_ip(this,&b);
  return *this;
}

ImageStateBase& ImageStateBase::operator/=(const ImageStateBase& b)
{
  static dispatch::binary_dispatch_ip<fnc_div_ip> f_div_ip;
  f_div_ip(this,&b);
  return *this;
}

}}}
