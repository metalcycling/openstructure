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

#include "paste_impl.hh"

#include "image_state/dispatch.hh"
#include "image.hh"

namespace ost { namespace img { namespace detail {

// define functors to be dispatched
namespace {
// paste
template<typename T1, class D1, typename T2, class D2>
struct fnc_paste_ip {
  void operator()(image_state::ImageStateImpl<T1,D1>* lhs,
                  const image_state::ImageStateImpl<T2,D2>* rhs) {
    Extent ov=Overlap(lhs->GetExtent(),rhs->GetExtent());
    for(ExtentIterator it(ov);!it.AtEnd();++it) {
      lhs->Value(it)=Val2Val<T2,T1>(rhs->Value(it));
    }
  }
};

// partial specialization for half_frequency rhs
template<typename T1, class D1>
struct fnc_paste_ip<T1,D1,Complex,image_state::HalfFrequencyDomain> {
  void operator()(image_state::ImageStateImpl<T1,D1>* lhs,
                  const image_state::ImageStateImpl<Complex,
                          image_state::HalfFrequencyDomain>* rhs)
  {
    Extent ov=Overlap(lhs->GetExtent(),rhs->GetLogicalExtent());
    for(ExtentIterator it(ov);!it.AtEnd();++it) {
      lhs->Value(it)=Val2Val<Complex,T1>(rhs->GetComplex(it));
    }
  }
};
} // anon ns

PasteFnc::PasteFnc():
  target_()
{}

PasteFnc::PasteFnc(const ImageHandle& h):
  target_(h)
{}

template<typename V, class D>
void PasteFnc::VisitState(const image_state::ImageStateImpl<V,D>& isi)
{
  static image_state::dispatch::binary_dispatch_ip<fnc_paste_ip> f_paste_ip;
  f_paste_ip(target_.ImageStatePtr().get(),&isi);
}


}

#ifdef NO_EXPL_INST
#else
template class image_state::ImageStateNonModAlgorithm<detail::PasteFnc>;
#endif

}} // ns
