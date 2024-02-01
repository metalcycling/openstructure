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

#include "correlate.hh"

#include <ost/img/image_state/dispatch.hh>

namespace ost { namespace img { namespace alg {

namespace {

template <typename X>
X Conj(const X& x);

template<>
Complex Conj(const Complex& x) {return std::conj(x);}

template <typename X>
X Conj(const X& x) {return x;}

// explicit correlation with wrap around
template<typename T1, class D1, typename T2, class D2>
struct fnc_expl_correlate_wrap_op {
  ImageStateBasePtr operator()(const ImageStateImpl<T1,D1>* lhs,
			       const ImageStateImpl<T2,D2>* rhs) {
    boost::shared_ptr<ImageStateImpl<T1,D1> > res = lhs->CloneState(false);

    Extent ex1=lhs->GetExtent();
    for(ExtentIterator it1(ex1);!it1.AtEnd(); ++it1) {
      T1 val(0); // default ctor sets value to zero
      for(ExtentIterator it2(rhs->GetExtent());!it2.AtEnd();++it2) {
	val+=lhs->Value(ex1.WrapAround(Point(it1)+Point(it2))) * Val2Val<T2,T1>(Conj(rhs->Value(it2)));
      }
      res->Value(it1) = val;
    }
    return res;
  }
};

// explicit correlation without wrap around
template<typename T1, class D1, typename T2, class D2>
struct fnc_expl_correlate_op {
  ImageStateBasePtr operator()(const ImageStateImpl<T1,D1>* lhs,
			       const ImageStateImpl<T2,D2>* rhs) {
    throw("non wrap-around correlate not yet implemented");
  }
};

} // anon ns

ExplicitCorrelateFnc::ExplicitCorrelateFnc():
  ref_()
{}

ExplicitCorrelateFnc::ExplicitCorrelateFnc(const ConstImageHandle& ref, bool wrap):
  ref_(ref), wrap_(wrap)
{}

template <typename V, class D>
ImageStateBasePtr ExplicitCorrelateFnc::VisitState(const ImageStateImpl<V,D>& isi) const 
{
  if(wrap_) {
    image_state::dispatch::binary_dispatch_op<fnc_expl_correlate_wrap_op> cnv;
    return cnv(&isi,ref_.ImageStatePtr().get());
  } else {
    image_state::dispatch::binary_dispatch_op<fnc_expl_correlate_op> cnv;
    return cnv(&isi,ref_.ImageStatePtr().get());
  }
}

}
template class TEMPLATE_DEF_EXPORT image_state::ImageStateConstModOPAlgorithm<alg::ExplicitCorrelateFnc>;
}} //  ns

