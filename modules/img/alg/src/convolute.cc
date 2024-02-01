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

#include "convolute.hh"

#include <ost/img/image_state/dispatch.hh>

namespace ost { namespace img { namespace alg {

namespace {

// explicit convolution with wrap around
template<typename T1, class D1, typename T2, class D2>
struct fnc_expl_convolute_wrap_op {
  ImageStateBasePtr operator()(const ImageStateImpl<T1,D1>* lhs,
			       const ImageStateImpl<T2,D2>* rhs) {
    boost::shared_ptr<ImageStateImpl<T1,D1> > res = lhs->CloneState(false);

    Extent ex1=lhs->GetExtent();
    for(ExtentIterator it1(ex1);!it1.AtEnd(); ++it1) {
      T1 val(0);
      for(ExtentIterator it2(rhs->GetExtent());!it2.AtEnd();++it2) {
	val+=lhs->Value(ex1.WrapAround(Point(it1)-Point(it2))) * Val2Val<T2,T1>(rhs->Value(it2));
      }
      res->Value(it1) = val;
    }
    return res;
  }
};

// explicit convolution without wrap around
template<typename T1, class D1, typename T2, class D2>
struct fnc_expl_convolute_op {
  ImageStateBasePtr operator()(const ImageStateImpl<T1,D1>* lhs,
			       const ImageStateImpl<T2,D2>* rhs) {
    boost::shared_ptr<ImageStateImpl<T1,D1> > res = lhs->CloneState(false);

    Extent ex1=lhs->GetExtent();
    for(ExtentIterator it1(ex1);!it1.AtEnd(); ++it1) {
      T1 val(0);
      for(ExtentIterator it2(rhs->GetExtent());!it2.AtEnd();++it2) {
	val+=lhs->GetCheckedValue(Point(it1)-Point(it2)) * Val2Val<T2,T1>(rhs->Value(it2));
      }
      res->Value(it1) = val;
    }
    return res;
  }
};

} // anon ns

ExplicitConvoluteFnc::ExplicitConvoluteFnc():
  ref_()
{}

ExplicitConvoluteFnc::ExplicitConvoluteFnc(const ConstImageHandle& ref, bool wrap):
  ref_(ref.ImageStatePtr()), wrap_(wrap)
{}

ExplicitConvoluteFnc::ExplicitConvoluteFnc(const ImageStateBasePtr& ref, bool wrap):
  ref_(ref), wrap_(wrap)
{}

template <typename V, class D>
ImageStateBasePtr ExplicitConvoluteFnc::VisitState(const ImageStateImpl<V,D>& isi) const 
{
  if(wrap_) {
    image_state::dispatch::binary_dispatch_op<fnc_expl_convolute_wrap_op> cnv;
    return cnv(&isi,ref_.get());
  } else {
    image_state::dispatch::binary_dispatch_op<fnc_expl_convolute_op> cnv;
    return cnv(&isi,ref_.get());
  }
}

}

namespace image_state {
  
template class TEMPLATE_DEF_EXPORT ImageStateConstModOPAlgorithm<alg::ExplicitConvoluteFnc>;  

}}} // ns

