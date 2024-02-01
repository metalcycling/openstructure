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
  Authors: Valerio Mariani
*/

#include <ost/img/image_state.hh>
#include <ost/img/alg/fill.hh>
#include <ost/img/alg/gaussian.hh>

#include "smooth_mask_image.hh"

namespace ost { namespace img { namespace alg {

SmoothMaskImageBase::SmoothMaskImageBase(){}; // To make the compiler happy


SmoothMaskImageBase::SmoothMaskImageBase(const MaskPtr& m, Real sigma, bool float_bg):
  mask_(m),sigma_(sigma),float_bg_(float_bg)
{}

template <typename T, class D>
void SmoothMaskImageBase::VisitState(ImageStateImpl<T,D>& is)
{
  T val(0);
  RealSpatialImageState mask_image_is(is.GetExtent(),PixelSampling());
  if(float_bg_){
    unsigned int count=0;
    for(ExtentIterator it(is.GetExtent());!it.AtEnd();++it) {
      if(mask_->IsInside(Point(it).ToVec2())) {
        val+=is.Value(it);
        ++count;
        mask_image_is.Value(it)=1.0;
      }else{
        mask_image_is.Value(it)=0.0;
      }
    }
    if(count>0){
      val/=static_cast<Real>(count);
    }
    GaussianFilter gaussian_filter_alg(sigma_);
    mask_image_is.ApplyIP(gaussian_filter_alg);
    is-=val;
    // weird upcast hack to make the compiler find the *= op
    static_cast<image_state::ImageStateBase&>(is).operator*=(mask_image_is);
    is+=val;
  }else{
    for(ExtentIterator it(is.GetExtent());!it.AtEnd();++it) {
      if(mask_->IsInside(Point(it).ToVec2())) {
        mask_image_is.Value(it)=1.0;
      }else{
        mask_image_is.Value(it)=0.0;
      }
    }
    GaussianFilter gaussian_filter_alg(sigma_);
    mask_image_is.ApplyIP(gaussian_filter_alg);
    // weird upcast hack to make the compiler find the *= op
    static_cast<image_state::ImageStateBase&>(is).operator*=(mask_image_is);
  }
    
}


}

namespace image_state {

template class TEMPLATE_DEF_EXPORT ImageStateModIPAlgorithm<alg::SmoothMaskImageBase>;


}}} // namespaces
