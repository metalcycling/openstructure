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
  Authors: Andreas Schenk, Ansgar Philippsen
*/

#include <ost/img/image_state.hh>

#include "mask_image.hh"

namespace ost { namespace img { namespace alg {

MaskImageBase::MaskImageBase(){} // To make the compiler happy

MaskImageBase::MaskImageBase(const MaskPtr& m, bool float_bg):
  mask_(m),float_bg_(float_bg)
{}

template <typename T, class D>
void MaskImageBase::VisitState(ImageStateImpl<T,D>& is) const
{
  T val(0);
  if(float_bg_){
    unsigned int count=0;
    for(ExtentIterator it(is.GetExtent());!it.AtEnd();++it) {
      if(mask_->IsInside(Point(it).ToVec2())) {
        val+=is.Value(it);
        ++count;
      }
    }
    if(count>0){
      val/=static_cast<Real>(count);
    }
  }
  for(ExtentIterator it(is.GetExtent());!it.AtEnd();++it) {
    if(mask_->IsOutside(Point(it).ToVec2())){
      is.Value(it)=val;
    }
  }    
}

}

namespace image_state {

template class TEMPLATE_DEF_EXPORT ImageStateConstModIPAlgorithm<alg::MaskImageBase>;
  
}}} // namespaces
