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
#include <iostream>

#include <ost/img/value_util.hh>

#include "local_sigma_threshold.hh"

namespace ost { namespace img { namespace alg {

LocalSigmaThresholdBase::LocalSigmaThresholdBase():
  size_(5),
  level_(1.0)
{}

LocalSigmaThresholdBase::LocalSigmaThresholdBase(int size, Real level):
  size_(size),
  level_(level)
{}

template <typename T, class D>
ImageStateBasePtr LocalSigmaThresholdBase::VisitState(const ImageStateImpl<T,D>& s) const
{
  Extent full_extent=s.GetExtent();
  Extent search_extent = Extent(full_extent.GetStart()+Point(size_,size_),
				full_extent.GetEnd()-Point(size_,size_));
  boost::shared_ptr<RealSpatialImageState> out_state(new RealSpatialImageState(search_extent,s.GetSampling()));

  Real level2=level_*level_;
  Real fac = 1.0/static_cast<Real>((size_*2+1)*(size_*2+1));
  for(ExtentIterator it(search_extent);!it.AtEnd();++it) {
    Point p(it);
    Real sum=0.0;
    Real sum2=0.0;
    for(int a=-size_;a<=size_;++a) {
      for(int b=-size_;b<=size_;++b) {
	Real val=Val2Val<T,Real>(s.Value(p+Point(a,b)));
	sum+=val;
	sum2+=val*val;
      }
    }
    Real mean=sum*fac;
    Real var=sum2*fac-mean*mean;
    out_state->Value(p) = var>level2 ? 1.0 : 0.0;
  }

  return out_state;
}


}

namespace image_state {
  
template class TEMPLATE_DEF_EXPORT ImageStateConstModOPAlgorithm<alg::LocalSigmaThresholdBase>;  

}}} // ns

