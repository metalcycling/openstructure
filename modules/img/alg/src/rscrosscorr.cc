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
#include <ost/img/image_state.hh>
#include <ost/img/alg/rscrosscorr.hh>

namespace ost { namespace img { namespace alg {

Real RealSpatialCrossCorrelation(const ConstImageHandle& image1,
                                 const ConstImageHandle& image2,
                                 const Extent& extent)
{
  img::Extent image1_overlap=Overlap(image1.GetExtent(),extent);
  img::Extent image2_overlap=Overlap(image2.GetExtent(),extent);

  if ( (image1.IsReal() && image2.IsReal())!=true ) {

    throw(Error("One of the input images is not of Real type"));

  } else if ( (image1.IsSpatial() && image2.IsSpatial()) !=true) {

    throw(Error("One of the input images is not in Spatial domain"));

  } else if (image1_overlap != image2_overlap) {

    throw(Error("Input images and extent are not compatible with each other"));

  } else {

    ImageStateBasePtr base_is_pointer_image1 = image1.ImageStatePtr();
    boost::shared_ptr<img::RealSpatialImageState> image1_pointer = boost::dynamic_pointer_cast<img::RealSpatialImageState>(base_is_pointer_image1);

    ImageStateBasePtr base_is_pointer_image2 = image2.ImageStatePtr();
    boost::shared_ptr<img::RealSpatialImageState> image2_pointer = boost::dynamic_pointer_cast<img::RealSpatialImageState>(base_is_pointer_image2);

    Real corr=0.0;
    Real sum1=0.0;
    Real sum2=0.0;

    // total of values
    for(ExtentIterator img_ov_it(image1_overlap);!img_ov_it.AtEnd();++img_ov_it) {
      Real v1=Val2Val<Real,Real>(image1_pointer->Value(img_ov_it));
      Real v2=Val2Val<Real,Real>(image2_pointer->Value(img_ov_it));
      corr += v1 * v2;
      sum1 += v1 * v1;
      sum2 += v2 * v2;  
    }
    
    return corr/(sqrt(sum1) * sqrt (sum2));

  }
}

Real RealSpatialCrossCorrelation(const ConstImageHandle& image1,
                                 const ConstImageHandle& image2)
{
  return RealSpatialCrossCorrelation(image1, image2, image2.GetExtent());
}

}}} // ns

