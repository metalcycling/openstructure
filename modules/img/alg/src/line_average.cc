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
#include "line_average.hh"
#include <ost/message.hh>
#include <ost/img/alg/line_iterator.hh>

namespace ost {namespace img { namespace alg {

template <typename T, class D>
ImageStateBasePtr LineAverageBase::VisitState(const ImageStateImpl<T,D>& isi) const
{
  Extent extent = isi.GetExtent();
  Size size = isi.GetSize();
  PixelSampling pixel_sampling = isi.GetSampling();

  if (isi.GetSize()[2] != 1)  {
      throw(Error("This algorithm does not support 3D images"));
  }
  if (geometry_ == Axis::X) {
    boost::shared_ptr<ImageStateImpl<T,D> > mean(new ImageStateImpl<T,D>(Extent(Point(extent.GetStart()),Size(1,size[1])),pixel_sampling));

    int curr_line =0;
    for (LineIterator li(isi.GetExtent(),Axis::X);!li.AtEnd();++li) {
      T sum(0);
      for(ExtentIterator ei(li.GetLine()); !ei.AtEnd(); ++ei)  {
        sum += isi.Value(ei);
        curr_line=Point(ei)[1];
      }
      sum /= static_cast<Real>(li.GetSize());
      mean->Value(Point(0,curr_line)) = sum;
    }
  return mean;
  } else if (geometry_ == Axis::Y) {
    boost::shared_ptr<ImageStateImpl<T,D> > mean(new ImageStateImpl<T,D>(Extent(Point(extent.GetStart()),Size(size[0],1)),pixel_sampling));

    int curr_col =0;
    for (LineIterator li(isi.GetExtent(),Axis::Y);!li.AtEnd();++li) {
      T sum(0);
      for(ExtentIterator ei(li.GetLine()); !ei.AtEnd(); ++ei)  {
        sum += isi.Value(ei);
        curr_col=Point(ei)[0];
      }
      sum /= static_cast<Real>(li.GetSize());
      mean->Value(Point(curr_col,0)) = sum;
    }
  return mean;
  } else {
    throw(Error("This algorithm does not support Z geometry"));
  }
  return ImageStateBasePtr(); // to make the compiler happy
}

} // ns

#ifdef NO_EXPL_INST
namespace {
  static alg::LinColMean x;
}
#else
template class TEMPLATE_DEF_EXPORT image_state::ImageStateConstModOPAlgorithm<alg::LineAverageBase>;
#endif

}} // ns
