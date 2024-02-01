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

#include <ost/message.hh>
#include <ost/log.hh>

#include "discrete_shrink.hh"

namespace ost { namespace img { namespace alg {

namespace {

template <typename T, class D>
void do_shrink(const ImageStateImpl<T,D>& isi, T& dest, const Extent& inner_ext);

template <>
void do_shrink(const ImageStateImpl<Word,SpatialDomain>& isi, Word& dest, const Extent& inner_ext)
{
  Real sum=0.0;
  Real count=0.0;
  for(ExtentIterator bit(inner_ext);!bit.AtEnd();++bit) {
    sum+=static_cast<Real>(isi.Value(bit));
    count+=1.0;
  }

  dest = static_cast<Word>(sum/count);
}

template <typename T, class D>
void do_shrink(const ImageStateImpl<T,D>& isi, T& dest, const Extent& inner_ext)
{
  T sum(0);
  Real count=0.0;
  for(ExtentIterator bit(inner_ext);!bit.AtEnd();++bit) {
    sum+=isi.Value(bit);
    count+=1.0;
  }

  dest = sum/count;
}

}

template <typename T, class D>
ImageStateBasePtr DiscreteShrinkFnc::VisitState(const ImageStateImpl<T,D>& isi) const
{
  Size size(isi.GetLogicalExtent().GetSize());
  Point start(isi.GetLogicalExtent().GetStart());
  Size newsize;

  div_t dt;
  dt = div(static_cast<int>(size[0]),bs_[0]); 
  newsize[0] = dt.quot + ((dt.rem>0) ? 1 : 0);
  dt = div(static_cast<int>(size[1]),bs_[1]);
  newsize[1] = dt.quot + ((dt.rem>0) ? 1 : 0);
  dt = div(static_cast<int>(size[2]),bs_[2]);
  newsize[2] = dt.quot + ((dt.rem>0) ? 1 : 0);
  Extent new_ext(newsize);
  Point newstart;

  dt = div(start[0],bs_[0]);  newstart[0] = dt.quot;
  dt = div(start[1],bs_[1]);  newstart[1] = dt.quot;
  dt = div(start[2],bs_[2]);  newstart[2] = dt.quot;
  
  new_ext.Shift(newstart);

  LOG_DEBUG("extent of shrunken image" << isi.GetExtent() << " " << new_ext<< new_ext.GetEnd());


  geom::Vec3 ao = isi.GetAbsoluteOrigin();
  PixelSampling new_ps = isi.GetSampling();
  new_ps.SetExtent(new_ext);

  typename ImageStateImpl<T,D>::SharedPtrType ni(new ImageStateImpl<T,D>(new_ext,new_ps));
  ni->SetAbsoluteOrigin(ao);
  ni->GetSampling().SetPixelSampling(CompMultiply(ni->GetSampling().GetPixelSampling(),Vec3(bs_[0],bs_[1],bs_[2])));

  for(ExtentIterator it(ni->GetExtent(),ni->GetDomain()); !it.AtEnd(); ++it) {
    Point t(it[0]*bs_[0],
	    it[1]*bs_[1],
	    it[2]*bs_[2]);
    Extent inner_ext = Overlap(Extent(t,bs_),isi.GetExtent());
    do_shrink(isi,ni->Value(it),inner_ext);
  }



  return ni;
}

}

template class TEMPLATE_DEF_EXPORT image_state::ImageStateConstModOPAlgorithm<alg::DiscreteShrinkFnc>;

}} // ns
