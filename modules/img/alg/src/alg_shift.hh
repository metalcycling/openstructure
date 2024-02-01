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

#ifndef IMG_ALG_TRANSFORM_SHIFT_HH
#define IMG_ALG_TRANSFORM_SHIFT_HH

#include <ost/img/image_state.hh>
#include <ost/img/value_util.hh>
#include <ost/img/alg/module_config.hh>

namespace {

unsigned int absmod(int x, unsigned int y)
{
  return x<0 ? y+std::div(x,y).rem : std::div(x,y).rem;
}

} // anon ns

namespace ost { namespace img { namespace alg {

class DLLEXPORT_IMG_ALG ShiftFnc {
public:
  ShiftFnc(): shift_() {}
  ShiftFnc(const Point& s): shift_(s) {}

  template <typename T, class D>
  ImageStateBasePtr VisitState(const ImageStateImpl<T,D>& in_state) const {
    boost::shared_ptr<ImageStateImpl<T,D> > isi = in_state.CloneState(false);

    unsigned int depth=in_state.GetExtent().GetDepth();
    unsigned int height=in_state.GetExtent().GetHeight();
    unsigned int width=in_state.GetExtent().GetWidth();
    
    unsigned int p0=::absmod(shift_[0],width);
    unsigned int p1=::absmod(shift_[1],height);
    unsigned int p2=::absmod(shift_[2],depth);

    for(unsigned int u=0;u<width;++u) {
      for(unsigned int v=0;v<height;++v) {
	for(unsigned int w=0;w<depth;++w) {
	  isi->Value(Index((u+p0)%width,(v+p1)%height,(w+p2)%depth))=in_state.Value(Index(u,v,w));
	}
      }
    }

    return isi;
  }

  static String GetAlgorithmName() {return "Shift";}

private:
  Point shift_;
};

typedef ImageStateConstModOPAlgorithm<ShiftFnc> Shift;

}}} // ns

#endif
