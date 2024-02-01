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

#ifndef IMG_ALG_TRANSFORM_MIRROR_HH
#define IMG_ALG_TRANSFORM_MIRROR_HH

#include <ost/img/image_state.hh>
#include <ost/img/value_util.hh>
#include <ost/img/alg/module_config.hh>
namespace ost { namespace img { namespace alg {

struct DLLEXPORT_IMG_ALG MirrorFnc {
  MirrorFnc(): planes_(0) {}
  MirrorFnc(int p): planes_(p) {}


  template <typename T, class D>
  ImageStateBasePtr VisitState(const ImageStateImpl<T,D>& in_state) const {
    boost::shared_ptr<ImageStateImpl<T,D> > out_state(new ImageStateImpl<T,D>(in_state.GetExtent().Mirror(this->planes_),in_state.GetSampling()));
    
    for(ExtentIterator it(in_state.GetExtent()); !it.AtEnd(); ++it) {
      Point p(it);
      out_state->Value(p.Mirror(this->planes_))=in_state.Value(p);
    }
    
    return out_state;
  }
  
  static String GetAlgorithmName() {return "Mirror";}
private:
  int planes_;
};

typedef ImageStateConstModOPAlgorithm<MirrorFnc> Mirror;

}}} // ns

// no template instantiation needed for completely inlined algorithm

#endif
