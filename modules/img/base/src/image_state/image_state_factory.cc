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
  image state factory

  Author: Ansgar Philippsen
*/

#include "image_state_factory.hh"

// single inclusion of all template instantation from here!
#include "image_state_inst.hh"
#include "image_state_base.hh"

#include <ost/img/extent.hh>
#include <ost/message.hh>

namespace ost { namespace img { namespace image_state {

ImageStateBasePtr CreateState(const Extent &e, DataType type, DataDomain dom)
{
  PixelSampling ps(Vec3(1.0,1.0,1.0),dom,e);
  if(type==REAL) {
    if(dom==SPATIAL) {
      return ImageStateBasePtr(new RealSpatialImageState(e,ps));
    } else if (dom==FREQUENCY) {
      return ImageStateBasePtr(new RealFrequencyImageState(e,ps));
    } else {
      throw Error("half-frequency encoding not supported for type REAL");
    }
  } else if (type==COMPLEX) {
    if(dom==HALF_FREQUENCY) {
      return ImageStateBasePtr(new ComplexHalfFrequencyImageState(e,ps));
    } else if(dom==SPATIAL) {
      return ImageStateBasePtr(new ComplexSpatialImageState(e,ps));
    } else {
      return ImageStateBasePtr(new ComplexFrequencyImageState(e,ps));
    }
  } else if (type==WORD) {
    if(dom==SPATIAL) {
      return ImageStateBasePtr(new WordSpatialImageState(e,ps));
    } else {
      throw Error("domain must be SPATIAL for data type WORD");
    }
  }
  throw Error("unexpected end of create state");
  return ImageStateBasePtr(); // to make some compilers happy
}
  
}}} // namespaces
