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
  concrete transformations

  Author: Ansgar Philippsen
*/

#ifndef IMG_TRANSFORMATION_IMPL_HH
#define IMG_TRANSFORMATION_IMPL_HH

#include <ost/img/alg/module_config.hh>
#include "transformation_base.hh"

namespace ost { namespace img { namespace alg {

//! rotate in XY plane
class DLLEXPORT_IMG_ALG Rotate2D : public Transformation {
public:
  //! initialize 2D rotation with angle in radians
  Rotate2D(Real angle);
};

//! rotate 90 degree CW in XY plane 
class DLLEXPORT_IMG_ALG Rotate2D90: public Transformation {
public:
  //! initialize amount of 90deg CW rotations
  /*!
    step parameter: 2 corresponds to a 180deg rotation, 3 
    to a 90deg CCW roation, 0 does nothing, all other values
    are wrapped to lie within the range 0-3 (e.g. -1 equals 3)
  */
  Rotate2D90(int step);
};

//! apply scale
class DLLEXPORT_IMG_ALG Scale: public Transformation {
public:
  //! initialize uniform scaling
  Scale(Real s);
};

//! apply translation
class DLLEXPORT_IMG_ALG Translate: public Transformation {
public:
  //! initialize translation along vector
  Translate(const Vec3& t);
};


}}} // ns

#endif
