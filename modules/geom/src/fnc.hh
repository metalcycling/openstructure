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
#ifndef GEOM_FNC_HH
#define GEOM_FNC_HH

#include "vec.hh"
#include "composite.hh"

namespace geom {

//! provides line line interaction results
struct LineLineInteraction {
  LineLineIteraction(const Line& l1, const Line& l2);
  Vec3 point1;
  Real coeff1;
  Vec3 point2;
  Real coeff2;
  Real distance;
};

//! returns the distance between two points
Real Distance(const Vec3& p1, const Vec3& p2);
//! returns the distance between a point and a line
Real Distance(const Vec3& p, const Line& l);
//! returns the distance between a point and a line
Real Distance(const Line& l, const Vec3& p);
//! returns the minimal distance between the points in two Vec3List
Real MinDistance(const Vec3List& l1, const Vec3List& l2);
} // ns geom

#endif
