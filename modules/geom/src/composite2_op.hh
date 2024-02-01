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
#ifndef GEOM_COMPOSITE2_OP_HH
#define GEOM_COMPOSITE2_OP_HH

#include "constants.hh"
#include "vec2.hh"
#include "composite2.hh"

namespace geom {


Real DLLEXPORT_OST_GEOM Angle(const Line2& l1, const Line2& l2);
Vec2 DLLEXPORT_OST_GEOM IntersectionPoint(const Line2& l1, const Line2& l2);

Real DLLEXPORT_OST_GEOM Distance(const Line2& l, const Vec2& v);

bool DLLEXPORT_OST_GEOM IsOnLine(const Line2& l, const Vec2& v, Real ephilon=EPSILON);
bool DLLEXPORT_OST_GEOM operator==(const Line2& l1, const Line2& l2);
bool DLLEXPORT_OST_GEOM Equal(const Line2& l1, const Line2& l2,Real ephilon=EPSILON);
bool DLLEXPORT_OST_GEOM EqualPosition(const Line2& l1, const Line2& l2,Real ephilon=EPSILON);

bool DLLEXPORT_OST_GEOM AreParallel(const Line2& l1, const Line2& l2,Real ephilon=EPSILON);

DLLEXPORT_OST_GEOM bool IsInPolygon(const Polygon2& p, const Vec2& v);
DLLEXPORT_OST_GEOM Polygon2 SimplifyPolygon(const Polygon2& p);
DLLEXPORT_OST_GEOM bool IsInRectangle(const Rectangle2& r, const Vec2& v);
DLLEXPORT_OST_GEOM bool IsInCircle(const Circle2& c, const Vec2& v);

} // ns


#endif
