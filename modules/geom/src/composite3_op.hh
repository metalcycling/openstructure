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
#ifndef GEOM_COMPOSITE3_OP_HH
#define GEOM_COMPOSITE3_OP_HH

#include "constants.hh"
#include "vec3.hh"
#include "composite3.hh"

namespace geom {

Real DLLEXPORT_OST_GEOM Angle(const Line3& l1, const Line3& l2);
Real DLLEXPORT_OST_GEOM Angle(const Plane& p1, const Plane& p2);
Real DLLEXPORT_OST_GEOM Angle(const Plane& p, const Line3& l);

Vec3 DLLEXPORT_OST_GEOM IntersectionPoint(const Line3& l1, const Line3& l2);
Vec3 DLLEXPORT_OST_GEOM IntersectionPoint(const Line3& l, const Plane& p);
Vec3 DLLEXPORT_OST_GEOM IntersectionPoint(const Plane& p, const Line3& l);
Line3 DLLEXPORT_OST_GEOM IntersectionLine(const Plane& p1, const Plane& p2);

Real DLLEXPORT_OST_GEOM Distance(const Line3& l1, const Line3& l2);
Real DLLEXPORT_OST_GEOM Distance(const Line3& l, const Vec3& v);
Real DLLEXPORT_OST_GEOM Distance(const Plane& p, const Vec3& v);

bool DLLEXPORT_OST_GEOM IsOnLine(const Line3& l, const Vec3& v, Real ephilon=EPSILON);
bool DLLEXPORT_OST_GEOM IsInPlane(const Plane& p,  const Vec3& v,Real ephilon=EPSILON);
bool DLLEXPORT_OST_GEOM IsInPlane(const Plane& p,  const Line3& l,Real ephilon=EPSILON);
int  DLLEXPORT_OST_GEOM RelativePosition(const Plane& p,  const Vec3& v,Real ephilon=EPSILON);

bool DLLEXPORT_OST_GEOM operator==(const Line3& l1, const Line3& l2);
bool DLLEXPORT_OST_GEOM Equal(const Line3& l1, const Line3& l2,Real ephilon=EPSILON);
bool DLLEXPORT_OST_GEOM EqualPosition(const Line3& l1, const Line3& l2,Real ephilon=EPSILON);
bool DLLEXPORT_OST_GEOM operator==(const Plane& p1, const Plane& p2);
bool DLLEXPORT_OST_GEOM Equal(const Plane& p1, const Plane& p2,Real ephilon=EPSILON);
bool DLLEXPORT_OST_GEOM EqualPosition(const Plane& p1, const Plane& p2,Real ephilon=EPSILON);

bool DLLEXPORT_OST_GEOM AreParallel(const Plane& p,  const Line3& l,Real ephilon=EPSILON);
bool DLLEXPORT_OST_GEOM AreParallel(const Plane& p1, const Plane& p2,Real ephilon=EPSILON);
bool DLLEXPORT_OST_GEOM AreParallel(const Line3& l1, const Line3& l2,Real ephilon=EPSILON);

bool DLLEXPORT_OST_GEOM AreIntersecting(const Line3& l1, const Line3& l2, Real ephilon=EPSILON);

DLLEXPORT bool IsInSphere(const Sphere& p,const Vec3& v);

} // ns


#endif
