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
#ifndef GEOM_QUAT_HH
#define GEOM_QUAT_HH

#include <iostream>
#include <boost/operators.hpp>

#include <ost/geom/module_config.hh>

/*
  Unit Quaternion

  math and code assembled from various places, among them:

  E.B. Dam, M. Koch, M. Lillholm, Technical Report DIKU-TR-98/5,
  Department of Computer Science, University of Copenhagen

  boost quaternion example, (C) Copyright Hubert Holin 2001.
  Distributed under the Boost Software License, Version 1.0.

  Python Computer Graphics Kit, original code by Matthias Baas (C) 2004
  Distributed under the GPL v2.0
*/

namespace geom { 

class Mat3;
class Vec3;

//! unit quaternion
/*!
  The four components of a unit quaternion can be seen as
  [cos(theta),sin(theta)*x,sin(theta)*y,sin(theta)*z],
  where (x,y,z) constitute the unit vector of the rotation axis,
  and theta the half-rotation angle around that axis.
*/
class DLLEXPORT_OST_GEOM Quat:
    private boost::additive<Quat>,
    private boost::multiplicative<Quat, Real>,
    private boost::multiplicative<Quat,Quat>
{
public:
  Quat();

  Quat(Real w, Real x, Real y, Real z);

  // initialize with a rotation (in rad) around a given axis
  Quat(Real angle, const geom::Vec3& axis);

  // initialize from a rotation matrix
  Quat(const Mat3& rotmat);

  // return 3x3 rotation matrix
  Mat3 ToRotationMatrix() const;

  Vec3 GetAxis() const;

  //get angle component
  Real GetAngle() const;

  // operators

  // negateable
  Quat operator-();
  // addable
  Quat& operator+=(const Quat& q);
  // subtractable
  Quat& operator-=(const Quat& q);
  // multipliable with scalar
  Quat& operator*=(Real s);
  // multipliable with other quat
  Quat& operator*=(const Quat& q);
  // dividable with scalar
  Quat& operator/=(Real s);
  // dividable with other quat
  Quat& operator/=(const Quat& q);
  // comparable
  bool operator==(const Quat& q) const;
  bool operator!=(const Quat& q) const { return !this->operator==(q); }

  // Apply rotation to vector.
  Vec3 Rotate(const Vec3& vec) const;

  
  // these are public for direct value access
  Real w,x,y,z;
};


Quat DLLEXPORT_OST_GEOM Conjugate(const Quat& q);

// inner product
Real DLLEXPORT_OST_GEOM Dot(const Quat& q0, const Quat& q1);

// spherical linear interpolation, with t in range [0,1]
Quat DLLEXPORT_OST_GEOM Slerp(const Quat& q0, const Quat& q1, Real t);

Quat DLLEXPORT_OST_GEOM Inv(const Quat& q);
Quat DLLEXPORT_OST_GEOM Exp(const Quat& q);
Quat DLLEXPORT_OST_GEOM Log(const Quat& q);

//normalize quaternion
Quat DLLEXPORT_OST_GEOM Normalize(const Quat& q);

DLLEXPORT_OST_GEOM std::ostream& operator<<(std::ostream& str, const Quat& q);

} // ns

#endif
