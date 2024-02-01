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
#ifndef GEOM_VECMAT2_OP_HH
#define GEOM_VECMAT2_OP_HH

#include "constants.hh"

#include <ost/geom/module_config.hh>

#include <ost/geom/vec2.hh>
#include <ost/geom/mat2.hh>

namespace geom {

class Vec2;
class Mat2;

//! returns squared length of vector
inline Real Length2(const Vec2& v)
{
  return v[0]*v[0]+v[1]*v[1];
}

//! returns length of vector
inline Real Length(const Vec2& v)
{
  return std::sqrt(Length2(v));
}

//! return true if components differ not more than ephilon
inline bool Equal(const Vec2& v1, const Vec2& v2, Real ephilon=EPSILON)
{
  return std::fabs(v1[0]-v2[0])<ephilon &&
    std::fabs(v1[1]-v2[1])<ephilon;
}

//! return true if components differ not more than ephilon
inline bool Equal(const Mat2& m1, const Mat2& m2, Real ephilon=EPSILON)
{
  return std::fabs(m1(0,0)-m2(0,0))<ephilon &&
    std::fabs(m1(0,1)-m2(0,1))<ephilon &&
    std::fabs(m1(1,0)-m2(1,0))<ephilon &&
    std::fabs(m1(1,1)-m2(1,1))<ephilon;
}

//! vector dot product
inline Real Dot(const Vec2& v1, const Vec2& v2)
{
    return v1[0]*v2[0]+v1[1]*v2[1];
}

//! Normalize vector
inline Vec2 Normalize(const Vec2& v)
{
  Real l=Length(v);
  if(l==0.0) {
    return v;
  }
  return v/l;
}

//! multiply each component of v1 with that of v2
inline Vec2 CompMultiply(const Vec2& v1, const Vec2& v2)
{
  Vec2 nrvo(v1[0]*v2[0],v1[1]*v2[1]);
  return nrvo;
}

//! divide each component of v1 by that of v2
inline Vec2 CompDivide(const Vec2& v1, const Vec2& v2)
{
  Vec2 nrvo(v1[0]/v2[0],v1[1]/v2[1]);
  return nrvo;
}

//! vector matrix multiplication
/*!
  This corresponds to

          |a b|
  (x,y) * |c d| = (ax+cy, bx+dy)
*/
inline Vec2 operator*(const Vec2& v,const Mat2& m)
{
  Vec2 nrvo(v[0]*m(0,0)+v[1]*m(1,0),
      v[0]*m(0,1)+v[1]*m(1,1));
  return nrvo;
}

//! vector matrix multiplication
/*!
  This corresponds to

   |a b|   (x)   (ax+by)
   |c d| * (y)  =(cx+dy)
*/
inline Vec2 operator*(const Mat2& m, const Vec2& v)
{
  Vec2 nrvo(v[0]*m(0,0)+v[1]*m(0,1),
      v[0]*m(1,0)+v[1]*m(1,1));
  return nrvo;
}

//! determinant
Real DLLEXPORT_OST_GEOM Det(const Mat2& m);
//! Transpose
Mat2 DLLEXPORT_OST_GEOM Transpose(const Mat2& m);

//! Matrix inversion
Mat2 DLLEXPORT_OST_GEOM Invert(const Mat2& m);

//! angle beetwen two vectors
Real DLLEXPORT_OST_GEOM Angle(const Vec2& v1, const Vec2& v2);
//! angle beetwen two vectors (honors sign)
Real  DLLEXPORT_OST_GEOM SignedAngle(const Vec2& v1, const Vec2& v2);

//! matrix matrix multiplication
inline Mat2 operator*(const Mat2& m1, const Mat2& m2)
{
  Mat2 nrvo(m1(0,0)*m2(0,0)+m1(0,1)*m2(1,0),
      m1(0,0)*m2(0,1)+m1(0,1)*m2(1,1),
      m1(1,0)*m2(0,0)+m1(1,1)*m2(1,0),
      m1(1,0)*m2(0,1)+m1(1,1)*m2(1,1));
  return nrvo;
}

inline Vec2 Min(const Vec2& v1, const Vec2& v2)
{
  Vec2 nrvo(std::min(v1[0],v2[0]),
            std::min(v1[1],v2[1]));
  return nrvo;
}

inline Vec2 Max(const Vec2& v1, const Vec2& v2)
{
  Vec2 nrvo(std::max(v1[0],v2[0]),
            std::max(v1[1],v2[1]));
  return nrvo;
}

//! Rotation
DLLEXPORT Vec2 Rotate(const Vec2& v,Real ang);

} // namespace geom

#endif
