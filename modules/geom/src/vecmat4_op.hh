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
#ifndef GEOM_VECMAT4_OP_HH
#define GEOM_VECMAT4_OP_HH

#include "constants.hh"

#include <ost/geom/module_config.hh>
#include <ost/geom/vec4.hh>
#include <ost/geom/vec3.hh>
#include <ost/geom/mat4.hh>
#include <ost/geom/mat3.hh>

namespace geom {


//! returns squared length of vector
inline Real Length2(const Vec4& v)
{
    return v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3];
}

//! returns length of vector
inline Real Length(const Vec4& v)
{
  return std::sqrt(Length2(v));
}

//! return true if components differ not more than ephilon
inline bool Equal(const Vec4& v1, const Vec4& v2, Real epsilon=EPSILON)
{
  return std::fabs(v1[0]-v2[0])<epsilon &&
    std::fabs(v1[1]-v2[1])<epsilon &&
    std::fabs(v1[2]-v2[2])<epsilon &&
    std::fabs(v1[3]-v2[3])<epsilon;
}

//! return true if components differ not more than ephilon
inline bool Equal(const Mat4& m1, const Mat4& m2, Real epsilon=EPSILON)
{
  return std::fabs(m1(0,0)-m2(0,0))<epsilon &&
    std::fabs(m1(0,1)-m2(0,1))<epsilon &&
    std::fabs(m1(0,2)-m2(0,2))<epsilon &&
    std::fabs(m1(0,3)-m2(0,3))<epsilon &&
    std::fabs(m1(1,0)-m2(1,0))<epsilon &&
    std::fabs(m1(1,1)-m2(1,1))<epsilon &&
    std::fabs(m1(1,2)-m2(1,2))<epsilon &&
    std::fabs(m1(1,3)-m2(1,3))<epsilon &&
    std::fabs(m1(2,0)-m2(2,0))<epsilon &&
    std::fabs(m1(2,1)-m2(2,1))<epsilon &&
    std::fabs(m1(2,2)-m2(2,2))<epsilon &&
    std::fabs(m1(2,3)-m2(2,3))<epsilon &&
    std::fabs(m1(3,0)-m2(3,0))<epsilon &&
    std::fabs(m1(3,1)-m2(3,1))<epsilon &&
    std::fabs(m1(3,2)-m2(3,2))<epsilon &&
    std::fabs(m1(3,3)-m2(3,3))<epsilon;
}

//! vector dot product
inline Real Dot(const Vec4& v1, const Vec4& v2)
{
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3];
}

inline Vec4 Normalize(const Vec4& v)
{
  Real l=Length(v);
  if(l==0.0) {
    return v;
  }
  return v/l;
}

//! multiply each component of v1 with that of v2
inline Vec4 CompMultiply(const Vec4& v1, const Vec4& v2)
{
  Vec4 nrvo(v1[0]*v2[0],v1[1]*v2[1],v1[2]*v2[2],v1[3]*v2[3]);
  return nrvo;
}

//! divide each component of v1 by that of v2
inline Vec4 CompDivide(const Vec4& v1, const Vec4& v2)
{
  Vec4 nrvo(v1[0]/v2[0],v1[1]/v2[1],v1[2]/v2[2],v1[3]/v2[3]);
  return nrvo;
}

//! vector matrix multiplication
inline Vec4 operator*(const Vec4& v,const Mat4& m)
{
  Vec4 nrvo(v[0]*m(0,0)+v[1]*m(1,0)+v[2]*m(2,0)+v[3]*m(3,0),
      v[0]*m(0,1)+v[1]*m(1,1)+v[2]*m(2,1)+v[3]*m(3,1),
      v[0]*m(0,2)+v[1]*m(1,2)+v[2]*m(2,2)+v[3]*m(3,2),
      v[0]*m(0,3)+v[1]*m(1,3)+v[2]*m(2,3)+v[3]*m(3,3));
  return nrvo;
}

//! vector matrix multiplication
inline Vec4 operator*(const Mat4& m, const Vec4& v)
{
  Vec4 nrvo(v[0]*m(0,0)+v[1]*m(0,1)+v[2]*m(0,2)+v[3]*m(0,3),
      v[0]*m(1,0)+v[1]*m(1,1)+v[2]*m(1,2)+v[3]*m(1,3),
      v[0]*m(2,0)+v[1]*m(2,1)+v[2]*m(2,2)+v[3]*m(2,3),
      v[0]*m(3,0)+v[1]*m(3,1)+v[2]*m(3,2)+v[3]*m(3,3));
  return nrvo;
}

// algebraic complement
Real DLLEXPORT_OST_GEOM Comp(const Mat4& m, unsigned int i, unsigned int j);

// minor
Real DLLEXPORT_OST_GEOM Minor(const Mat4& m, unsigned int i, unsigned int j);

Real DLLEXPORT_OST_GEOM Det(const Mat4& m);
Mat4 DLLEXPORT_OST_GEOM Transpose(const Mat4& m);
Mat4 DLLEXPORT_OST_GEOM Invert(const Mat4& m);
Mat4 DLLEXPORT_OST_GEOM operator*(const Mat4& m1, const Mat4& m2);
Real DLLEXPORT_OST_GEOM Angle(const Vec4& v1, const Vec4& v2);

//! returns std::min of each component
inline Vec4 Min(const Vec4& v1, const Vec4& v2)
{
  Vec4 nrvo(std::min(v1[0],v2[0]),
            std::min(v1[1],v2[1]),
            std::min(v1[2],v2[2]),
            std::min(v1[3],v2[3]));
  return nrvo;
}

//! returns std::max of each component
inline Vec4  Max(const Vec4& v1, const Vec4& v2)
{
  Vec4 nrvo(std::max(v1[0],v2[0]),
            std::max(v1[1],v2[1]),
            std::max(v1[2],v2[2]),
            std::max(v1[3],v2[3]));
  return nrvo;
}

} // ns

#endif
