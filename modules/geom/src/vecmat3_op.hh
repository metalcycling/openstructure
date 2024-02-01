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
#ifndef GEOM_VECMAT3_OP_HH
#define GEOM_VECMAT3_OP_HH

#include <ostream>
#include "constants.hh"

#include <ost/geom/module_config.hh>
#include <ost/geom/vec3.hh>
#include <ost/geom/mat3.hh>

namespace geom {


//! returns squared length of vector
inline Real Length2(const Vec3& v)
{
  return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
}

//! returns length of vector
inline Real Length(const Vec3& v)
{
  return std::sqrt(Length2(v));
}

//! return true if components differ not more than epsilon
inline bool Equal(const Vec3& v1, const Vec3& v2, 
                  Real ephilon=EPSILON)
{
  return std::fabs(v1[0]-v2[0])<ephilon &&
    std::fabs(v1[1]-v2[1])<ephilon &&
    std::fabs(v1[2]-v2[2])<ephilon;
}

//! return true if components differ not more than ephilon
inline bool Equal(const Mat3& m1, const Mat3& m2, 
                  Real ephilon=EPSILON)
{
  return std::fabs(m1(0,0)-m2(0,0))<ephilon &&
    std::fabs(m1(0,1)-m2(0,1))<ephilon &&
    std::fabs(m1(0,2)-m2(0,2))<ephilon &&
    std::fabs(m1(1,0)-m2(1,0))<ephilon &&
    std::fabs(m1(1,1)-m2(1,1))<ephilon &&
    std::fabs(m1(1,2)-m2(1,2))<ephilon &&
    std::fabs(m1(2,0)-m2(2,0))<ephilon &&
    std::fabs(m1(2,1)-m2(2,1))<ephilon &&
    std::fabs(m1(2,2)-m2(2,2))<ephilon;
}

//! vector dot product
inline Real Dot(const Vec3& v1, const Vec3& v2)
{
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

//! Normalize
inline Vec3 Normalize(const Vec3& v)
{
  Real l=Length(v);
  if(l==0.0) {
    return v;
  }
  return v/l;
}

//! vector cross product
inline Vec3 Cross(const Vec3& v1, const Vec3& v2)
{
  Vec3 nrvo(v1[1]*v2[2]-v2[1]*v1[2],
      v1[2]*v2[0]-v2[2]*v1[0],
      v1[0]*v2[1]-v2[0]*v1[1]);
  return nrvo;
}

//! multiply each component of v1 with that of v2
inline Vec3 CompMultiply(const Vec3& v1, const Vec3& v2)
{
  Vec3 nrvo(v1[0]*v2[0],v1[1]*v2[1],v1[2]*v2[2]);
  return nrvo;
}

//! divide each component of v1 by that of v2
inline Vec3 CompDivide(const Vec3& v1, const Vec3& v2)
{
  Vec3 nrvo(v1[0]/v2[0],v1[1]/v2[1],v1[2]/v2[2]);
  return nrvo;
}

//! vector matrix multiplication
inline Vec3 operator*(const Vec3& v,const Mat3& m)
{
  Vec3 nrvo(v[0]*m(0,0)+v[1]*m(1,0)+v[2]*m(2,0),
      v[0]*m(0,1)+v[1]*m(1,1)+v[2]*m(2,1),
      v[0]*m(0,2)+v[1]*m(1,2)+v[2]*m(2,2));
  return nrvo;
}

//! vector matrix multiplication
inline Vec3 operator*(const Mat3& m, const Vec3& v)
{
  Vec3 nrvo(v[0]*m(0,0)+v[1]*m(0,1)+v[2]*m(0,2),
      v[0]*m(1,0)+v[1]*m(1,1)+v[2]*m(1,2),
      v[0]*m(2,0)+v[1]*m(2,1)+v[2]*m(2,2));
  return nrvo;
}

//! matrix matrix multiplication
inline Mat3 operator*(const Mat3& m1, const Mat3& m2)
{
  Mat3 nrvo(m1(0,0)*m2(0,0)+m1(0,1)*m2(1,0)+m1(0,2)*m2(2,0),
      m1(0,0)*m2(0,1)+m1(0,1)*m2(1,1)+m1(0,2)*m2(2,1),
      m1(0,0)*m2(0,2)+m1(0,1)*m2(1,2)+m1(0,2)*m2(2,2),

      m1(1,0)*m2(0,0)+m1(1,1)*m2(1,0)+m1(1,2)*m2(2,0),
      m1(1,0)*m2(0,1)+m1(1,1)*m2(1,1)+m1(1,2)*m2(2,1),
      m1(1,0)*m2(0,2)+m1(1,1)*m2(1,2)+m1(1,2)*m2(2,2),

      m1(2,0)*m2(0,0)+m1(2,1)*m2(1,0)+m1(2,2)*m2(2,0),
      m1(2,0)*m2(0,1)+m1(2,1)*m2(1,1)+m1(2,2)*m2(2,1),
      m1(2,0)*m2(0,2)+m1(2,1)*m2(1,2)+m1(2,2)*m2(2,2));
  return nrvo;
}

Mat3 DLLEXPORT_OST_GEOM Invert(const Mat3& m);
Mat3 DLLEXPORT_OST_GEOM Transpose(const Mat3& m);

// algebraic complement
Real DLLEXPORT_OST_GEOM Comp(const Mat3& m, unsigned int i, unsigned int j);

// minor
Real DLLEXPORT_OST_GEOM Minor(const Mat3& m, unsigned int i, unsigned int j);

// determinant
Real DLLEXPORT_OST_GEOM Det(const Mat3& m);

// angle between two vectors
Real DLLEXPORT_OST_GEOM Angle(const Vec3& v1, const Vec3& v2);

// signed angle between two vectors, based on a reference normal
Real DLLEXPORT_OST_GEOM SignedAngle(const Vec3& v1, const Vec3& v2, const Vec3& ref);

Mat3 DLLEXPORT_OST_GEOM EulerTransformation(Real theta, Real phi, Real xi);

Mat3 DLLEXPORT_OST_GEOM AxisRotation(const Vec3& axis, Real angle);

/// \brief get arbitrary vector orthogonal to axis
/// 
/// The returned vector is of length 1
Vec3 DLLEXPORT_OST_GEOM OrthogonalVector(const Vec3& axis);

/// \brief Get dihedral angle for p1-p2-p3-p4
inline Real DihedralAngle(const Vec3& p1, const Vec3& p2,
                          const Vec3& p3, const Vec3& p4) {
  const Vec3 r1 = p2-p1;
  const Vec3 r2 = p3-p2;
  const Vec3 r3 = p4-p3;
  const Vec3 r12cross = Cross(r1, r2);
  const Vec3 r23cross = Cross(r2, r3);
  return std::atan2(Dot(r1*Length(r2), r23cross), Dot(r12cross, r23cross));
}

//! returns std::min of each component
inline Vec3 Min(const Vec3& v1, const Vec3& v2)
{
  Vec3 nrvo(std::min(v1[0],v2[0]),
            std::min(v1[1],v2[1]),
            std::min(v1[2],v2[2]));
  return nrvo;
}

//! returns std::max of each component
inline Vec3 Max(const Vec3& v1, const Vec3& v2)
{
  Vec3 nrvo(std::max(v1[0],v2[0]),
            std::max(v1[1],v2[1]),
            std::max(v1[2],v2[2]));
  return nrvo;
}

//! returns the distance between two points
inline Real Distance(const Vec3& p1, const Vec3& p2)
{
    return Length(p1-p2);
}


//! return the squared distance between two points with periodic boundaries in x,y,z given in ucell_size
inline Real Distance2WithPBC(const Vec3& v1, const Vec3& v2, const Vec3& ucell_size){
  Vec3 v;
  v=v1-v2;
  for (int i=0; i<3; i++) {
    if (std::fabs(v[i])>ucell_size[i]/2.){ 
      v[i]=std::fabs(v[i])-ucell_size[i]*int(std::fabs(v[i])/ucell_size[i]+0.5);
    }
  }
  return Length2(v);
}
//! return the distance between two points with periodic boundaries in x,y,z given in ucell_size
inline Real DistanceWithPBC(const Vec3& v1, const Vec3& v2, const Vec3& ucell_size){
  return sqrt(Distance2WithPBC(v1, v2, ucell_size));
}
//! returns the minimal distance between the points in two Vec3List
Real DLLEXPORT_OST_GEOM MinDistance(const Vec3List& l1, const Vec3List& l2);
//! returns the minimal distance between the points in two Vec3List 
//  with periodic boundaries in x,y,z given in ucell_size
Real DLLEXPORT_OST_GEOM MinDistanceWithPBC(const Vec3List& l1, const Vec3List& l2, Vec3& ucell_size);
//! returns the indices index1, index2 corresponding to the points in
//! the Vec3List l1 and l2 having the minimal distance.
std::vector<unsigned int> DLLEXPORT_OST_GEOM MinDistanceIndices(const Vec3List& l1, const Vec3List& l2);
//! Calculates the Unit Cell Vectors from their sizes and angles (given as Vec3(gamma,beta,alpha)).
Vec3List DLLEXPORT_OST_GEOM CalculateUnitCellVectors(const Vec3& ucell_size, const Vec3& ucell_angles);
//!wraps a vector in a box with periodic boundaries
Vec3 DLLEXPORT_OST_GEOM WrapVec3(const Vec3& v1,const Vec3& box_center,const Vec3& ucell_size);
//!wraps all the vectors in a Vec3List in a box with periodic boundaries
Vec3List DLLEXPORT_OST_GEOM WrapVec3List(const Vec3List& vl,const Vec3& box_center,const Vec3& ucell_size);
//!wraps a vector in a non-rothogonal box with periodic boundaries
Vec3 DLLEXPORT_OST_GEOM WrapVec3(const Vec3& v1,const Vec3& box_center,const Vec3& ucell_size,const Vec3& ucell_angles);
//!wraps all the vectors in a Vec3List in a non-rothogonal box with periodic boundaries
Vec3List DLLEXPORT_OST_GEOM WrapVec3List(const Vec3List& vl,const Vec3& box_center,const Vec3& ucell_size,const Vec3& ucell_angles);
  
} // ns

#endif
