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
#ifndef GEOM_COMPOSITE3_HH
#define GEOM_COMPOSITE3_HH

#include <iostream>
#include <vector>
#include <boost/filesystem/path.hpp>

#include "vec3.hh"
#include "mat3.hh"
#include "quat.hh"

#include "module_config.hh"

/*
  composite classes in 3D space
*/

namespace geom {

//! Line3
class DLLEXPORT_OST_GEOM Line3 {
public:
  Line3();
  Line3(const Vec3& from, const Vec3& to);

  Vec3 At(Real r) const;
  Vec3 GetOrigin() const;
  Vec3 GetDirection() const;
private:
  Vec3 ori_,dir_;
};

DLLEXPORT_OST_GEOM std::ostream& operator<<(std::ostream& s, const Line3& l);

class DLLEXPORT_OST_GEOM Plane {
public:
  Plane():n_(Vec3(0,0,1)),p_(0){};
  //! initialize based on 3 points defining a plane
  Plane(const Vec3& p1, const Vec3& p2, const Vec3& p3);
  //! initialize based on a point and a normal
  Plane(const Vec3& p, const Vec3& n);
  //! initialize based on a line and point
  Plane(const Line3& l, const Vec3& p);
  //! initialize a,b,c,d (ax+by+cz+d=0)
  Plane(Real a, Real b, Real c, Real d);
  //! initialize based x-,y-,z- intercepts
  Plane(Real x, Real y, Real z);

  Vec3 GetOrigin() const { return -p_*n_; }
  Vec3 GetNormal() const;
  void SetNormal(const Vec3& n);
  Real GetP()    const;
  void SetP(Real p) { p_=p; }
  Vec3 At(Real x, Real y) const;

  enum PLANE_TYPE {
    YZ=0x1, ZY=0x1, XZ=0x2, ZX=0x2, XY=0x4, YX=0x4
  };
  
private:
  // planes are internally specified in HNF 
  Vec3 n_;
  Real p_;
};

class DLLEXPORT_OST_GEOM Sphere {
public:
  Sphere();
  Sphere(const Vec3& origin, Real r);
  Vec3 GetOrigin() const;
  Real GetRadius() const;
  void SetOrigin(const Vec3& v);
  void SetRadius(Real r);

private:
  Vec3    origin_;
  Real radius_;
};

/// \brief cuboid axis defined by a normalized direction vector and a
///     half extent
/// 
/// \relates Cuboid
class DLLEXPORT_OST_GEOM CuboidAxis {
public:
  CuboidAxis(): axis_(), half_extent_(0.0)
  { }
  CuboidAxis(const Vec3& axis, Real half_extent);
  const Vec3& GetVector() const { return axis_; }
  Real GetHalfExtent() const { return half_extent_; }
  Real GetExtent() const { return 2.0*half_extent_; }
private:
  Vec3       axis_;
  Real     half_extent_;
};

/// \brief arbitrary oriented bounding cuboid
class DLLEXPORT_OST_GEOM Cuboid {
public:
  Cuboid();
  Cuboid(const Vec3& center, const CuboidAxis& a, 
              const CuboidAxis& b, const CuboidAxis& c);

  Vec3 GetHalfExtents() const 
  { 
    return Vec3(axes_[0].GetHalfExtent(), axes_[1].GetHalfExtent(), 
                      axes_[2].GetHalfExtent());
  }
  
  Vec3 GetCenter() const 
  {
    return center_;
  }
  const Vec3& GetVecA() const { return axes_[0].GetVector(); }
  const Vec3& GetVecB() const { return axes_[1].GetVector(); }  
  const Vec3& GetVecC() const { return axes_[2].GetVector(); }
  
  const CuboidAxis& GetAxisA() const { return axes_[0]; }
  const CuboidAxis& GetAxisB() const { return axes_[1]; }  
  const CuboidAxis& GetAxisC() const { return axes_[2]; }  
private:
  Vec3   center_;
  CuboidAxis      axes_[3];
};

class DLLEXPORT_OST_GEOM Rotation3
{
public:
  Rotation3();
  Rotation3(Real phi, Real theta, Real psi,
            const Vec3& origin=Vec3(0.0,0.0,0.0));
  Rotation3(const Vec3& axis, Real angle,
            const Vec3& origin=Vec3(0.0,0.0,0.0));
  Rotation3(const Line3& line, Real angle);
  Rotation3(const Mat3& rot, const Vec3& origin=Vec3(0.0,0.0,0.0));
  Rotation3(const Quat& q, const Vec3& origin=Vec3(0.0,0.0,0.0));
  Vec3 GetOrigin() const;
  Real GetPhi() const;
  Real GetTheta() const;
  Real GetPsi() const;
  Quat GetQuat() const;
  Vec3 GetRotationAxis() const;
  Real GetRotationAngle() const;
  Mat3 GetRotationMatrix() const;
  void SetOrigin(const Vec3& o);
  void SetPhi(Real phi);
  void SetTheta(Real theta);
  void SetPsi(Real psi);
  void SetQuat(const Quat& q);
  void SetRotationAxis(const Vec3& v);
  void SetRotationAngle(Real angle);
  void SetRotationMatrix(const Mat3& rot);
  Vec3 Apply(const Vec3&) const;
  bool operator==(const Rotation3& rhs) const;

private:
  Quat generate_from_eulers(Real psi, Real theta, Real phi);
  Quat generate_from_axis_angle(const Vec3& axis, Real angle);
  Vec3 find_invariant_vector(Mat3 rot);
  Vec3 find_orthogonal_vector(const Vec3& xyz);
  Vec3 find_vector_for_BOD(const Vec3& xyz, const Vec3& uvw);
  Quat generate_from_matrix(const Mat3& rot);
  Real reduce_angle(Real ang) const;
  Quat q_;
  Vec3 origin_;
};

typedef std::vector<Rotation3> Rotation3List;

DLLEXPORT_OST_GEOM Rotation3List ImportEulerAngles (const boost::filesystem::path& loc);
DLLEXPORT_OST_GEOM void ExportEulerAngles (const Rotation3List& rot_list,
                                  const boost::filesystem::path& loc);

} // ns

#endif
