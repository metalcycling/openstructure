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
#ifndef GEOM_VEC4_H
#define GEOM_VEC4_H

#include <stdexcept>
#include <cassert>
#include <cstddef> // for size_t
#include <ostream>

#include <boost/operators.hpp>

#include <ost/geom/module_config.hh>
#include <ost/config.hh>

namespace geom {

// fw decl
class Vec2;
class Vec3;

/*
  Four dimensional, homogeneous vector class, using Real precision.
*/
class DLLEXPORT Vec4:
    private boost::equality_comparable<Vec4>,
    private boost::additive<Vec4>,
    private boost::additive<Vec4, Real>,
    private boost::multiplicative<Vec4, Real>
{
public:
  //! Default initialization, all components are set to zero
  Vec4(): x(0), y(0), z(0), w(1) { }

  //! Initialization with x, y and z component
  Vec4(Real px, Real py, Real pz, Real pw): x(px), y(py), z(pz), w(pw) { }

  //! copy ctor
  Vec4(const Vec4& v): x(v[0]), y(v[1]), z(v[2]), w(v[3]) { }

  //! (implicit) initialization with 2D vector
  Vec4(const Vec2& v);

  //! (implicit) initialization with 3D vector
  Vec4(const Vec3& v);

  //! explicit initialization with an array of floats
  explicit Vec4(const float v[4]): x(v[0]), y(v[1]), z(v[2]), w(v[3]) { }
  
  //! explicit initialization with an array of doubles  
  explicit Vec4(const double v[4]): x(v[0]), y(v[1]), z(v[2]), w(v[3]) { }
  //! assignement op
  Vec4& operator=(const Vec4& v)
  {
    x=v.x;
    y=v.y;
    z=v.z;
    w=v.w;
    return *this;
  }

  //! comparable
  bool operator==(const Vec4& rhs) const
  {
    return x==rhs.x && y==rhs.y && z==rhs.z && w==rhs.w;
  }

  //! element access
  Real& operator[](std::size_t indx)
  {
    assert(indx<4);
    return (&x)[indx];
  }
  
  //! const element access
  const Real& operator[](std::size_t indx) const
  {
    assert(indx<4);
    return (&x)[indx];
  }

  Real& At(size_t indx) {
    if (indx>3) {
      throw std::out_of_range("index must be smaller than 4");
    }
    return (&x)[indx];
  }

  const Real& At(size_t indx) const {
    if (indx>3) {
      throw std::out_of_range("index must be smaller than 4");
    }
    return (&x)[indx];
  }

  //! element access
  Real GetX() const { return x; }
  Real GetY() const { return y; }
  Real GetZ() const { return z; }
  Real GetW() const { return w; }
  
  void SetX(Real v) { x=v; }
  void SetY(Real v) { y=v; }
  void SetZ(Real v) { z=v; }
  void SetW(Real v) { w=v; }  

  //! addable op
  Vec4& operator+=(const Vec4& rhs)
  {
    x+=rhs.x;
    y+=rhs.y;
    z+=rhs.z;
    w+=rhs.w;
    return *this;
  }
  
  Vec4& operator+=(Real d)
  {
    x+=d;
    y+=d;
    z+=d;
    w+=d;
    return *this;
  }
  
  //! subtractable op
  Vec4& operator-=(const Vec4& rhs)
  {
    x-=rhs.x;
    y-=rhs.y;
    z-=rhs.z;
    w-=rhs.w;
    return *this;
  }
  
  Vec4& operator-=(Real d)
  {
    x-=d;
    y-=d;
    z-=d;
    w-=d;
    return *this;
  }
  //! negateable
  Vec4 operator-() const
  {
    return Vec4(-x, -y, -z, -w);
  }
  
  //! multipliable
  Vec4& operator*=(Real d)
  {
    x*=d;
    y*=d;
    z*=d;
    w*=d;
    return *this;
  }
  
  //! dividable
  Vec4& operator/=(Real d)
  {
    Real one_over_d=Real(1.0)/d;
    x*=one_over_d;
    y*=one_over_d;
    z*=one_over_d;
    w*=one_over_d;
    return *this;
  }

  Real* Data() {return &x;}
  const Real* Data() const {return &x;}

  Real x;
  Real y;
  Real z;
  Real w;
};

inline Vec4 operator/(Real d, const Vec4& v)
{
  Vec4 nrvo(d/v[0],d/v[1],d/v[2],d/v[3]);
  return nrvo;
}

inline std::ostream& operator<<(std::ostream& os, const Vec4& v)
{
  os << "(" << v[0] << "," << v[1] << "," << v[2] << "," << v[3] << ")";
  return os;
}
}

#include <ost/geom/vec2.hh>
#include <ost/geom/vec3.hh>

namespace geom {
  
//! (implicit) initialization with 2D vector
inline Vec4::Vec4(const Vec2& v): x(v.x), y(v.y), z(0), w(1) { }

//! (implicit) initialization with 3D vector
inline Vec4::Vec4(const Vec3& v): x(v.x), y(v.y), z(v.z), w(1.0) { }

} // namespace geom

# endif
