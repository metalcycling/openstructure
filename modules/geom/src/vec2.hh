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
#ifndef GEOM_VEC2_H
#define GEOM_VEC2_H

#include <stdexcept>
#include <cassert>
#include <cstddef> // for size_t
#include <ostream>
#include <vector>
#include <boost/operators.hpp>


#include <ost/config.hh>
#include <ost/geom/module_config.hh>
#include <ost/geom/exc.hh>
namespace geom {

// fw decl
class Vec3;
class Vec4;

/*
  Two dimensional vector class, using Real precision.
*/
class DLLEXPORT Vec2:
    private boost::equality_comparable<Vec2>,
    private boost::additive<Vec2>,
    private boost::additive<Vec2, Real>,
    private boost::multiplicative<Vec2, Real>
{
public:
  //! Default initialization, all components are set to zero
  Vec2(): x(0), y(0) { }

  //! Initialization with x, y and z component
  Vec2(Real px, Real py): x(px), y(py) { }

  //! copy ctor
  Vec2(const Vec2& v): x(v.x), y(v.y) { }

  //! explicit initialization with 3D vector
  explicit Vec2(const Vec3& v);

  //! explicit initialization with a 4D homogeneous vector
  explicit Vec2(const Vec4& v);

  //! explicit initialization with an array of doubles
  explicit Vec2(const float v[2]): x(v[0]), y(v[1]) { }


  //! explicit initialization with an array of floats
  explicit Vec2(const double v[2]): x(v[0]), y(v[1]) { }


  Real GetX() const { return x; }
  Real GetY() const { return y; }

  void SetX(Real d) { x=d; }
  void SetY(Real d) { y=d; }

  //! comparable
  bool operator==(const Vec2& rhs) const
  {
    return x==rhs.x && y==rhs.y;
  }

  //! element access
  Real& operator[](std::size_t indx)
  {
    assert(indx<2);
    return (&x)[indx];
  }
  
  //! const element access
  const Real& operator[](std::size_t indx) const
  {
    assert(indx<2);
    return (&x)[indx];
  }

  Real& At(size_t indx) {
    if (indx>1) {
      throw std::out_of_range("index must be smaller than 2");
    }
    return (&x)[indx];
  }

  const Real& At(size_t indx) const {
    if (indx>1) {
      throw std::out_of_range("index must be smaller than 2");
    }
    return (&x)[indx];
  }
  
  //! addable op
  Vec2& operator+=(const Vec2& rhs)
  {
    x+=rhs.x;
    y+=rhs.y;
    return *this;
  }
  
  Vec2& operator+=(Real d)
  {
    x+=d;
    y+=d;
    return *this;
  }
  
  //! subtractable op
  Vec2& operator-=(const Vec2& rhs)
  {
    x-=rhs.x;
    y-=rhs.y;
    return *this;
  }
  
  Vec2& operator-=(Real d)
  {
    x-=d;
    y-=d;
    return *this;
  }
  //! negateable
  Vec2 operator-() const
  {
    return Vec2(-x, -y);
  }

  //! multipliable
  Vec2& operator*=(Real d)
  {
    x*=d;
    y*=d;
    return *this;
  }

  //! dividable
  Vec2& operator/=(Real d)
  {
    Real one_over_d=Real(1.0)/d;
    x*=one_over_d;
    y*=one_over_d;
    return *this;
  }

  Real* Data() {return &x;}
  const Real* Data() const { return &x; }

  Real x;
  Real y;
};

inline Vec2 operator/(Real d, const Vec2& v)
{
  return Vec2(d/v.x, d/v.y);
}

inline std::ostream& operator<<(std::ostream& os, const Vec2& v)
{
  os << "(" << v[0] << "," << v[1] << ")";
  return os;  
}

}

#include <ost/geom/vec3.hh>
#include <ost/geom/vec4.hh>

namespace geom {

inline Vec2::Vec2(const Vec3& v): x(v.x), y(v.y) { }


inline Vec2::Vec2(const Vec4& v): x(v.x), y(v.y) 
{
  if (std::fabs(v.w)<1e-10) {
    throw DivideByZeroException();
  }
  x/=v.w;
  y/=v.w;
}


typedef std::vector<Vec2> Vec2List;

} // namespace geom

# endif
