//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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

/*
  integer triplet based point

  Authors: Ansgar Philippsen, Andreas Schenk
*/

#ifndef IMG_POINT_H
#define IMG_POINT_H

#include <vector>
#include <iosfwd>
#include <boost/operators.hpp>

#include <ost/img/module_config.hh>
#include "vecmat.hh"

namespace ost { namespace img {

//fw decl
class Size;

//! class encapsulating 1D to 3D point
class DLLEXPORT_OST_IMG_BASE Point: boost::additive<Point,
                                    boost::additive2<Point, Size,
                                    boost::multipliable2<Point,int,
                                    boost::less_than_comparable<Point,
                                    boost::equality_comparable<Point> > > > >{
public:
  Point():
    x(),
    y(),
    z()
  {
  }
  
  Point(const Point &p):
    x(p.x),
    y(p.y),
    z(p.z)
  {
  } 


  //! 1D int constructor
  Point(int a):
    x(a),
    y(),
    z()
  {
  } 

  //! 2D int constructor
  Point(int a, int b):
    x(a),
    y(b),
    z()
  {
  } 

  //! 3D int constructor
  Point(int a, int b, int c):
    x(a),
    y(b),
    z(c)
  {
  }
  
  //! conversion from Vec2
  explicit Point(const Vec2& v):
    x(static_cast<int>(round(v[0]))),
    y(static_cast<int>(round(v[1]))),
    z()
  {
  } 


  //! conversion from Vec3
  explicit Point(const Vec3& v):
    x(static_cast<int>(round(v[0]))),
    y(static_cast<int>(round(v[1]))),
    z(static_cast<int>(round(v[2])))
  {
  } 


  //! conversion from Vec4 (normalization)
  explicit Point(const Vec4& v):
    x(static_cast<int>(round(v[0]))),
    y(static_cast<int>(round(v[1]))),
    z(static_cast<int>(round(v[2])))
  {
    if(std::abs(v[3])<1e-100) {
      throw geom::OutOfRangeException("4th element of Vec4 is too close to zero for normalization");
    } else {
      Real sf = 1.0/v[3];
      x = static_cast<int>(x*sf);
      y = static_cast<int>(y*sf);
      z = static_cast<int>(z*sf);
    }
  }

  //! (implicit) conversion of size to point
  Point(const Size& size);

  //! return mirror point according to planes
  Point Mirror(int planes);

  // operators
  int& operator[](unsigned int index) 
  {
    assert(index<=2);
    return (&x)[index];
  }
  
  int operator[](unsigned int index) const 
  {
    assert(index<=2);
    return (&x)[index];
  }
  // operators

  int& At(unsigned int index) 
  {
    if(index>2) throw std::range_error("Point index out of range");
    return (&x)[index];
  }
  
  int At(unsigned int index) const 
  {
    if(index>2) throw geom::OutOfRangeException("Point index out of range");
    return (&x)[index];
  }


  Point& operator=(const Point& p);

  Point& operator+=(const Point& p);

  Point& operator-=(const Point& p);

  Point operator-() const;

  Point& operator+=(const Size& p);

  Point& operator-=(const Size& p);

  bool operator==(const Point &p) const;

  bool operator<(const Point &p) const;

  Point& operator*=(int s) {x*=s;y*=s;z*=s;return *this;}
  
  // conversion to vectors
  Vec2 ToVec2() const;
  Vec3 ToVec3() const;
  Vec4 ToVec4() const;

  int x;
  int y;
  int z;
};


DLLEXPORT_OST_IMG_BASE std::ostream& operator<<(std::ostream& os, const Point &p);


}} // namespace img


#endif



