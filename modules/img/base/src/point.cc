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
  point

  Authors: Ansgar Philippsen, Andreas Schenk
*/

#include <ost/message.hh>
#include <iostream>
#include <sstream>

#include "point.hh"
#include "size.hh"
#include "vecmat.hh"

namespace ost { namespace img {

Point::Point(const Size& size):
  x(size[0]),
  y(size[1]),
  z(size[2])
{
}

Point Point::Mirror(int planes)
{
  return Point(planes & Plane::YZ ? -x : x,
               planes & Plane::XZ ? -y : y,
               planes & Plane::XY ? -z : z);
}

Point& Point::operator+=(const Point& p) 
{
  x+=p.x;   
  y+=p.y;
  z+=p.z;
  return *this;
}

Point& Point::operator-=(const Point& p) 
{
  x-=p.x;   
  y-=p.y;
  z-=p.z;
  return *this;
}

Point Point::operator-() const 
{
  return Point(-x,-y,-z);
}

Point& Point::operator=(const Point& p) 
{
  x=p.x;
  y=p.y;
  z=p.z;
  return *this;
}


bool Point::operator==(const Point &p) const 
{
  return  (x==p.x && 
           y==p.y && 
           z==p.z);
}

bool Point::operator<(const Point &p) const 
{
  return ( z!=p.z ? z<p.z : (y!=p.y ?  y<p.y : x<p.x) ); 
}


Point& Point::operator+=(const Size& s) 
{
  x+=s[0];   
  y+=s[1];
  z+=s[2];
  return *this;
}

Point& Point::operator-=(const Size& s) 
{
  x-=s[0];   
  y-=s[1];
  z-=s[2];
  return *this;
}

Vec2 Point::ToVec2() const 
{
  Vec2 nrvo(static_cast<Real>(x),
            static_cast<Real>(y));
  return nrvo;
}

Vec3 Point::ToVec3() const 
{
  Vec3 nrvo(static_cast<Real>(x),
      static_cast<Real>(y),
      static_cast<Real>(z));
  return nrvo;
}

Vec4 Point::ToVec4() const 
{
  Vec4 nrvo(static_cast<Real>(x),
      static_cast<Real>(y),
      static_cast<Real>(z),
      1.0);
  return nrvo;
}

// global funcs

std::ostream& operator<<(std::ostream& os, const Point &p) 
{
  os << "(" << p[0] << "," << p[1] << "," << p[2] << ")";
  return os;
}


}} // namespace
