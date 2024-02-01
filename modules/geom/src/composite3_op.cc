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
#include <cmath>
#include <iostream>
#include "exc.hh"
#include "composite3_op.hh"
#include "vecmat3_op.hh"

namespace geom {


Real Angle(const Line3& l1, const Line3& l2)
{
  return Angle(l1.GetDirection(),l2.GetDirection());
}
Real Angle(const Plane& p1, const Plane& p2)
{
  return Angle(p1.GetNormal(),p2.GetNormal());
}
Real Angle(const Plane& p, const Line3& l)
{
  return M_PI/2.0-Angle(p.GetNormal(),l.GetDirection());
}


Vec3 IntersectionPoint(const Line3& l1, const Line3& l2)
{
  if(!AreIntersecting(l1,l2)) throw GeomException("Lines are not intersecting!");
  Vec3 d=l2.GetOrigin()-l1.GetOrigin();
  Vec3 d2=Cross(l1.GetDirection(),l2.GetDirection());
  return l1.GetOrigin()+l1.GetDirection()* Dot(Cross(d,l2.GetDirection()),d2) / (Length(d2)*Length(d2));
}

Vec3 IntersectionPoint(const Line3& l, const Plane& p)
{
  if(AreParallel(p,l)) throw GeomException("Plane and line are parallel!");
  Real scal=Dot(p.GetNormal(),l.GetDirection());
  Real t=Dot(p.GetNormal(),p.GetOrigin()-l.GetOrigin())/scal;
  /*
    How many times the direction vector of the line has to 
    be stretched to meet the plane
  */
  return l.At(t);
}

Vec3 IntersectionPoint(const Plane& p, const Line3& l)
{
  return IntersectionPoint(l,p);
}

Line3 IntersectionLine(const Plane& p1, const Plane& p2)
{
  if(AreParallel(p1,p2)) throw GeomException("Lines are parallel!");
  Vec3 dir=Cross(p1.GetNormal(),p2.GetNormal());
  dir=dir/Length(dir);
  Vec3 tmp=Cross(p1.GetNormal(),dir);
  Vec3 ori=IntersectionPoint(Line3(p1.GetOrigin(),p1.GetOrigin()+tmp),p2);
  return Line3(ori,ori+dir);
}

Real Distance(const Line3& l1, const Line3& l2)
{
  if(AreParallel(l1,l2)) return 0;
  Vec3 d=l2.GetOrigin()-l1.GetOrigin();
  return std::fabs(Length(d)*Length(Cross(l1.GetDirection(),l2.GetDirection()))/Length(Cross(l1.GetDirection(),l2.GetDirection())));
}
Real Distance(const Line3& l, const Vec3& v)
{
  if(Length2(l.GetDirection())==0.0) return 0.0;
  return Length(Cross(l.GetDirection(),v-l.GetOrigin()))/Length(l.GetDirection());
}

Real Distance(const Plane& p, const Vec3& v)
{
  return Dot(p.GetNormal(),v)+p.GetP();
}

bool IsOnLine(const Line3& l, const Vec3& v, Real ephilon)
{
  return Length2(Cross(v-l.GetOrigin(), l.GetDirection()))<(ephilon*ephilon);
}

bool IsInPlane(const Plane& p,  const Line3& l,Real ephilon)
{
  return Distance(p,l.GetOrigin())<ephilon && AreParallel(p,l,ephilon);  
}
bool IsInPlane(const Plane& p,  const Vec3& v,Real ephilon)
{
  return std::fabs(Distance(p,v))<ephilon;  
}
int RelativePosition(const Plane& p,  const Vec3& v,Real ephilon)
{
  // 1=above,0=in plane,-1=below
  Real d=Distance(p,v);
  if (std::fabs(d)<ephilon) return 0;
  if (std::fabs(d)<ephilon) return 1;
  return -1;  
}


bool operator==(const Plane& p1, const Plane& p2)
{
  return Equal(p1,p2);
}
bool Equal(const Plane& p1, const Plane& p2,Real ephilon)
{
  return std::fabs(p1.GetP()-p2.GetP())<ephilon && 
         Equal(p1.GetNormal(),p2.GetNormal(),ephilon);
}
bool EqualPosition(const Plane& p1, const Plane& p2,Real ephilon)
{
  return  Equal(p1,p2,ephilon) || 
          (std::fabs(p1.GetP()+p2.GetP())<ephilon && Equal(p1.GetNormal(),
                                                           -p2.GetNormal(),
                                                           ephilon));
}
bool operator==(const Line3& l1, const Line3& l2)
{
  return Equal(l1,l2);
}

bool Equal(const Line3& l1, const Line3& l2,Real ephilon)
{
  return Equal(l1.GetDirection(), l2.GetDirection(), ephilon) &&
         Equal(l1.GetOrigin(), l2.GetOrigin(), ephilon);
}

bool EqualPosition(const Line3& l1, const Line3& l2,Real ephilon)
{
  return IsOnLine(l1,l2.GetOrigin(), ephilon) && AreParallel(l1,l2, ephilon);
}


bool AreParallel(const Plane& p,  const Line3& l,Real ephilon)
{
  return std::fabs(Dot(l.GetDirection(),p.GetNormal()))<ephilon;
}

bool AreParallel(const Plane& p1, const Plane& p2,Real ephilon)
{
  return std::fabs(1-Dot(p1.GetNormal(),p2.GetNormal()))<ephilon;
}

bool AreParallel(const Line3& l1, const Line3& l2, Real ephilon)
{
  return std::fabs(1.0-std::fabs(Dot(l1.GetDirection(), 
                                     l2.GetDirection())))<ephilon;
}


bool AreIntersecting(const Line3& l1, const Line3& l2, Real ephilon)
{
  if (AreParallel(l1, l2, ephilon)) {
    return false;
  }
  Plane plane(l1.GetOrigin(), Cross(l1.GetDirection(), 
                                    l2.GetDirection()));
  return IsInPlane(plane, l2.GetOrigin(), ephilon);
}

bool IsInSphere(const Sphere& s, const Vec3& v){
  return Length(s.GetOrigin()-v)<=s.GetRadius();
}

} // ns

