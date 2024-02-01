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
#include <ost/geom/geom.hh>


#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "helper.hh"
using namespace geom;

BOOST_AUTO_TEST_SUITE( geom );

BOOST_AUTO_TEST_CASE(init_composite3)
{
  
  // Vec3 v1(3,4,5);
  // Vec3 v2(8,10,-2);
  // Vec3 v3(-5,3,8);
  
  Vec3 v1(0,5,0);
  Vec3 v2(1,5,0);
  Vec3 v3(0,5,1);
  //------------- Line3 --------------
  //default
  Line3 l;
  BOOST_CHECK(match(l.GetOrigin(),0.0,0.0,0.0));
  BOOST_CHECK(match(l.GetDirection(),1.0,0.0,0.0));
      
  //2 Points
  l=Line3(v1,v2);
  BOOST_CHECK(Normalize(v2-v1)==l.GetDirection());

  //------------- Plane --------------
  Plane p;
  BOOST_CHECK(match(p.GetOrigin(),0.0,0.0,0.0));
  BOOST_CHECK(match(p.GetNormal(),0.0,0.0,1.0));
  BOOST_CHECK(p.GetP()==0.0);
  
  
  p=Plane(v1,v2,v3);
  BOOST_CHECK(dequal(Length(p.GetNormal()), 1.0));
  BOOST_CHECK(IsInPlane(p, p.GetOrigin()));
  BOOST_CHECK(IsInPlane(p, v1));
  BOOST_CHECK(IsInPlane(p, v2));
  BOOST_CHECK(IsInPlane(p, v3));
  BOOST_CHECK(std::fabs(Dot(v2-v1,p.GetNormal()))<EPSILON);
  BOOST_CHECK(std::fabs(Dot(v3-v1,p.GetNormal()))<EPSILON);
  BOOST_CHECK(std::fabs(Dot(v3-v2,p.GetNormal()))<EPSILON);

  BOOST_CHECK(std::fabs(Dot(v1-p.GetOrigin(),p.GetNormal()))<EPSILON);
  BOOST_CHECK(std::fabs(Dot(v2-p.GetOrigin(),p.GetNormal()))<EPSILON);
  BOOST_CHECK(std::fabs(Dot(v3-p.GetOrigin(),p.GetNormal()))<EPSILON);
  
  //initialize based on a point and a normal
  p=Plane(v1,v1);
  BOOST_CHECK(v1==p.GetOrigin());
  BOOST_CHECK(v1/Length(v1)==p.GetNormal());
  
  // initialize based on a line and point
  p=Plane(Line3(v1,v2),v3);
  BOOST_CHECK(std::fabs(Dot(v2-v1,p.GetNormal()))<EPSILON);
  BOOST_CHECK(std::fabs(Dot(v3-v1,p.GetNormal()))<EPSILON);
  BOOST_CHECK(std::fabs(Dot(v3-v2,p.GetNormal()))<EPSILON);
  BOOST_CHECK(std::fabs(Dot(p.GetOrigin()-v1,p.GetNormal()))<EPSILON);
  BOOST_CHECK(std::fabs(Dot(p.GetOrigin()-v2,p.GetNormal()))<EPSILON);
  BOOST_CHECK(std::fabs(Dot(p.GetOrigin()-v3,p.GetNormal()))<EPSILON);
  
  // initialize a,b,c,d (ax+by+cz+d=0)
  Real a=UniformRandom();
  Real b=UniformRandom();
  Real c=UniformRandom();
  Real d=UniformRandom();
  p=Plane(a,b,c,d);
  BOOST_CHECK(std::fabs(d/Length(Vec3(a,b,c))+Dot(p.GetNormal(),p.GetOrigin()))<EPSILON);
  BOOST_CHECK(Length(Normalize(Vec3(a,b,c))-p.GetNormal())<EPSILON);
  
  
  //initialize based x-,y-,z- intercepts
  p=Plane(a,b,c);
  Plane p2(Vec3(a,0.0,0.0),Vec3(0.0,b,0.0),Vec3(0.0,0.0,c));
  BOOST_CHECK(Length(p.GetOrigin()-p2.GetOrigin())<EPSILON);
  BOOST_CHECK(Length(p.GetNormal()-p2.GetNormal())<EPSILON);
  
}

BOOST_AUTO_TEST_CASE(operators_composite3)
{
  Vec3 v1(2,3,4);
  Vec3 v2(4,3,2);
  Line3 l(v1,v2);
  Line3 l2(v1,v2);
  BOOST_CHECK(l==l2);
  
  Plane p(v1,v2);
  Plane p2(v1,v2);
  BOOST_CHECK(p==p2);
  
}

BOOST_AUTO_TEST_CASE(func_composite3)
{
  Vec3 v1(2,2,0);
  Vec3 v2(2,2,2);
  Vec3 v3(3,3,3);
  Vec3 v4(1,3,3);

  Plane p1(v1,v2);  
  Plane p2(v3,v4);
  Plane p3(v1,v2,v3);
  Plane p4(v1,v2,v4);
  Line3 l1(v1,v2);
  Line3 l2(v3,v4);
  Line3 l3(v1,v3);
  Line3 l4(v1+Vec3(2.0,2.0,2.0),v3+Vec3(2.0,2.0,2.0));
  
  Vec3 tmp1=v2-v1;
  Vec3 tmp2=v4-v3;
  Real ang=std::fabs(acos(Dot(tmp1,tmp2)/Length(tmp1)/Length(tmp2)));
  
  // angle between two lines
  BOOST_CHECK(dequal(Angle(l1,l2), ang));

  // angle between two planes
  BOOST_CHECK(dequal(std::fabs(acos(Dot(Normalize(v2),Normalize(v4)))),Angle(p1,p2)));
  
  // angle between line and plane
  BOOST_CHECK(dequal(M_PI/2.0-std::fabs(acos(Dot(Normalize(v2),
                                                 Normalize(tmp2)))),
                                              Angle(p1,l2)));

  BOOST_CHECK(Equal(IntersectionPoint(l2,l3), v3));
  
  BOOST_CHECK(!AreIntersecting(l3, l4));
  try{
    BOOST_CHECK_THROW(IntersectionPoint(l3,l4), GeomException);
  }catch(...)
  {
    BOOST_ERROR("Failed to catch GeomException");
  }
  BOOST_CHECK(Equal(IntersectionPoint(p1,l3), v1));
  BOOST_CHECK(Equal(IntersectionPoint(l3,p1), v1));
  
  BOOST_CHECK(EqualPosition(IntersectionLine(p3,p4),l1));

  Real d=Length(Dot(v3-v2,l1.GetDirection())*l1.GetDirection()-(v3-v2));
  BOOST_CHECK(dequal(Distance(l1,v3),d));
  
  BOOST_CHECK(dequal(Distance(p1,v3),Dot(p1.GetNormal(),v3)+p1.GetP()));


  BOOST_CHECK(!Equal(l1,Line3(l1.GetOrigin()+2.0*l1.GetDirection(),
                              l1.GetOrigin()+3.0*l1.GetDirection())));
  BOOST_CHECK(EqualPosition(l1,Line3(l1.GetOrigin()+2.0*l1.GetDirection(),
                                     l1.GetOrigin()-3.0*l1.GetDirection())));

  BOOST_CHECK(EqualPosition(l1,Line3(l1.GetOrigin()-2.0*l1.GetDirection(),
                                    l1.GetOrigin()+3.0*l1.GetDirection())));
  BOOST_CHECK(Equal(p3,Plane(v1,p3.GetNormal())));
  BOOST_CHECK(Equal(p3,Plane(v2,p3.GetNormal())));

  BOOST_CHECK(EqualPosition(p3,Plane(p3.GetOrigin(),-p3.GetNormal())));
  
  BOOST_CHECK(AreParallel(l1,l2)==false);
  BOOST_CHECK(AreParallel(l1,Line3(v1+Vec3(2.0,0,0),v2+Vec3(2.0,0,0)))==true);
  BOOST_CHECK(AreParallel(p1,p2)==false);
  BOOST_CHECK(AreParallel(p1, Plane(v3, p1.GetNormal()))==true);


  BOOST_CHECK(IsOnLine(l1,v1));
  BOOST_CHECK(IsOnLine(l1,v2));
  BOOST_CHECK(IsInPlane(p3,v1));
  BOOST_CHECK(!IsInPlane(p3,v4));
  BOOST_CHECK(IsInPlane(p3,l1));

  BOOST_CHECK(AreIntersecting(l1, l3));
  BOOST_CHECK(!AreIntersecting(l1, l2));

  BOOST_CHECK(RelativePosition(p3, v1)==0);
  BOOST_CHECK(RelativePosition(p3, v4)==-1);

}

BOOST_AUTO_TEST_CASE(rotation3)
{
  Vec3 v(1,0,0);
  Rotation3 r(Vec3(0,1,0), 30.0*M_PI/180.0);
  Vec3 vrot=r.Apply(v);
  BOOST_CHECK_CLOSE(Real(cos(30.0*M_PI/180.0)), Real(vrot[0]), Real(1e-5));
  BOOST_CHECK_SMALL(vrot[1],Real(1e-5));
  BOOST_CHECK_CLOSE(Real(-sin(30.0*M_PI/180.0)), Real(vrot[2]), Real(1e-5));

  Rotation3 r2(Vec3(0,1,0), 30.0*M_PI/180.0,Vec3(2,0,0));
  Vec3 vrot2=r2.Apply(v);
  BOOST_CHECK_CLOSE(Real(2.0)-Real(cos(30.0*M_PI/180.0)), Real(vrot2[0]), Real(1e-5));
  BOOST_CHECK_SMALL(vrot2[1],Real(1e-5));
  BOOST_CHECK_CLOSE(Real(sin(30.0*M_PI/180.0)), Real(vrot2[2]), Real(1e-5));
}

BOOST_AUTO_TEST_SUITE_END();
