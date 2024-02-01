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


#include <boost/python.hpp>
using namespace boost::python;

#include <ost/geom/geom.hh>
using namespace geom;
  


Real (*Line3Angle)(const Line3& l1, const Line3& l2)                            = &Angle;
Vec3   (*Line3IntersectionPoint)(const Line3& l1, const Line3& l2)                = &IntersectionPoint;
Real (*Line3Distance)(const Line3& l, const Vec3& v)                            = &Distance;
bool   (*Line3IsOnLine)(const Line3& l, const Vec3& v, Real epsilon)            = &IsOnLine;
bool   (*Line3Equal)(const Line3& l1, const Line3& l2,Real epsilon)             = &Equal;
bool   (*Line3EqualPosition)(const Line3& l1, const Line3& l2,Real epsilon)     = &EqualPosition;
bool   (*Line3AreParallel)(const Line3& l1, const Line3& l2,Real epsilon)       = &AreParallel;
Real (*PlanePlaneAngle)(const Plane& p1, const Plane& p2)                       = &Angle;
Real (*PlaneAngle)(const Plane& p, const Line3& l)                              = &Angle;
Real (*Vec3Angle)(const Vec3& p, const Vec3& l)                                 = &Angle;
bool   (*PlaneVec3IsInPlane)(const Plane& p,  const Vec3& v,Real epsilon)       = &IsInPlane;
bool   (*PlaneLine3IsInPlane)(const Plane& p,  const Line3& l,Real epsilon)     = &IsInPlane;
int    (*PlaneRelativePosition)(const Plane& p,  const Vec3& v,Real epsilon)    = &RelativePosition;
bool   (*PlaneEqual)(const Plane& p1, const Plane& p2,Real epsilon)             = &Equal;
bool   (*PlaneEqualPosition)(const Plane& p1, const Plane& p2,Real epsilon)     = &EqualPosition;
bool   (*PlaneAreParallel)(const Plane& p,  const Line3& l,Real epsilon)        = &AreParallel;
bool   (*PlanePlaneAreParallel)(const Plane& p1, const Plane& p2,Real epsilon)  = &AreParallel;
Real (*PlaneDistance)(const Plane& p, const Vec3& v)                            = &Distance;
Vec3   (*Line3PlaneIntersectionPoint)(const Line3& l, const Plane& p)             = &IntersectionPoint;
Vec3   (*PlaneIntersectionPoint)(const Plane& p, const Line3& l)                  = &IntersectionPoint;
Line3  (*Line3IntersectionLine)(const Plane& p1, const Plane& p2)                 = &IntersectionLine;
Real (*Line3Line3Distance)(const Line3& l1, const Line3& l2)                    = &Distance;
bool   (*Line3AreIntersecting)( const Line3& l1, const Line3& l2, Real epsilon) = &AreIntersecting;


void export_Composite3_op()
{
  using namespace geom;
  
  def("Angle",Line3Angle);
  def("IntersectionPoint",Line3IntersectionPoint);
  def("Distance",Line3Distance);
  def("IsOnLine",Line3IsOnLine);
  def("Equal",Line3Equal);
  def("EqualPosition",Line3EqualPosition);
  def("AreParallel",Line3AreParallel);
  def("Angle",PlanePlaneAngle);
  def("Angle",PlaneAngle);
  def("Angle",Vec3Angle);
  def("IsInPlane",PlaneVec3IsInPlane);
  def("IsInPlane",PlaneLine3IsInPlane);
  def("PlaneRelativePosition",RelativePosition);
  def("Equal",PlaneEqual);
  def("EqualPosition",PlaneEqualPosition);
  def("AreParallel",PlaneAreParallel);
  def("AreParallel",PlanePlaneAreParallel);
  def("Distance",PlaneDistance);
  def("IntersectionPoint",Line3PlaneIntersectionPoint);
  def("IntersectionPoint",PlaneIntersectionPoint);
  def("IntersectionLine",Line3IntersectionLine);
  def("Distance",Line3Line3Distance);
  def("AreIntersecting",Line3AreIntersecting);
  def("IsInSphere",&IsInSphere);
}
