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
  


Real (*Line2Angle)(const Line2& l1, const Line2& l2)                        = &Angle;
Vec2   (*Line2IntersectionPoint)(const Line2& l1, const Line2& l2)            = &IntersectionPoint;
Real (*Line2Distance)(const Line2& l, const Vec2& v)                        = &Distance;
bool   (*Line2IsOnLine)(const Line2& l, const Vec2& v, Real epsilon)        = &IsOnLine;
bool   (*Line2Equal)(const Line2& l1, const Line2& l2,Real epsilon)         = &Equal;
bool   (*Line2EqualPosition)(const Line2& l1, const Line2& l2,Real epsilon) = &EqualPosition;
bool   (*Line2AreParallel)(const Line2& l1, const Line2& l2,Real epsilon)   = &AreParallel;


void export_Composite2_op()
{
  using namespace geom;
  
  def("Angle",Line2Angle);
  def("IntersectionPoint",Line2IntersectionPoint);
  def("Distance",Line2Distance);
  def("IsOnLine",Line2IsOnLine);
  def("Equal",Line2Equal);
  def("EqualPosition",Line2EqualPosition);
  def("AreParallel",Line2AreParallel);

  def("SimplifyPolygon",&SimplifyPolygon);
  def("IsInRectangle",&IsInRectangle);
  def("IsInCircle",&IsInCircle);
}
