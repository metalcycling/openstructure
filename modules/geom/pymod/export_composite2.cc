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

geom::Vec2 Rectangle2_getitem(const geom::Rectangle2& r, unsigned int i)
{
  return r[i];
}

void Rectangle2_setitem(geom::Rectangle2& r,unsigned int i,
                        const  geom::Vec2& val)
{
  r[i]=val;
}

void export_Composite2()
{
  using namespace geom;

  class_<Line2>("Line2",init<>())
    .def(init<const Vec2&,const Vec2&>())
    .def(self == self)
    .def("At",&Line2::At)
    .def("GetOrigin",&Line2::GetOrigin)
    .def("GetDirection",&Line2::GetDirection)
    .add_property("direction", &Line2::GetDirection)
    .add_property("origin", &Line2::GetOrigin)
  ;

  class_<Polygon2>("Polygon2",init<>())
    .def(init<const std::vector<Vec2>& >())
    .def("GetNode",&Polygon2::GetNode)
    .def("AddNode",&Polygon2::AddNode)
    .def("GetNodeCount",&Polygon2::GetNodeCount)
    .def("GetArea",&Polygon2::GetArea)
    .def("GetCentroid",&Polygon2::GetCentroid)
    .def("Expand",&Polygon2::Expand)
    .def(self + Vec2())
    .def(self * Real())
  ;
  class_<Rectangle2>("Rectangle2",init<>())
    .def(init<Vec2,Vec2 >())
    .def("GetWidth",&Rectangle2::GetWidth)
    .def("GetHeight",&Rectangle2::GetHeight)
    .def("GetArea",&Rectangle2::GetArea)
    .def("GetStart",&Rectangle2::GetStart)
    .def("GetEnd",&Rectangle2::GetEnd)
    .def("SetStart",&Rectangle2::SetStart)
    .def("SetEnd",&Rectangle2::SetEnd)
    .def("__getitem__",Rectangle2_getitem)
    .def("__setitem__",Rectangle2_setitem)
  ;
  class_<Ellipse2>("Ellipse2",init<>())
   .def(init<Vec2,Real,Real,Real >())
   .def("At",&Ellipse2::At)
   .def("AtAngle",&Ellipse2::AtAngle)
   .def("GetBoundingBox",&Ellipse2::GetBoundingBox)
   .def("GetA",&Ellipse2::GetA)
   .def("GetB",&Ellipse2::GetB)
   .def("GetGamma",&Ellipse2::GetGamma)
   .def("GetArea",&Ellipse2::GetArea)
   .def("GetOrigin",&Ellipse2::GetOrigin)
   .def("SetA",&Ellipse2::SetA)
   .def("SetB",&Ellipse2::SetB)
   .def("SetGamma",&Ellipse2::SetGamma)
   .def("SetOrigin",&Ellipse2::SetOrigin)
 ;
  class_<Circle2>("Circle2",init<>())
   .def(init<Vec2&,Real>())
.def("GetRadius",&Circle2::GetRadius)
.def("SetRadius",&Circle2::SetRadius)
.def("GetCenter",&Circle2::GetCenter)
.def("SetCenter",&Circle2::SetCenter)
.def("GetCircumference",&Circle2::GetCircumference)
.def("GetArea",&Circle2::GetArea)
;

}
