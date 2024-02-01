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
  Author: Ansgar Philippsen
*/

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;

#include <ost/img/point_list.hh>
#include <ost/img/size.hh>

namespace ost { namespace img {


int point_len(const Point& p) {return 3;}
int point_getitem(const Point& p, int i) {return p.At(i);}
void point_setitem(Point& p, int i, int v) {p.At(i)=v;}

}}

void export_Point()
{
  using namespace ost::img;
  // point wrapper
  class_<Point>("Point", init<optional<int,int,int> >(args("x","y","z")))
    .def(init<const Size&>(args("size")))
    .def(init<const Vec2&>(args("vec")))
    .def(init<const Vec3&>(args("vec")))
    .def(init<const Vec4&>(args("vec")))
    .def(init<const Point&>(args("point")))
    .def_readwrite("x", &Point::x)
    .def_readwrite("y", &Point::y)
    .def_readwrite("z", &Point::z)
    .def("Mirror",&Point::Mirror)
    .def("ToVec2",&Point::ToVec2)
    .def("ToVec3",&Point::ToVec3)
    .def("ToVec4",&Point::ToVec4)
    .def("__len__",point_len)
    .def("__getitem__",point_getitem)
    .def("__setitem__",point_setitem)
    .def(self += self)
    .def(self + self)
    .def(self -= self)
    .def(self - self)
    .def(self == self)
    .def(self != self)
    .def(self < self)
    .def(self <= self)
    .def(self > self)
    .def(self >= self)
    .def(- self)
    .def(self_ns::str(self))
    ;
}

void export_PointList()
{
  using namespace ost::img;
  class_<PointList>("PointList")
    .def(init<PeakList>())
    .def(vector_indexing_suite<PointList>())
    ;
}
