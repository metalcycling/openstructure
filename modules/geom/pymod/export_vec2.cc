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
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>


#include <ost/geom/geom.hh>

using namespace boost::python;

const Real Vec2_getitem(const geom::Vec2& v, int i) {
  return v.At(i);
}
void Vec2_setitem(geom::Vec2& v,const  int i,const  Real val) {
  v.At(i)=val;
}


String vec2_repr(const geom::Vec2& v)
{
  std::stringstream ss;
  ss << "geom.Vec2(" << v[0] << ", " << v[1] << ")";
  return ss.str();
}

list vec2_data(const geom::Vec2& v)
{
  list nrvo;
  for(size_t k=0;k<2;++k) {
    nrvo.append(v.Data()[k]);
  }
  return nrvo;
}

void export_Vec2()
{
  using namespace geom;

  class_<Vec2>("Vec2",init<>())
    .def(init<Real,Real>())
    .def(init<const Vec2&>())
    .def(init<const Vec3&>())
    .def(init<const Vec4&>())
    .def(self *= Real())
    .def(self /= Real())
    .def(self += Real())
    .def(self += self)
    .def(self -= self)
    .def(-self)
    .def(self * Real())
    .def(self + Real())
    .def(Real() + self)
    .def(Real() * self)
    .def(self * Mat2())
    .def(self / Real())
    .def(self + self)
    .def(self - self)
    .def(self == self)
    .def(self != self)
    .def("__repr__", vec2_repr)
    .def(self_ns::str(self))
    .def("__getitem__",Vec2_getitem)
    .def("__setitem__",Vec2_setitem)
    .def("GetX", &Vec2::GetX)
    .def("GetY", &Vec2::GetY)
    .add_property("x", &Vec2::GetX, &Vec2::SetX)
    .add_property("y", &Vec2::GetY, &Vec2::SetY)
    .add_property("data",vec2_data)
  ;
  class_<Vec2List>("Vec2List", init<>())
    .def(vector_indexing_suite<Vec2List>())
  ;  
}
