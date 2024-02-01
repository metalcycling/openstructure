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



const Real Vec4_getitem(const geom::Vec4& v, int i) {
  return v.At(i);
}

void Vec4_setitem(geom::Vec4& v,const  int i,const  Real val) {
  v.At(i)=val;
}

String vec4_repr(const geom::Vec4& v)
{
  std::stringstream ss;
  ss << "geom.Vec4(" << v[0] << ", " << v[1] << "," << v[2]
     << ", " << v[3] << ")";
  return ss.str();
}

list vec4_data(const geom::Vec4& v)
{
  list nrvo;
  for(size_t k=0;k<4;++k) {
    nrvo.append(v.Data()[k]);
  }
  return nrvo;
}

void export_Vec4()
{
  using namespace geom;

  class_<Vec4>("Vec4",init<>())
    .def(init<Real,Real,Real,Real>())
    .def(init<const Vec2&>())
    .def(init<const Vec3&>())
    .def(init<const Vec4&>())
    .def(self *= Real())
    .def(self /= Real())
    .def(self += Real())
    .def(self += self)
    .def("__repr__", vec4_repr)
    .def(self -= self)
    .def(-self)
    .def(self * Real())
    .def(self + Real())
    .def(Real() + self)
    .def(Real() * self)
    .def(self * Mat4())
    .def(self / Real())
    .def(self + self)
    .def(self - self)
    .def(self == self)
    .def(self != self)
    .def(self_ns::str(self))
    .def("__getitem__",Vec4_getitem)
    .def("__setitem__",Vec4_setitem)
    .add_property("x", &Vec4::GetX, &Vec4::SetX)
    .add_property("y", &Vec4::GetY, &Vec4::SetY)
    .add_property("z", &Vec4::GetZ, &Vec4::SetZ)
    .add_property("w", &Vec4::GetW, &Vec4::SetW)
    .add_property("data",vec4_data)
  ;

}
