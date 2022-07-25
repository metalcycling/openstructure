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
#include <ost/geom/vec3.hh>
#include <ost/geom/geom.hh>
#include <ost/geom/export_helper/vector.hh>

using namespace boost::python;


const Real Vec3_getitem(const geom::Vec3& v, int i) {
  return v.At(i);
}
void Vec3_setitem(geom::Vec3& v,const  int i,const  Real val) {
  v.At(i)=val;
}

geom::Vec3 NormalizeV3(const geom::Vec3& v) {
  return geom::Normalize(v);
}

String vec3_repr(const geom::Vec3& v)
{
  std::stringstream ss;
  ss << "geom.Vec3(" << v[0] << ", " << v[1] << "," << v[2] << ")";
  return ss.str();
}

list vec3_data(const geom::Vec3& v)
{
  list nrvo;
  for(size_t k=0;k<3;++k) {
    nrvo.append(v.Data()[k]);
  }
  return nrvo;
}

boost::python::tuple wrap_FitCylinder(const geom::Vec3List& vl,const geom::Vec3& initial_direction) {
  std::pair<geom::Line3,Real> pair=vl.FitCylinder(initial_direction);
  return boost::python::make_tuple<geom::Line3,Real>(pair.first,pair.second);
}  

void export_Vec3()
{
  using namespace geom;
  
  class_<Vec3>("Vec3",init<>())
    .def(init<Real,Real,Real>())
    .def(init<const Vec2&>())
    .def(init<const Vec4&>())
    .def(init<const Vec3&>())
    .def(self *= Real())
    .def(self /= Real())
    .def(self += Real())
    .def(self += self)
    .def(self -= self)
    .def(-self)
    .def(self * Real())
    .def(Real() * self)
    .def(self * Mat3())
    .def(self / Real())
    .def(self + self)
    .def(self + Real())
    .def(Real() + self)
    .def(self - self)
    .def(self == self)
    .def(self != self)
    .def(self_ns::str(self))
    .def("__getitem__",Vec3_getitem)
    .def("__setitem__",Vec3_setitem)
    .def("__repr__", vec3_repr)
    .def("GetX", &Vec3::GetX)
    .def("GetY", &Vec3::GetY)
    .def("GetZ", &Vec3::GetZ)
    .add_property("x", &Vec3::GetX, &Vec3::SetX)
    .add_property("y", &Vec3::GetY, &Vec3::SetY)
    .add_property("z", &Vec3::GetZ, &Vec3::SetZ)
    .add_property("data",vec3_data)
  ;
  
  def("Normalize", &NormalizeV3);
  def("Cross", &Cross);
  
  class_<Vec3List>("Vec3List", init<>())
    .def(vector_indexing_suite<Vec3List>())
    .def(geom::VectorAdditions<Vec3List>())
    .def(self *= Real())
    .def(self /= Real())
    .def(self += Real())
    .def(self += self)
    .def(self -= self)
    //.def(-self)
    .def(self * Real())
    .def(Real() * self)
    .def(self / Real())
    .def(self + self)
    .def(self + Real())
    .def(Real() + self)
    .def(self - self)
    .def(self == self)
    .def(self != self)
    .add_property("center", &Vec3List::GetCenter)
    .add_property("inertia", &Vec3List::GetInertia)
    .add_property("principal_axes", &Vec3List::GetPrincipalAxes)
    .def("GetODRLine", &Vec3List::GetODRLine)
    .def("FitCylinder", wrap_FitCylinder,(arg("direction initial guess")))
    .def("ApplyTransform", &Vec3List::ApplyTransform,(arg("transform")))
    .def("GetSummedSquaredDistances", &Vec3List::GetSummedSquaredDistances,(arg("other")))
    .def("GetRMSD", &Vec3List::GetRMSD,(arg("other")))
    .def("GetGDTHA", &Vec3List::GetGDTHA, (arg("other"), arg("norm")=true))
    .def("GetGDTTS", &Vec3List::GetGDTTS, (arg("other"), arg("norm")=true))
  ;
}
