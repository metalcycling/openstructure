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
#include <boost/python/slice.hpp>
using namespace boost::python;

#include <ost/geom/geom.hh>
using namespace geom;



const Real Mat3_getitem(const geom::Mat3& m, tuple i) {
  int a = extract<int> (i[0]);
  int b = extract<int> (i[1]);
  return m.At(a, b);
}
void Mat3_setitem(geom::Mat3& m,const  tuple i,const  Real val) { 
  int a = extract<int> (i[0]);
  int b = extract<int> (i[1]);
  m.At(a, b) = val;
}

Mat2 Mat3_getslice(const geom::Mat3& m, slice s) {
  tuple start=extract<tuple> (s.start());
  tuple end=extract<tuple> (s.stop());
  int start0=extract<int> (start[0]);
  int start1=extract<int> (start[1]);
  int end0=extract<int> (end[0]);
  int end1=extract<int> (end[1]);
  if(end0-start0!=1 || end1-start1!=1) throw GeomException("Invalid slice");
  return Mat2(m(start0,start1),m(start0,start1+1),m(start0+1,start1),m(start0+1,start1+1));
}
void Mat3_setslice(geom::Mat3& m,const  slice s,const  Mat2& m2)
{
  tuple start=extract<tuple> (s.start());
  tuple end=extract<tuple> (s.stop());
  int start0=extract<int> (start[0]);
  int start1=extract<int> (start[1]);
  int end0=extract<int> (end[0]);
  int end1=extract<int> (end[1]);
  if(end0-start0!=1 || end1-start1!=1) throw GeomException("Invalid slice");
  m(start0,start1)=m2(0,0);
  m(start0,start1+1)=m2(0,1);
  m(start0+1,start1)=m2(1,0);
  m(start0+1,start1+1)=m2(1,1);
}

String mat3_repr(const geom::Mat3& m)
{
  std::stringstream ss;

  ss << "geom.Mat3(" << m(0,0) << ", " << m(0,1) << ", " << m(0,2) << ", "
     << m(1, 0) << ", " << m(1,1) << ", " << m(1, 2) << ", "
     << m(2, 0) << "," << m(2, 1) << ", " << m(2, 2) << ")";
  return ss.str();
}

list mat3_data(const geom::Mat3& m)
{
  list nrvo;
  for(size_t k=0;k<9;++k) {
    nrvo.append(m.Data()[k]);
  }
  return nrvo;
}

void export_Mat3()
{
  class_<Mat3>("Mat3",init<>())
    .def(init<Real,Real,Real,Real,Real,Real,Real,Real,Real>())
    .def(init<const Mat2&>())
    .def(init<Real,Real,Real>())
    .def(self += self)
    .def(self -= self)
    .def("__repr__", mat3_repr)
    .def(self + self)
    .def(self - self)
    .def(self *= Real())
    .def(self /= Real())
    .def(self * Real())
    .def(self * Vec3())
    .def(self * self)
    .def(self *= self)
    .def(self / Real())
    .def(self == self)
    .def(self != self)
    .def(self_ns::str(self))
    .def("__getitem__",Mat3_getitem)
    .def("__getitem__",Mat3_getslice)
    .def("__setitem__",Mat3_setitem)
    .def("__setitem__",Mat3_setslice)
    .def("GetCol", &Mat3::GetCol)
    .def("GetRow", &Mat3::GetRow)
    .add_property("data",mat3_data)
  ;
}
