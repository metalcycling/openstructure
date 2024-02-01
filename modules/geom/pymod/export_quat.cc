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

geom::Quat normalize (const geom::Quat& q)
{
  return geom::Normalize(q);
}

void export_Quat()
{
  using namespace geom;

  class_<Quat>("Quat",init<>())
    .def(init<Real, Real,Real,Real>())
    .def(init<const Quat&>())
    .def(init<Real, const Vec3&>())
    .def(init<const Mat3&>())
    .def(self *= Real())
    .def(self /= Real())
    .def(self += self)
    .def(self -= self)
    .def(self == self)
    .def(self != self)
    .def(-self)
    .def(self * Real())
    .def(self * Quat())
    .def(self / Real())
    .def(self + self)
    .def(self - self)
    .def(self_ns::str(self))
    .def("ToRotationMatrix", &Quat::ToRotationMatrix)
    .def("Rotate", &Quat::Rotate)
    .def("GetAxis", &Quat::GetAxis)
    .def("GetAngle", &Quat::GetAngle)
    .def_readwrite("w",&Quat::w)
    .def_readwrite("x",&Quat::x)
    .def_readwrite("y",&Quat::y)
    .def_readwrite("z",&Quat::z)
  ;
  def("Conjugate",&Conjugate);
  def("Slerp",&Slerp);
  def("Normalize",normalize);
}

