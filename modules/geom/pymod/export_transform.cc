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

#include <ost/geom/transform.hh>

using namespace boost::python;
using namespace geom;

void export_Transform()
{
  class_<Transform>("Transform", init<>())
    .def(init<const Transform&>()) // shouldn't this be there automatically ?
    .def("GetMatrix",&Transform::GetMatrix)
    .def("SetMatrix",&Transform::SetMatrix)
    .add_property("matrix",&Transform::GetMatrix,&Transform::SetMatrix)
    .def("GetTransposedMatrix",&Transform::GetTransposedMatrix)
    .add_property("tmatrix",&Transform::GetTransposedMatrix)
    .add_property("transposed_matrix",&Transform::GetTransposedMatrix)
    .def("GetInvertedMatrix",&Transform::GetInvertedMatrix)
    .add_property("inverted_matrix",&Transform::GetInvertedMatrix)
    .def("SetTrans",&Transform::SetTrans)
    .def("GetTrans",&Transform::GetTrans)
    .add_property("trans",&Transform::GetTrans,&Transform::SetTrans)
    .def("SetCenter",&Transform::SetCenter)
    .def("GetCenter",&Transform::GetCenter)
    .add_property("center",&Transform::GetCenter,&Transform::SetCenter)
    .def("SetRot",&Transform::SetRot)
    .def("GetRot",&Transform::GetRot)
    .add_property("rot",&Transform::GetRot,&Transform::SetRot)
    .def("ApplyXAxisRotation",&Transform::ApplyXAxisRotation)
    .def("ApplyYAxisRotation",&Transform::ApplyYAxisRotation)
    .def("ApplyZAxisRotation",&Transform::ApplyZAxisRotation)
    .def("ApplyXAxisTranslation",&Transform::ApplyXAxisTranslation)
    .def("ApplyYAxisTranslation",&Transform::ApplyYAxisTranslation)
    .def("ApplyZAxisTranslation",&Transform::ApplyZAxisTranslation)
    .def("ApplyAxisRotation",&Transform::ApplyAxisRotation)
    ;
}
