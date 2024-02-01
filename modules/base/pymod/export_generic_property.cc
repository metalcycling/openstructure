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

#include <ost/generic_property.hh>

using namespace ost;

void export_GenericProp()
{
  String (GenericPropContainer::* get_string1)(const String&) const = &GenericPropContainer::GetStringProp;
  String (GenericPropContainer::* get_string2)(const String&, const String&) const = &GenericPropContainer::GetStringProp;
  Real (GenericPropContainer::* get_float1)(const String&) const = &GenericPropContainer::GetFloatProp;
  Real (GenericPropContainer::* get_float2)(const String&, Real) const = &GenericPropContainer::GetFloatProp;
  int (GenericPropContainer::* get_int1)(const String&) const = &GenericPropContainer::GetIntProp;
  int (GenericPropContainer::* get_int2)(const String&, int) const = &GenericPropContainer::GetIntProp;
  bool (GenericPropContainer::* get_bool1)(const String&) const = &GenericPropContainer::GetBoolProp;
  bool (GenericPropContainer::* get_bool2)(const String&, bool) const = &GenericPropContainer::GetBoolProp;
  geom::Vec3 (GenericPropContainer::* get_vec3)(const String&) const = &GenericPropContainer::GetVec3Prop;

  class_<GenericPropContainer, boost::noncopyable>("GenericPropContainer",no_init)
    .def("HasProp",&GenericPropContainer::HasProp)
    .def("SetStringProp",&GenericPropContainer::SetStringProp)
    .def("GetStringProp",get_string1)
    .def("GetStringProp",get_string2)
    .def("SetFloatProp",&GenericPropContainer::SetFloatProp)
    .def("GetFloatProp",get_float1)
    .def("GetFloatProp",get_float2)
    .def("SetIntProp",&GenericPropContainer::SetIntProp)
    .def("GetIntProp",get_int1)
    .def("GetIntProp",get_int2)
    .def("SetBoolProp",&GenericPropContainer::SetBoolProp)
    .def("GetBoolProp",get_bool1)
    .def("GetBoolProp",get_bool2)
    .def("SetVec3Prop",&GenericPropContainer::SetVec3Prop)
    .def("GetVec3Prop",get_vec3)
    .def("ClearProps",&GenericPropContainer::ClearProps)
    .def("GetPropAsString",&GenericPropContainer::GetPropAsString)
    ;
}
