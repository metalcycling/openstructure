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
#include <iostream>
extern void export_Vec2();
extern void export_VecMat2_op();
extern void export_Vec3();
extern void export_VecMat3_op();
extern void export_Vec4();
extern void export_VecMat4_op();
extern void export_Mat2();
extern void export_Mat3();
extern void export_Mat4();
extern void export_Composite2();
extern void export_Composite3();
extern void export_Composite2_op();
extern void export_Composite3_op();
extern void export_Quat();
extern void export_Transform();

namespace details {
#if (defined(OST_STATIC_PROPERTY_WORKAROUND))
struct Axis_ {
  int X() { return geom::Axis::X; }
  int Y() { return geom::Axis::Y; }
  int Z() { return geom::Axis::Z; }    
};
#else
struct Axis_ {};
#endif
}

BOOST_PYTHON_MODULE(_ost_geom)
{
  export_Vec2();
  export_VecMat2_op();
  export_Vec3();
  export_VecMat3_op();
  export_Vec4();
  export_VecMat4_op();
  export_Mat2();
  export_Mat3();
  export_Mat4();
  export_Composite2();
  export_Composite3();
  export_Composite2_op();
  export_Composite3_op();
  export_Quat();
  export_Transform();

#if (defined(OST_STATIC_PROPERTY_WORKAROUND))
  // workaround for a problem with boost python and python 2.6.3/4
  class_<details::Axis_> axis("Axis_");
  axis.add_property("X", &details::Axis_::X);
  axis.add_property("Y", &details::Axis_::Y);
  axis.add_property("Z", &details::Axis_::Z);
 ;
 scope().attr("Axis")=axis();
#else
  class_<details::Axis_>("Axis")
    .def_readonly("X", geom::Axis::X)
    .def_readonly("Y", geom::Axis::Y)
    .def_readonly("Z", geom::Axis::Z)
  ;
#endif
}
