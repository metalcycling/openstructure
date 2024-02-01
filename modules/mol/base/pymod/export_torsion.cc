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
#include <ost/mol/torsion_handle.hh>

using namespace boost::python;
using namespace ost;
using namespace ost::mol;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(X_rotate_angle_overloads, 
                                       TorsionHandle::RotateAngle, 1, 2)
                                       
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(X_set_angle_overloads, 
                                      TorsionHandle::SetAngle, 1, 2)                                       
void export_Torsion()
{
  class_<TorsionHandle>("TorsionHandle", init<>())
    .def("IsValid", &TorsionHandle::IsValid)
    .def("GetName", &TorsionHandle::GetName,
         return_value_policy<copy_const_reference>())
    .def("Name", 
         make_function(&TorsionHandle::GetName,
                       return_value_policy<copy_const_reference>()))
    .def("GetPos", &TorsionHandle::GetPos)
    .add_property("pos", &TorsionHandle::GetPos)
    .def("GetFirst", &TorsionHandle::GetFirst)
    .def("GetSecond", &TorsionHandle::GetSecond)
    .def("GetThird", &TorsionHandle::GetThird)
    .def("GetFourth", &TorsionHandle::GetFourth)
    .add_property("first", &TorsionHandle::GetFirst)
    .add_property("second", &TorsionHandle::GetSecond)
    .add_property("third", &TorsionHandle::GetThird)
    .add_property("fourth", &TorsionHandle::GetFourth)    
    .def("IsValid", &TorsionHandle::IsValid)
    .def("GetAngle", &TorsionHandle::GetAngle)    
    .add_property("angle", &TorsionHandle::GetAngle)
    .def(self_ns::str(self))

  ;
  class_<TorsionHandleList>("TorsionHandleList", no_init)
    .def(vector_indexing_suite<TorsionHandleList>())  
  ;
  
  def ("DihedralAngle", &DihedralAngle);
}
