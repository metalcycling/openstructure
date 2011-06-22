//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2011 by the OpenStructure authors
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

using namespace boost::python;


#include <ost/mol/mol.hh>

using namespace ost;
using namespace ost::mol;

#include <ost/export_helper/generic_property_def.hh>
#include <ost/export_helper/vector.hh>

void export_Bond()
{
  class_<BondHandle> bond_handle("BondHandle", init<>());
  bond_handle
    .def("GetPos", &BondHandle::GetPos)
    .add_property("pos", &BondHandle::GetPos)
    .add_property("first", 
                  make_function(&BondHandle::GetFirst))
    .add_property("second", 
                  make_function(&BondHandle::GetSecond))
    .add_property("length", 
                  &BondHandle::GetLength)
    .add_property("bond_order",
                  &BondHandle::GetBondOrder,
                  &BondHandle::SetBondOrder)
    .def("GetFirst", &BondHandle::GetFirst)    
    .def("GetSecond",&BondHandle::GetSecond)
    .def("GetOther",&BondHandle::GetOther)
    .def("GetLength",&BondHandle::GetLength)
    .def("GetBondOrder",&BondHandle::GetBondOrder)
    .def("SetBondOrder",&BondHandle::SetBondOrder)
    .def("IsValid", &BondHandle::IsValid)       
    .def(self == self)
    .def(self != self)
    .def(self_ns::str(self))
  ;    
  generic_prop_def<BondHandle>(bond_handle);
  class_<BondHandleList>("BondHandleList", no_init)
    .def(vector_indexing_suite<BondHandleList>())
    .def(ost::VectorAdditions<BondHandleList>())
  ;  
  def("BondExists", &BondExists);
}

