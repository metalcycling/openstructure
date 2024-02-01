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

using namespace boost::python;

#include <ost/mol/mol.hh>

using namespace ost;
using namespace ost::mol;

void export_PropertyID()
{
  {    
    scope prop_scope=class_<Prop>("Prop", no_init)
      .def("GetName", &Prop::GetName)
      .def("GetTypeName", &Prop::GetTypeName)
      .def_readwrite("id", &Prop::id)
      .def_readwrite("level", &Prop::level)
      .def_readwrite("type", &Prop::type)    
      .add_property("name", &Prop::GetName)
      .add_property("type_name", &Prop::GetTypeName)
    ;
    enum_<Prop::ID>("ID")
      .value("UNDEF", Prop::UNDEF) 
      .value("RNAME", Prop::RNAME) 
      .value("ANAME", Prop::ANAME) 
      .value("CNAME", Prop::CNAME) 
      .value("ELE", Prop::ELE) 
      .value("RNUM", Prop::RNUM) 
      .value("ANUM", Prop::ANUM) 
      .value("AX", Prop::AX) 
      .value("AY", Prop::AY) 
      .value("AZ", Prop::AZ) 
      .value("OCC", Prop::OCC) 
      .value("CUSTOM", Prop::CUSTOM) 
      .value("RTYPE", Prop::RTYPE) 
      .value("ISHETATM", Prop::ISHETATM) 
      .value("RBFAC", Prop::RBFAC) 
      .value("ABFAC", Prop::ABFAC) 
      .value("PEPTIDE", Prop::PEPTIDE) 
      .value("ACHARGE", Prop::ACHARGE)
    ;
    enum_<Prop::Level>("Level")
      .value("ATOM", Prop::ATOM)
      .value("RESIDUE", Prop::RESIDUE)
      .value("CHAIN", Prop::CHAIN)
      .value("UNSPECIFIED", Prop::UNSPECIFIED)
    ;
  }
  
  def("PropertyIDFromString", &PropertyIDFromString);
  def("PropertyFromString", &PropertyFromString);  
}
