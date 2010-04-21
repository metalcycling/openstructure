//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2010 by the OpenStructure authors
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
#include "generic_property_def.hh"

#include "generic_property_def.hh"

namespace {
  typedef EntityView (ResidueHandle::*QueryMethod)(const Query&, uint) const;
  typedef EntityView (ResidueHandle::*StringMethod)(const String&, uint) const;
  QueryMethod select_query=&ResidueHandle::Select;
  StringMethod select_string=&ResidueHandle::Select;

  void set_sec_struct1(ResidueBase* b, const SecStructure& s) {b->SetSecStructure(s);}
  void set_sec_struct2(ResidueBase* b, char c) {b->SetSecStructure(SecStructure(c));}

}

//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(X_insert_overloads,
//                                       ResidueHandle::InsertAtom, 2, 3)
void export_Residue()
{
  class_<ResNum>("ResNum", init<int>(args("num")))
    .def(init<int,char>(args("num", "ins_code")))
    .def("GetNum", &ResNum::GetNum)
    .def("GetInsCode", &ResNum::GetInsCode)
    .add_property("num", &ResNum::GetNum)
    .add_property("ins_code", &ResNum::GetInsCode)
    .def("__str__", &ResNum::AsString)
    .def("__repr__", &ResNum::AsString)
    .def(self<self)
    .def(self>self)
    .def(self>=self)
    .def(self<=self)
    .def(self==self)
    .def(self!=self)    
    .def(self+=int())
    .def(self-=int())
    .def(self+int())
    .def(self-int())    
  ;
  {
    scope sec_struct_scope=class_<SecStructure>("SecStructure", init<>())
      .def(init<char>())
      .def(init<SecStructure::Type>())
      .def("IsHelical", &SecStructure::IsHelical)
      .def("IsExtended", &SecStructure::IsExtended)
      .def("IsCoil", &SecStructure::IsCoil)        
    ;
    enum_<SecStructure::Type>("Type")
      .value("EXTENDED", SecStructure::EXTENDED)
      .value("COIL", SecStructure::COIL)
      .value("THREE_TEN_HELIX", SecStructure::THREE_TEN_HELIX)
      .value("ALPHA_HELIX", SecStructure::ALPHA_HELIX)
      .value("BETA_BRIDGE", SecStructure::BETA_BRIDGE)
      .value("BEND", SecStructure::BEND)
      .value("PI_HELIX", SecStructure::PI_HELIX)
      .value("TURN", SecStructure::TURN)
      .export_values()
    ;
  }
  class_<ResidueBase> residue_base("ResidueBase", no_init);
  residue_base
    .def("GetSecStructure", &ResidueBase::GetSecStructure)
    .def("SetSecStructure", set_sec_struct1)
    .def("SetSecStructure", set_sec_struct2)
    .def("GetPhiTorsion", &ResidueBase::GetPhiTorsion)
    .def("GetPsiTorsion", &ResidueBase::GetPsiTorsion)
    .def("GetOmegaTorsion", &ResidueBase::GetOmegaTorsion)
    .def("IsValid", &ResidueBase::IsValid)
    .def(self_ns::str(self))
    .def("GetOneLetterCode", &ResidueBase::GetOneLetterCode)
    .def("SetOneLetterCode", &ResidueBase::SetOneLetterCode)
    .add_property("one_letter_code", &ResidueBase::GetOneLetterCode, 
                 &ResidueBase::SetOneLetterCode)  
    .def("GetQualifedName", &ResidueBase::GetQualifiedName)
    .def("IsPeptideLinking", &ResidueBase::IsPeptideLinking)
    .def("GetKey", &ResidueBase::GetKey,
         return_value_policy<copy_const_reference>())
     .def("GetName", &ResidueBase::GetName,
         return_value_policy<copy_const_reference>())
    .def("GetNumber", &ResidueBase::GetNumber,
         return_value_policy<copy_const_reference>())
    .add_property("number",
                   make_function(&ResidueBase::GetNumber,
                                 return_value_policy<copy_const_reference>()))
    .add_property("sec_structure", &ResidueBase::GetSecStructure,
                  &ResidueBase::SetSecStructure)
    .add_property("key",
                   make_function(&ResidueBase::GetKey,
                                 return_value_policy<copy_const_reference>()))
    .add_property("name",
                   make_function(&ResidueBase::GetName,
                                 return_value_policy<copy_const_reference>()))
    .add_property("qualified_name", &ResidueBase::GetQualifiedName)
  ;
  generic_prop_def(residue_base);

  generic_prop_def(residue_base);

  class_<ResidueHandle, bases<ResidueBase> >("ResidueHandle", init<>())
    .def("GetChain",&ResidueHandle::GetChain)
    .add_property("chain", &ResidueHandle::GetChain)
    .add_property("entity", &ResidueHandle::GetEntity)
    .def("GetAtomList", &ResidueHandle::GetAtomList)
    .def("GetIndex", &ResidueHandle::GetIndex)
    .def("GetNext", &ResidueHandle::GetNext)
    .def("GetPrev", &ResidueHandle::GetPrev)
    .def("GetTorsionList", &ResidueHandle::GetTorsionList)
    .add_property("prev", &ResidueHandle::GetPrev)
    .add_property("next", &ResidueHandle::GetNext)
    .def("HasAltAtoms", &ResidueHandle::HasAltAtoms)
    .def("GetHandle", &ResidueHandle::GetHandle)
    .add_property("handle", &ResidueHandle::GetHandle)    
    .def("HasAltAtomGroup", &ResidueHandle::HasAltAtomGroup)
    .def("GetCurrentAltGroupName", &ResidueHandle::GetCurrentAltGroupName,
         return_value_policy<copy_const_reference>())
    .def("SwitchAtomPos", &ResidueHandle::SwitchAtomPos)
    .add_property("atoms", &ResidueHandle::GetAtomList)
    .def("FindAtom", &ResidueHandle::FindAtom, args("atom_name"))
    .def("FindTorsion", &ResidueHandle::FindTorsion)
    .def("Apply", &ResidueHandle::Apply, args("visitor"))
    .def("GetAtomCount", &ResidueHandle::GetAtomCount)
    .def("GetBondCount", &ResidueHandle::GetBondCount)
    .add_property("atom_count", &ResidueHandle::GetAtomCount)
    .add_property("index", &ResidueHandle::GetIndex)
    .def("Select", select_string, arg("flags")=0)
    .def("Select", select_query, arg("flags")=0)
    .def("GetMass", &ResidueHandle::GetMass)
    .def("GetCenterOfMass", &ResidueHandle::GetCenterOfMass)
    .def("GetCenterOfAtoms", &ResidueHandle::GetCenterOfAtoms)
    .def("GetGeometricCenter", &ResidueHandle::GetGeometricCenter)
    .add_property("geometric_center", &ResidueHandle::GetGeometricCenter)
    .def("GetGeometricStart", &ResidueHandle::GetGeometricStart)
    .def("GetGeometricEnd", &ResidueHandle::GetGeometricEnd)
    .def(self==self)
    .def(self!=self)
    .def("__hash__", &ResidueHandle::GetHashCode)
  ;

  class_<ResidueHandleList>("ResidueHandleList", no_init)
    .def(vector_indexing_suite<ResidueHandleList>())
  ;
}
