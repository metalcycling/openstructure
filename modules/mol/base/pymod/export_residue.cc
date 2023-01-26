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
#include <ost/mol/chem_class.hh>
#include <ost/mol/chem_type.hh>
#include <ost/mol/mol.hh>
#include <ost/geom/export_helper/vector.hh>
using namespace ost;
using namespace ost::mol;

#include <ost/export_helper/generic_property_def.hh>
#include "bounds.hh"
namespace {
  
  
  String sec_structure_as_string(SecStructure& s)
  {
    return String(1, char(s));
  }

  boost::python::list get_alt_group_names(const ost::mol::ResidueHandle& res) {
    std::vector<String>  v_alt_group_names = res.GetAltAtomGroupNames();
    boost::python::list alt_group_names;
    for(std::vector<String>::iterator it = v_alt_group_names.begin();
        it != v_alt_group_names.end(); ++it) {
      alt_group_names.append(*it);
    }
    return alt_group_names;
  }

  typedef EntityView (ResidueHandle::*QueryMethod)(const Query&, uint) const;
  typedef EntityView (ResidueHandle::*StringMethod)(const String&, uint) const;
  QueryMethod select_query=&ResidueHandle::Select;
  StringMethod select_string=&ResidueHandle::Select;

  void set_sec_struct1(ResidueBase* b, const SecStructure& s) {b->SetSecStructure(s);}
  void set_sec_struct2(ResidueBase* b, char c) {b->SetSecStructure(SecStructure(c));}

  void set_chemclass(ResidueBase* b, object po)
  {
    extract<ChemClass> ex1(po);
    if(ex1.check()) {
      b->SetChemClass(ex1());
    }
    extract<char> ex2(po);
    if(ex2.check()) {
      b->SetChemClass(ChemClass(ex2()));
    }
    std::string st=extract<std::string>(po);
    if(st.empty()) {
      throw Error("expected non-empty string as chem class");
    }
    b->SetChemClass(ChemClass(st[0]));
  }
}

void export_Residue()
{
  class_<ChemClass>("ChemClass", init<char>(args("chem_class")))
    .def(self!=self)
    .def(self==self)
    .def(self_ns::str(self))
    .def("IsPeptideLinking", &ChemClass::IsPeptideLinking)
    .def("IsNucleotideLinking", &ChemClass::IsNucleotideLinking)
  ;
  implicitly_convertible<char, ChemClass>();
  
  object ct_class = class_<ChemType>("ChemType", init<char>(args("chem_type")))
    .def(self!=self)
    .def(self==self)
    .def(self_ns::str(self))
    .def("IsIon", &ChemType::IsIon)
    .def("IsNucleotide", &ChemType::IsNucleotide)
    .def("IsSaccharide", &ChemType::IsSaccharide)
    .def("IsAminoAcid", &ChemType::IsAminoAcid)
    .def("IsCoenzyme", &ChemType::IsCoenzyme)
    .def("IsDrug", &ChemType::IsDrug)
    .def("IsNonCanonical", &ChemType::IsNonCanonical)
    .def("IsKnown", &ChemType::IsKnown)
    .def("IsWater", &ChemType::IsWater)
  ;
  implicitly_convertible<char, ChemType>();
  ct_class.attr("IONS")=char(ChemType::IONS);
  ct_class.attr("NONCANONICALMOLS")=char(ChemType::NONCANONICALMOLS);
  ct_class.attr("SACCHARIDES")=char(ChemType::SACCHARIDES);
  ct_class.attr("NUCLEOTIDES")=char(ChemType::NUCLEOTIDES);
  ct_class.attr("AMINOACIDS")=char(ChemType::AMINOACIDS);
  ct_class.attr("COENZYMES")=char(ChemType::COENZYMES);
  ct_class.attr("WATERCOORDIONS")=char(ChemType::WATERCOORDIONS);
  ct_class.attr("DRUGS")=char(ChemType::DRUGS);
  ct_class.attr("WATERS")=char(ChemType::WATERS);
  ct_class.attr("UNKNOWN")=char(ChemType::UNKNOWN);

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
    .def(self+=self)
    .def(self-=self)
    .def(self+self)
    .def(self-self)
    .def(self+=int())
    .def(self-=int())
    .def(self+int())
    .def(self-int())    
  ;
  class_<ResNumList>("ResNumList", init<>())
    .def(vector_indexing_suite<ResNumList>())
    .def(geom::VectorAdditions<ResNumList>())
  ;
  implicitly_convertible<int, ResNum>();
  
  scope().attr("PEPTIDE_LINKING")=char(ChemClass::PEPTIDE_LINKING);
  scope().attr("D_PEPTIDE_LINKING")=char(ChemClass::D_PEPTIDE_LINKING);
  scope().attr("L_PEPTIDE_LINKING")=char(ChemClass::L_PEPTIDE_LINKING);
  scope().attr("RNA_LINKING")=char(ChemClass::RNA_LINKING);
  scope().attr("DNA_LINKING")=char(ChemClass::DNA_LINKING);
  scope().attr("NON_POLYMER")=char(ChemClass::NON_POLYMER);
  scope().attr("L_SACCHARIDE")=char(ChemClass::L_SACCHARIDE);
  scope().attr("D_SACCHARIDE")=char(ChemClass::D_SACCHARIDE);
  scope().attr("SACCHARIDE")=char(ChemClass::SACCHARIDE);
  scope().attr("WATER")=char(ChemClass::WATER);
  scope().attr("UNKNOWN")=char(ChemClass::UNKNOWN);
  

  {
    scope sec_struct_scope=class_<SecStructure>("SecStructure", init<>())
      .def(init<char>())
      .def(init<SecStructure::Type>())
      .def("IsHelical", &SecStructure::IsHelical)
      .def("IsExtended", &SecStructure::IsExtended)
      .def("IsCoil", &SecStructure::IsCoil)
      .def("__str__", &sec_structure_as_string)
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
    .add_property("phi_torsion", &ResidueBase::GetPhiTorsion)
    .add_property("psi_torsion", &ResidueBase::GetPsiTorsion)
    .add_property("omega_torsion", &ResidueBase::GetOmegaTorsion)
    .def(self_ns::str(self))
    .def("GetOneLetterCode", &ResidueBase::GetOneLetterCode)
    .def("SetOneLetterCode", &ResidueBase::SetOneLetterCode)
    .add_property("one_letter_code", &ResidueBase::GetOneLetterCode, 
                 &ResidueBase::SetOneLetterCode)  
    .def("GetQualifiedName", &ResidueBase::GetQualifiedName)
    .add_property("qualified_name", &ResidueBase::GetQualifiedName)
    .def("IsPeptideLinking", &ResidueBase::IsPeptideLinking)
    .add_property("peptide_linking", &ResidueBase::IsPeptideLinking)
    .def("IsNucleotideLinking", &ResidueBase::IsNucleotideLinking)
    .add_property("nucleotide_linking", &ResidueBase::IsNucleotideLinking)
    
    .def("GetCentralAtom", &ResidueBase::GetCentralAtom)
    .def("SetCentralAtom", &ResidueBase::SetCentralAtom)
    .add_property("central_atom", &ResidueBase::GetCentralAtom, &ResidueBase::SetCentralAtom)
    .add_property("central_normal", &ResidueBase::GetCentralNormal)

    .def("GetKey", &ResidueBase::GetKey,
         return_value_policy<copy_const_reference>())
    .def("GetName", &ResidueBase::GetName,
         return_value_policy<copy_const_reference>())
    .def("GetNumber", &ResidueBase::GetNumber,
         return_value_policy<copy_const_reference>())
    .def("GetChemClass", &ResidueBase::GetChemClass)
    .add_property("chem_class", &ResidueBase::GetChemClass, set_chemclass)
    .def("SetChemClass", set_chemclass)
    .def("GetChemType", &ResidueBase::GetChemType)
    .add_property("chem_type", &ResidueBase::GetChemType)
    .add_property("is_ligand", &ResidueBase::IsLigand, &ResidueBase::SetIsLigand)
    .def("IsLigand", &ResidueBase::IsLigand)
    .def("SetIsLigand", &ResidueBase::SetIsLigand)
    .add_property("is_protein", &ResidueBase::IsProtein,
                                &ResidueBase::SetIsProtein)
    .def("IsProtein", &ResidueBase::IsProtein)
    .def("SetIsProtein", &ResidueBase::SetIsProtein)
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
    .def("IsValid", &ResidueBase::IsValid)
    .add_property("valid", &ResidueBase::IsValid)
  ;
  generic_prop_def<ResidueBase>(residue_base);

  class_<ResidueHandle, bases<ResidueBase> >("ResidueHandle", init<>())
    .def("GetChain",&ResidueHandle::GetChain)
    .add_property("chain", &ResidueHandle::GetChain)
    .def("GetEntity",&ResidueHandle::GetEntity)
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
    .def("GetAltAtomGroupNames", &get_alt_group_names)
    .def("GetCurrentAltGroupName", &ResidueHandle::GetCurrentAltGroupName,
         return_value_policy<copy_const_reference>())
    .def("SwitchAtomPos", &ResidueHandle::SwitchAtomPos)
    .add_property("atoms", &ResidueHandle::GetAtomList)
    .def("FindAtom", &ResidueHandle::FindAtom, args("atom_name"))
    .def("FindTorsion", &ResidueHandle::FindTorsion)
    .def("GetAtomCount", &ResidueHandle::GetAtomCount)
    .def("GetBondCount", &ResidueHandle::GetBondCount)
    .add_property("atom_count", &ResidueHandle::GetAtomCount)
    .add_property("bond_count", &ResidueHandle::GetBondCount)
    .add_property("index", &ResidueHandle::GetIndex)
    .def("Select", select_string, arg("flags")=0)
    .def("Select", select_query, arg("flags")=0)    
    .def("GetMass", &ResidueHandle::GetMass)
    .def("GetCenterOfMass", &ResidueHandle::GetCenterOfMass)
    .def("GetCenterOfAtoms", &ResidueHandle::GetCenterOfAtoms)
    .def("GetGeometricCenter", geom_center<ResidueHandle>)
    .add_property("mass", &ResidueHandle::GetMass)
    .add_property("center_of_mass", &ResidueHandle::GetCenterOfMass)
    .add_property("center_of_atoms", &ResidueHandle::GetCenterOfAtoms)  
    .add_property("geometric_center", geom_center<ResidueHandle>)
    .def("GetGeometricStart", geom_start<ResidueHandle>)
    .def("GetGeometricEnd", geom_end<ResidueHandle>)
    .def(self==self)
    .def(self!=self)
    .def("__hash__", &ResidueHandle::GetHashCode)
    .def("GetHashCode", &ResidueHandle::GetHashCode)
    .def("GetBounds", &ResidueHandle::GetBounds)
    .add_property("bounds", &ResidueHandle::GetBounds)
    .add_property("hash_code", &ResidueHandle::GetHashCode)    
  ;

  class_<ResidueHandleList>("ResidueHandleList", no_init)
    .def(vector_indexing_suite<ResidueHandleList>())
    .def(geom::VectorAdditions<ResidueHandleList>())    
  ;
  
  def("InSequence", &InSequence);
}
