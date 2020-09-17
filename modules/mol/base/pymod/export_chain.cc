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
#include <ost/mol/chain_type.hh>
#include <ost/geom/export_helper/vector.hh>
using namespace ost;
using namespace ost::mol;
#include <ost/export_helper/generic_property_def.hh>
#include "bounds.hh"

namespace {
  typedef void (ChainHandle::*RnumMethod)(const ResNum&);
  typedef void (ChainHandle::*RhandleMethod)(const ResidueHandle&);  
  
  
  typedef ResidueHandle (ChainHandle::*SingleArgMethod)(const ResidueKey&);  
  typedef ResidueHandle (ChainHandle::*DoubleArgMethod)(const ResidueKey&, const ResNum&);    
 
  typedef EntityView (ChainHandle::*QueryMethod)(const Query&, uint) const;
  typedef EntityView (ChainHandle::*StringMethod)(const String&, uint) const;  
  QueryMethod select_query=&ChainHandle::Select;
  StringMethod select_string=&ChainHandle::Select;
  ChainType (*ChainTypeFromStringPtr)(const String& identifier) =
    &ChainTypeFromString;
}

void export_Chain()
{
  class_<ChainBase> chain_base("ChainBase", no_init);
  chain_base
    .def("GetName", &ChainBase::GetName)
    .add_property("name", &ChainBase::GetName)
    .def(self_ns::str(self))
    .add_property("valid", &ChainBase::IsValid)
    .def("IsValid", &ChainBase::IsValid)
    .def("GetType", &ChainBase::GetType)
    .def("GetDescription", &ChainBase::GetDescription)
    .def("IsPolypeptide", &ChainBase::IsPolypeptide)
    .def("IsPolynucleotide", &ChainBase::IsPolynucleotide)
    .def("IsPolysaccharide", &ChainBase::IsPolysaccharide)
    .def("IsOligosaccharide", &ChainBase::IsOligosaccharide)
    .def("IsPolymer", &ChainBase::IsPolymer)
    .add_property("is_polypeptide", &ChainBase::IsPolypeptide)
    .add_property("is_polynucleotide", &ChainBase::IsPolynucleotide)
    .add_property("is_polysaccharide", &ChainBase::IsPolysaccharide)
    .add_property("is_oligosaccharide", &ChainBase::IsOligosaccharide)
    .add_property("is_polymer", &ChainBase::IsPolymer)
    .add_property("type", &ChainBase::GetType)
    .add_property("description", &ChainBase::GetDescription)
  ;
  generic_prop_def<ChainBase>(chain_base);
  class_<ChainHandle, bases<ChainBase> >("ChainHandle", init<>())
    .def("GetAtomList", &ChainHandle::GetAtomList)
    .add_property("atoms", &ChainHandle::GetAtomList)
    .def("GetResidueList", &ChainHandle::GetResidueList)
    .add_property("residues", &ChainHandle::GetResidueList)   
    .add_property("entity", &ChainHandle::GetEntity) 
    .def("GetEntity", &ChainHandle::GetEntity)
    .def("GetResidueByIndex", &ChainHandle::GetResidueByIndex)
    .def("GetHandle", &ChainHandle::GetHandle)
    .add_property("handle", &ChainHandle::GetHandle)    
    .def("FindResidue", &ChainHandle::FindResidue, arg("residue_number"))    
    .def("FindAtom", &ChainHandle::FindAtom, args("residue_number", "atom_name"))        
    //.def("AppendResidue", append_one_arg, args("residue_key"))
    //.def("AppendResidue", append_two_arg, args("residue_key", 
    //                                            "residue_number"))
    .def("GetNext", &ChainHandle::GetNext)
    .def("GetPrev", &ChainHandle::GetPrev)
    //.def("InsertResidueBefore", &ChainHandle::InsertResidueBefore, 
    //     args("residue_number", "key"))
    .def("GetPrev", &ChainHandle::GetPrev)
    //.def("DeleteResidue", delete_by_number, args("residue_number")) 
    //.def("DeleteResidue", delete_by_handle, args("residue_handle"))        
    .def("GetResidueCount", &ChainHandle::GetResidueCount)
    .def("GetAtomCount", &ChainHandle::GetAtomCount)
    .def("GetBondCount", &ChainHandle::GetBondCount)   
    .add_property("residue_count", &ChainHandle::GetResidueCount)
    .add_property("atom_count", &ChainHandle::GetAtomCount)
    .add_property("chain_type", &ChainHandle::GetType)
    .add_property("description", &ChainHandle::GetDescription)
    .def("InSequence", &ChainHandle::InSequence)
    .def("Select", select_string, arg("flags")=0)
    .def("Select", select_query, arg("flags")=0)
    .def("GetMass", &ChainHandle::GetMass)
    .def("GetCenterOfMass", &ChainHandle::GetCenterOfMass)
    .def("GetCenterOfAtoms", &ChainHandle::GetCenterOfAtoms)
    .def("GetGeometricCenter", geom_center<ChainHandle>)
    .add_property("geometric_center", geom_center<ChainHandle>)
    .add_property("mass", &ChainHandle::GetMass)
    .add_property("center_of_mass", &ChainHandle::GetCenterOfMass)
    .add_property("center_of_atoms", &ChainHandle::GetCenterOfAtoms)  
    .add_property("in_sequence", &ChainHandle::InSequence)  
    .def("GetBounds", &ChainHandle::GetBounds)
    .add_property("bounds", &ChainHandle::GetBounds)
    .def("GetGeometricStart", geom_start<ChainHandle>)
    .def("GetGeometricEnd", geom_end<ChainHandle>)
    .def(self==self)
    .def(self!=self)
  ;
  
  class_<ChainHandleList>("ChainHandleList", no_init)
    .def(vector_indexing_suite<ChainHandleList>())
    .def(geom::VectorAdditions<ChainHandleList>())    
  ;

  {
    enum_<ChainType>("ChainType")
      .value("CHAINTYPE_POLY",                  CHAINTYPE_POLY)
      .value("CHAINTYPE_NON_POLY",              CHAINTYPE_NON_POLY)
      .value("CHAINTYPE_WATER",                 CHAINTYPE_WATER)
      .value("CHAINTYPE_POLY_PEPTIDE_D",        CHAINTYPE_POLY_PEPTIDE_D)
      .value("CHAINTYPE_POLY_PEPTIDE_L",        CHAINTYPE_POLY_PEPTIDE_L)
      .value("CHAINTYPE_POLY_DN",               CHAINTYPE_POLY_DN)
      .value("CHAINTYPE_POLY_RN",               CHAINTYPE_POLY_RN)
      .value("CHAINTYPE_POLY_SAC_D",            CHAINTYPE_POLY_SAC_D)
      .value("CHAINTYPE_POLY_SAC_L",            CHAINTYPE_POLY_SAC_L)
      .value("CHAINTYPE_POLY_DN_RN",            CHAINTYPE_POLY_DN_RN)
      .value("CHAINTYPE_UNKNOWN",               CHAINTYPE_UNKNOWN)
      .value("CHAINTYPE_MACROLIDE",             CHAINTYPE_MACROLIDE)
      .value("CHAINTYPE_CYCLIC_PSEUDO_PEPTIDE", CHAINTYPE_CYCLIC_PSEUDO_PEPTIDE)
      .value("CHAINTYPE_POLY_PEPTIDE_DN_RN",    CHAINTYPE_POLY_PEPTIDE_DN_RN)
      .value("CHAINTYPE_BRANCHED",              CHAINTYPE_BRANCHED)
      .value("CHAINTYPE_OLIGOSACCHARIDE",       CHAINTYPE_OLIGOSACCHARIDE)
      .value("CHAINTYPE_N_CHAINTYPES",          CHAINTYPE_N_CHAINTYPES)
      .export_values()
    ;
  }

  def("ChainTypeFromString", ChainTypeFromStringPtr);
  def("StringFromChainType", &StringFromChainType);
}
