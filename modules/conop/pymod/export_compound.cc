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
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;

#include <ost/message.hh>
#include <ost/conop/compound.hh>
#include <ost/conop/compound_lib.hh>
using namespace ost::mol;
using namespace ost::conop;

namespace {

Compound::Dialect tr_dialect(const String& dialect)
{
  if (dialect=="PDB") {
    return Compound::PDB;
  }
  if (dialect=="CHARMM") {
    return Compound::CHARMM;
  }
  if (dialect=="OPLS") {
    return Compound::OPLS;
  }
  if (dialect=="AMBER") {
    return Compound::AMBER;
  }
  std::stringstream ss;
  ss << "unknown compound dialect '" << dialect << "'";
  throw ost::Error(ss.str());
}

void set_dialect(CompoundPtr compound, const String& dialect)
{
  compound->SetDialect(tr_dialect(dialect));
}

char get_chemclass(CompoundPtr compound)
{
  return char(compound->GetChemClass());
}


void set_chemclass(CompoundPtr compound, char cc)
{
  compound->SetChemClass(ChemClass(cc));
}

char get_chemtype(CompoundPtr compound)
{
  return char(compound->GetChemType());
}

CompoundPtr find_compound(CompoundLibPtr comp_lib, 
                          const String& tlc, const String& dialect)
{
  return comp_lib->FindCompound(tlc, tr_dialect(dialect));
}

bool is_residue_complete(CompoundLibPtr comp_lib,
                         const ost::mol::ResidueHandle& res,
                         bool check_hydrogens, const String& dialect)
{
  return comp_lib->IsResidueComplete(res, check_hydrogens, tr_dialect(dialect));
}

String get_creation_date(CompoundLibPtr comp_lib)
{
  return comp_lib->GetCreationDate().ToString();
}

}
void export_Compound() {

  class_<Compound, CompoundPtr>("Compound", no_init)
    .def("GetID", &Compound::GetID, 
         return_value_policy<copy_const_reference>())
    .def("SetOneLetterCode", &Compound::SetOneLetterCode)
    .def("GetOneLetterCode", &Compound::GetOneLetterCode)

    .add_property("three_letter_code", make_function(&Compound::GetID, return_value_policy<copy_const_reference>()))
    .add_property("name", 
                  make_function(&Compound::GetName, 
                 return_value_policy<copy_const_reference>()))
    .add_property("id", make_function(&Compound::GetID, return_value_policy<copy_const_reference>()))    
    .add_property("one_letter_code", &Compound::GetOneLetterCode, 
                  &Compound::SetOneLetterCode)                  
    .def("GetAtomSpecs", &Compound::GetAtomSpecs,
         return_value_policy<copy_const_reference>())
    .add_property("bond_specs", make_function(&Compound::GetBondSpecs,
         return_value_policy<copy_const_reference>()))
    .add_property("atom_specs", make_function(&Compound::GetAtomSpecs,
         return_value_policy<copy_const_reference>()))
    .def("AddAtom", &Compound::AddAtom)
    .def("AddBond", &Compound::AddBond)
    .def("IsPeptideLinking", &Compound::IsPeptideLinking)
    .add_property("chem_class", &get_chemclass,
                  &set_chemclass)
    .add_property("chem_type", &get_chemtype)
    .add_property("formula",make_function(&Compound::GetFormula, 
                  return_value_policy<copy_const_reference>()),
                  &Compound::SetFormula)
    .add_property("dialect", &Compound::GetDialectAsString, 
                  &set_dialect)
    .add_property("inchi",
                  make_function(&Compound::GetInchi, 
                                return_value_policy<copy_const_reference>()))
    .add_property("inchi_key",
                  make_function(&Compound::GetInchiKey, 
                                return_value_policy<copy_const_reference>()))
    .add_property("smiles",
                  make_function(&Compound::GetSMILES,
                                return_value_policy<copy_const_reference>()))
  ;
  
  class_<AtomSpec>("AtomSpec", no_init)
    .def_readonly("element", &AtomSpec::element)
    .def_readonly("name", &AtomSpec::name)
    .def_readonly("alt_name", &AtomSpec::alt_name)
    .def_readonly("is_leaving", &AtomSpec::is_leaving)
    .def_readonly("is_aromatic", &AtomSpec::is_aromatic)
    .def_readonly("ordinal", &AtomSpec::ordinal)
    .def_readonly("charge", &AtomSpec::charge)
  ;
  
  class_<BondSpec>("BondSpec", no_init)
    .def_readonly("atom_one", &BondSpec::atom_one)
    .def_readonly("atom_two", &BondSpec::atom_two)    
    .def_readonly("order", &BondSpec::order)

  ;  

  class_<CompoundLib>("CompoundLib", no_init)
    .def("Load", &CompoundLib::Load, arg("readonly")=true).staticmethod("Load")
    .def("FindCompound", &find_compound, 
         (arg("tlc"), arg("dialect")="PDB"))
    .def("IsResidueComplete", &is_residue_complete, (arg("residue"), 
                                                     arg("check_hydrogens")=false,
                                                     arg("dialect")="PDB"))
    .def("ClearCache", &CompoundLib::ClearCache)
    .def("GetOSTVersionUsed", &CompoundLib::GetOSTVersionUsed)
    .def("GetCreationDate", &get_creation_date, (arg("comp_lib")))
  ;
  
  class_<AtomSpecList>("AtomSpecList", init<>())
    .def(vector_indexing_suite<AtomSpecList>())
  ;  
  
  class_<BondSpecList>("BondSpecList", init<>())
    .def(vector_indexing_suite<BondSpecList>())
  ;  
  register_ptr_to_python<CompoundLibPtr>();    
}
