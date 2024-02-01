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
#include <boost/shared_ptr.hpp>
#include <ost/mol/mm/forcefield.hh>

using namespace boost::python;

namespace {

ost::mol::mm::InteractionPtr GetLJOneType(ost::mol::mm::ForcefieldPtr p,
                                          String type) {
  return p->GetLJ(type);
}

ost::mol::mm::InteractionPtr GetLJTwoTypes(ost::mol::mm::ForcefieldPtr p,
                                           String type1, String type2,
                                           bool pair) {
  return p->GetLJ(type1, type2, pair);
}

boost::python::list WrapGetAtomRenamingRules(ost::mol::mm::ForcefieldPtr p,
                                             const String& res_name) {
  boost::python::list result;
  typedef ost::mol::mm::Forcefield::AtomRenamingType AtomRenamingType;
  const AtomRenamingType& rules = p->GetAtomRenamingRules(res_name);
  for (uint i = 0; i < rules.size(); ++i) {
    result.append(boost::python::make_tuple(rules[i].first, rules[i].second));
  }
  return result;
}

boost::python::list WrapGetBuildingBlockNames(ost::mol::mm::ForcefieldPtr p) {
  boost::python::list result;
  std::vector<String> v_result = p->GetBuildingBlockNames();
  for(std::vector<String>::iterator it = v_result.begin();
      it != v_result.end(); ++it) {
    result.append(*it);
  }
  return result;
}

}


void export_Forcefield()
{

  class_<ost::mol::mm::Forcefield>("Forcefield", init<>())
    .def("Load",&ost::mol::mm::Forcefield::Load,(arg("filename"))).staticmethod("Load")
    .def("Save",&ost::mol::mm::Forcefield::Save,(arg("filename")))
    .def("GetBuildingBlock",&ost::mol::mm::Forcefield::GetBuildingBlock,(arg("name")))
    .def("GetBlockModifier",&ost::mol::mm::Forcefield::GetBlockModifier,(arg("name")))
    .def("GetBuildingBlockNames",&WrapGetBuildingBlockNames)
    .def("GetAtomType",&ost::mol::mm::Forcefield::GetAtomType,(arg("res_name"),arg("atom_name")))
    .def("GetHydrogenConstructor",&ost::mol::mm::Forcefield::GetHydrogenConstructor,(arg("name")))
    .def("GetNTerModifier",&ost::mol::mm::Forcefield::GetNTerModifier,(arg("residue_name"),arg("ter_name")=""))
    .def("GetCTerModifier",&ost::mol::mm::Forcefield::GetCTerModifier,(arg("residue_name"),arg("ter_name")=""))
    .def("GetBond",&ost::mol::mm::Forcefield::GetBond,(arg("atom_type_1"),arg("atom_type_2")))
    .def("GetAngle",&ost::mol::mm::Forcefield::GetAngle,(arg("atom_type_1"),arg("atom_type_2"),arg("atom_type_3")))
    .def("GetDihedrals",&ost::mol::mm::Forcefield::GetDihedrals,(arg("atom_type_1"),arg("atom_type_2"),arg("atom_type_3"),arg("atom_type_4")))
    .def("GetImpropers",&ost::mol::mm::Forcefield::GetImpropers,(arg("atom_type_1"),arg("atom_type_2"),arg("atom_type_3"),arg("atom_type_4")))
    .def("GetCMap",&ost::mol::mm::Forcefield::GetCMap,(arg("atom_type_1"),arg("atom_type_2"),arg("atom_type_3"),arg("atom_type_4"),arg("atom_type_5")))
    .def("GetImplicitGenborn",&ost::mol::mm::Forcefield::GetImplicitGenborn,(arg("atom_type")))
    .def("GetLJ",&GetLJOneType,(arg("atom_type")))
    .def("GetLJ",&GetLJTwoTypes,(arg("atom_type_1"),arg("atom_type_2")))
    .def("GetConstraint",&ost::mol::mm::Forcefield::GetConstraint,(arg("atom_type_1"),arg("atom_type_2")))
    .def("GetMass",&ost::mol::mm::Forcefield::GetMass,(arg("atom_type")))
    .def("GetFudgeLJ",&ost::mol::mm::Forcefield::GetFudgeLJ)
    .def("GetFudgeQQ",&ost::mol::mm::Forcefield::GetFudgeQQ)
    .def("AddBuildingBlock",&ost::mol::mm::Forcefield::AddBuildingBlock,(arg("name"),arg("block")))
    .def("AddBond",&ost::mol::mm::Forcefield::AddBond,(arg("bond"))) 
    .def("AddAngle",&ost::mol::mm::Forcefield::AddAngle,(arg("angle"))) 
    .def("AddDihedral",&ost::mol::mm::Forcefield::AddDihedral,(arg("dihedral"))) 
    .def("AddImproper",&ost::mol::mm::Forcefield::AddImproper,(arg("improper"))) 
    .def("AddCMap",&ost::mol::mm::Forcefield::AddCMap,(arg("cmap"))) 
    .def("AddImplicitGenborn",&ost::mol::mm::Forcefield::AddImplicitGenborn,(arg("implicit_genborn"))) 
    .def("AddLJ",&ost::mol::mm::Forcefield::AddLJ,(arg("lj"))) 
    .def("AddLJPair",&ost::mol::mm::Forcefield::AddLJPair,(arg("lj_pair"))) 
    .def("AddMass",&ost::mol::mm::Forcefield::AddMass,(arg("mass"))) 
    .def("AddResidueRenamingRule",&ost::mol::mm::Forcefield::AddResidueRenamingRule,(arg("name"),arg("main_name"),arg("n_ter_name"),arg("c_ter_name"),arg("two_ter_name"))) 
    .def("AddAtomRenamingRule",&ost::mol::mm::Forcefield::AddAtomRenamingRule,(arg("res_name"),arg("old_atom_name"),arg("new_atom_name"))) 
    .def("AddHydrogenConstructor",&ost::mol::mm::Forcefield::AddHydrogenConstructor,(arg("name"),arg("hydrogen_constructor"))) 
    .def("AddBlockModifier",&ost::mol::mm::Forcefield::AddBlockModifier,(arg("name"),arg("modifier"))) 
    .def("SetStandardCTer",&ost::mol::mm::Forcefield::SetStandardCTer,(arg("res_name"),arg("ter_name")))
    .def("SetStandardNTer",&ost::mol::mm::Forcefield::SetStandardNTer,(arg("res_name"),arg("ter_name")))
    .def("SetFudgeLJ",&ost::mol::mm::Forcefield::SetFudgeLJ,(arg("fudge_factor"))) 
    .def("SetFudgeQQ",&ost::mol::mm::Forcefield::SetFudgeQQ,(arg("fudge_factor"))) 
    .def("SetGenPairs",&ost::mol::mm::Forcefield::SetGenPairs,(arg("set_gen_pairs_flag"))) 
    .def("GetResidueRenamingMain",&ost::mol::mm::Forcefield::GetResidueRenamingMain,(arg("res_name")))
    .def("GetResidueRenamingNTer",&ost::mol::mm::Forcefield::GetResidueRenamingNTer,(arg("res_name")))
    .def("GetResidueRenamingCTer",&ost::mol::mm::Forcefield::GetResidueRenamingCTer,(arg("res_name")))
    .def("GetResidueRenamingTwoTer",&ost::mol::mm::Forcefield::GetResidueRenamingTwoTer,(arg("res_name")))
    .def("GetAtomRenaming",&ost::mol::mm::Forcefield::GetAtomRenaming,(arg("res_name"),arg("atom_name")))
    .def("AssignFFSpecificNames",&ost::mol::mm::Forcefield::AssignFFSpecificNames,(arg("ent"),arg("reverse")=false))
    .def("HasAtomRenamingRules", &ost::mol::mm::Forcefield::HasAtomRenamingRules,
         (arg("res_name")))
    .def("GetAtomRenamingRules", &WrapGetAtomRenamingRules, (arg("res_name")))
  ;

  boost::python::register_ptr_to_python<ost::mol::mm::ForcefieldPtr>();
  
}
