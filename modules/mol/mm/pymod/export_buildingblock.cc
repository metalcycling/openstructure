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

#include <limits>
#include <boost/python.hpp>
#include <ost/mol/mm/buildingblock.hh>
#include <ost/log.hh>
#include <vector>

using namespace boost::python;

namespace{

template<typename T>
boost::python::list VecToList(std::vector<T>& vec){
  boost::python::list l;
  for(typename std::vector<T>::iterator it=vec.begin();it!=vec.end();++it){
    l.append(*it);
  }
  return l;
}

bool WrapMatch(ost::mol::mm::BuildingBlockPtr p, const ost::mol::ResidueHandle& res, bool match_connectivity){
  String match_fail_info;
  bool match = p->Match(res,match_connectivity,match_fail_info);
  if(!match){
    std::stringstream ss;
    ss << "Residue "<< res.GetQualifiedName() << " does not match the building";
    ss << "block. "<<match_fail_info;
    LOG_INFO(ss.str());
  }
  return match;
}

boost::python::list WrapGetAtoms(ost::mol::mm::BuildingBlockPtr p){
  std::vector<String> atoms = p->GetAtoms();
  return VecToList<String>(atoms);
}

boost::python::list WrapGetTypes(ost::mol::mm::BuildingBlockPtr p){
    std::vector<String> types = p->GetTypes();
  return VecToList<String>(types);
}

boost::python::list WrapGetCharges(ost::mol::mm::BuildingBlockPtr p){
  std::vector<Real> charges = p->GetCharges();
  return VecToList<Real>(charges);
}

boost::python::list WrapGetMasses(ost::mol::mm::BuildingBlockPtr p){
  std::vector<Real> masses = p->GetMasses();
  return VecToList<Real>(masses);
}


}

void export_Buildingblock()
{
  class_<ost::mol::mm::BuildingBlock>("BuildingBlock",init<>())
    .def("Match",&WrapMatch,(arg("res"),arg("match_connectivity")=true))
    .def("Connect",&ost::mol::mm::BuildingBlock::Connect,(arg("residue"),arg("xcs_editor")))
    .def("GetAtoms",&WrapGetAtoms)
    .def("GetTypes",&WrapGetTypes)
    .def("GetCharges",&WrapGetCharges)
    .def("GetMasses",&WrapGetMasses)
    .def("GetType",&ost::mol::mm::BuildingBlock::GetType,(arg("name")))
    .def("GetCharge",&ost::mol::mm::BuildingBlock::GetCharge,(arg("name")))
    .def("GetMass",&ost::mol::mm::BuildingBlock::GetMass,(arg("name")))
    .def("GetBonds",&ost::mol::mm::BuildingBlock::GetBonds)
    .def("GetAngles",&ost::mol::mm::BuildingBlock::GetAngles)
    .def("GetDihedrals",&ost::mol::mm::BuildingBlock::GetDihedrals)
    .def("GetImpropers",&ost::mol::mm::BuildingBlock::GetImpropers)
    .def("GetCMaps",&ost::mol::mm::BuildingBlock::GetCMaps)
    .def("GetExclusions",&ost::mol::mm::BuildingBlock::GetExclusions)
    .def("GetConstraints",&ost::mol::mm::BuildingBlock::GetConstraints)
    .def("AddAtom",&ost::mol::mm::BuildingBlock::AddAtom,(arg("name"),arg("type"),arg("charge"),arg("mass")=std::numeric_limits<Real>::quiet_NaN()))
    .def("AddBond",&ost::mol::mm::BuildingBlock::AddBond,(arg("bond"),arg("replace_existing")=false))
    .def("AddAngle",&ost::mol::mm::BuildingBlock::AddAngle,(arg("angle"),arg("replace_existing")=false))
    .def("AddDihedral",&ost::mol::mm::BuildingBlock::AddDihedral,(arg("dihedral"),arg("replace_existing")=false))
    .def("AddImproper",&ost::mol::mm::BuildingBlock::AddImproper,(arg("improper"),arg("replace_existing")=false))
    .def("AddExclusion",&ost::mol::mm::BuildingBlock::AddExclusion,(arg("exclusion"),arg("replace_existing")=false))
    .def("AddCMap",&ost::mol::mm::BuildingBlock::AddCMap,(arg("cmap"),arg("replace_existing")=false))
    .def("AddConstraint",&ost::mol::mm::BuildingBlock::AddConstraint,(arg("constraint"),arg("replace_existing")=false))
    .def("RemoveAtom",&ost::mol::mm::BuildingBlock::RemoveAtom,(arg("name")))
    .def("ReplaceAtom",&ost::mol::mm::BuildingBlock::ReplaceAtom,(arg("name"),arg("new_name"),arg("new_type"),arg("new_charge"),arg("new_mass")=std::numeric_limits<Real>::quiet_NaN()))
    .def("RemoveInteractionsToNext",&ost::mol::mm::BuildingBlock::RemoveInteractionsToNext)
    .def("RemoveInteractionsToPrev",&ost::mol::mm::BuildingBlock::RemoveInteractionsToPrev)
  ;

  boost::python::register_ptr_to_python<ost::mol::mm::BuildingBlockPtr>();

}
