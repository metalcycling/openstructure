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
#include <ost/mol/mm/block_modifiers.hh>
#include <ost/mol/mm/gromacs_block_modifiers.hh>
#include <ost/mol/mm/heuristic_block_modifiers.hh>

using namespace boost::python;

namespace{

template <typename Iterator>
void AppendToList(Iterator first, Iterator last, boost::python::list& l) {
  for (Iterator i = first; i != last; ++i) {
    l.append(*i);
  }
}

template <typename Vec>
void AppendVectorToList(const Vec& v, boost::python::list& l) {
  AppendToList(v.begin(), v.end(), l);
}

template<typename T>
std::vector<T> ListToVec(const boost::python::list& l){
  std::vector<T> vec;
  for (int i = 0; i < boost::python::len(l); ++i){
    vec.push_back(boost::python::extract<T>(l[i]));
  }
  return vec;
}

void WrapAddAddRule(ost::mol::mm::GromacsBlockModifierPtr p,
                    int number, int method, 
                    const boost::python::list& atom_names,
                    const boost::python::list& anchors, 
                    const String& type, Real new_charge){
  std::vector<String> v_atom_names = ListToVec<String>(atom_names);
  std::vector<String> v_anchors = ListToVec<String>(anchors);

  p->AddAddRule(number,method,v_atom_names,v_anchors,type,new_charge);

}

void WrapAddHydrogenRule(ost::mol::mm::GromacsHydrogenConstructorPtr p,
                         int number, int method, 
                         const boost::python::list& hydrogen_names, 
                         const boost::python::list& anchors){
  std::vector<String> v_hydrogen_names = ListToVec<String>(hydrogen_names);
  std::vector<String> v_anchors = ListToVec<String>(anchors);
  p->AddHydrogenRule(number,method,v_hydrogen_names,v_anchors);
}

boost::python::list
WrapGetHydrogenRules(ost::mol::mm::GromacsHydrogenConstructorPtr p) {
  // for data extraction: get all rules as list of tuples
  // -> data: (number, method, names, anchors)
  boost::python::list result;
  for (uint i = 0; i < p->GetNumHydrogenRules(); ++i) {
    int number = p->GetHydrogenRuleNumber(i);
    int method = p->GetHydrogenRuleMethod(i);
    boost::python::list names;
    AppendVectorToList(p->GetHydrogenRuleNames(i), names);
    boost::python::list anchors;
    AppendVectorToList(p->GetHydrogenRuleAnchors(i), anchors);
    result.append(boost::python::make_tuple(number, method, names, anchors));
  }
  return result;
}

boost::python::list WrapGetBonds(ost::mol::mm::GromacsBlockModifierPtr p) {
  boost::python::list result;
  AppendVectorToList(p->GetBonds(), result);
  return result;
}

boost::python::list WrapGetAngles(ost::mol::mm::GromacsBlockModifierPtr p) {
  boost::python::list result;
  AppendVectorToList(p->GetAngles(), result);
  return result;
}

boost::python::list WrapGetDihedrals(ost::mol::mm::GromacsBlockModifierPtr p) {
  boost::python::list result;
  AppendVectorToList(p->GetDihedrals(), result);
  return result;
}

boost::python::list WrapGetImpropers(ost::mol::mm::GromacsBlockModifierPtr p) {
  boost::python::list result;
  AppendVectorToList(p->GetImpropers(), result);
  return result;
}

boost::python::list WrapGetCmaps(ost::mol::mm::GromacsBlockModifierPtr p) {
  boost::python::list result;
  AppendVectorToList(p->GetCmaps(), result);
  return result;
}

boost::python::list WrapGetDeleteAtoms(ost::mol::mm::GromacsBlockModifierPtr p) {
  boost::python::list result;
  AppendVectorToList(p->GetDeleteAtoms(), result);
  return result;
}

boost::python::list
WrapGetReplaceRules(ost::mol::mm::GromacsBlockModifierPtr p) {
  // for data extraction: get all rules as list of tuples
  // -> data: (name, new_name, new_type, new_charge)
  boost::python::list result;
  for (uint i = 0; i < p->GetNumReplaceRules(); ++i) {
    result.append(boost::python::make_tuple(p->GetReplaceRuleName(i),
                                            p->GetReplaceRuleNewName(i),
                                            p->GetReplaceRuleNewType(i),
                                            p->GetReplaceRuleNewCharge(i)));
  }
  return result;
}

boost::python::list
WrapGetAddRules(ost::mol::mm::GromacsBlockModifierPtr p) {
  // for data extraction: get all rules as list of tuples
  // -> data: (number, method, names, anchors, type, charge)
  boost::python::list result;
  for (uint i = 0; i < p->GetNumAddRules(); ++i) {
    boost::python::list names;
    AppendVectorToList(p->GetAddRuleNames(i), names);
    boost::python::list anchors;
    AppendVectorToList(p->GetAddRuleAnchors(i), anchors);
    result.append(boost::python::make_tuple(p->GetAddRuleNumber(i),
                                            p->GetAddRuleMethod(i),
                                            names, anchors,
                                            p->GetAddRuleType(i),
                                            p->GetAddRuleCharge(i)));
  }
  return result;
}

}


void export_BlockModifiers()
{

  class_<ost::mol::mm::HydrogenConstructor, boost::noncopyable>("HydrogenConstructor",no_init);
  boost::python::register_ptr_to_python<ost::mol::mm::HydrogenConstructorPtr>();

  class_<ost::mol::mm::BlockModifier, boost::noncopyable>("BlockModifier",no_init);
  boost::python::register_ptr_to_python<ost::mol::mm::BlockModifierPtr>();

  class_<ost::mol::mm::GromacsHydrogenConstructor, bases<ost::mol::mm::HydrogenConstructor> >("GromacsHydrogenConstructor", init<>())
    .def("ApplyOnBuildingBlock",&ost::mol::mm::GromacsHydrogenConstructor::ApplyOnBuildingBlock,(arg("block")))
    .def("ApplyOnResidue",&ost::mol::mm::GromacsHydrogenConstructor::ApplyOnResidue,(arg("residue"),arg("editor")))
    .def("AddHydrogenRule",&WrapAddHydrogenRule,(arg("number"),arg("method"),arg("hydrogen_names"),arg("anchors")))
    .def("GetHydrogenRules",&WrapGetHydrogenRules)
  ;

  class_<ost::mol::mm::GromacsBlockModifier, bases<ost::mol::mm::BlockModifier> >("GromacsBlockModifier", init<>())
    .def("ApplyOnBuildingBlock",&ost::mol::mm::GromacsBlockModifier::ApplyOnBuildingBlock,(arg("block")))
    .def("ApplyOnResidue",&ost::mol::mm::GromacsBlockModifier::ApplyOnResidue,(arg("residue"),arg("editor")))
    .def("AddReplaceRule",&ost::mol::mm::GromacsBlockModifier::AddReplaceRule,(arg("name"),arg("new_name"),arg("new_type"),arg("new_charge")))
    .def("AddAddRule",&WrapAddAddRule,(arg("number"),arg("method"),arg("atom_names"),arg("anchors"),arg("type"),arg("charge")))
    .def("AddBond",&ost::mol::mm::GromacsBlockModifier::AddBond,(arg("bond")))
    .def("AddAngle",&ost::mol::mm::GromacsBlockModifier::AddAngle,(arg("angle")))
    .def("AddDihedral",&ost::mol::mm::GromacsBlockModifier::AddDihedral,(arg("dihedral")))
    .def("AddImproper",&ost::mol::mm::GromacsBlockModifier::AddImproper,(arg("improper")))
    .def("AddCMap",&ost::mol::mm::GromacsBlockModifier::AddCMap,(arg("cmap")))
    .def("AddDeleteAtom",&ost::mol::mm::GromacsBlockModifier::AddDeleteAtom,(arg("name")))
    .def("GetBonds", &WrapGetBonds)
    .def("GetAngles", &WrapGetAngles)
    .def("GetDihedrals", &WrapGetDihedrals)
    .def("GetImpropers", &WrapGetImpropers)
    .def("GetCmaps", &WrapGetCmaps)
    .def("GetDeleteAtoms", &WrapGetDeleteAtoms)
    .def("GetReplaceRules", &WrapGetReplaceRules)
    .def("GetAddRules", &WrapGetAddRules)
  ;
 
  class_<ost::mol::mm::HeuristicHydrogenConstructor, bases<ost::mol::mm::HydrogenConstructor> >("HeuristicHydrogenConstructor", init<ost::mol::mm::BuildingBlockPtr>())
  .def("ApplyOnBuildingBlock",&ost::mol::mm::HeuristicHydrogenConstructor::ApplyOnBuildingBlock,(arg("block")))
  .def("ApplyOnResidue",&ost::mol::mm::HeuristicHydrogenConstructor::ApplyOnResidue,(arg("residue"),arg("editor")))
  ;

  boost::python::register_ptr_to_python<ost::mol::mm::GromacsHydrogenConstructorPtr>();
  boost::python::register_ptr_to_python<ost::mol::mm::GromacsBlockModifierPtr>();
  boost::python::register_ptr_to_python<ost::mol::mm::HeuristicHydrogenConstructorPtr>();
  
}

