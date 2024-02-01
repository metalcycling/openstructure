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

#include <ost/mol/mm/interaction.hh>

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

template<typename T>
std::vector<T> ListToVec(boost::python::list& l){
  std::vector<T> vec;
  for (int i = 0; i < boost::python::len(l); ++i){
    vec.push_back(boost::python::extract<T>(l[i]));
  }
  return vec;
}

boost::python::list WrapGetNames(ost::mol::mm::InteractionPtr p){
  std::vector<String> names = p->GetNames();
  return VecToList<String>(names);
}

boost::python::list WrapGetTypes(ost::mol::mm::InteractionPtr p){
    std::vector<String> types = p->GetTypes();
  return VecToList<String>(types);
}

boost::python::list WrapGetParam(ost::mol::mm::InteractionPtr p){
  std::vector<Real> param = p->GetParam();
  return VecToList<Real>(param);
}

void WrapSetNames(ost::mol::mm::InteractionPtr p, boost::python::list l){
  std::vector<String> names = ListToVec<String>(l);
  p->SetNames(names); 
}

void WrapSetTypes(ost::mol::mm::InteractionPtr p, boost::python::list l){
  std::vector<String> types = ListToVec<String>(l);
  p->SetTypes(types); 
}

void WrapSetParam(ost::mol::mm::InteractionPtr p, boost::python::list l){
  std::vector<Real> param = ListToVec<Real>(l);
  p->SetParam(param); 
}

}

void export_Interaction()
{

  enum_<ost::mol::mm::FuncType>("FuncType")
    .value("Invalid", ost::mol::mm::Invalid)
    .value("HarmonicBond", ost::mol::mm::HarmonicBond)
    .value("HarmonicAngle", ost::mol::mm::HarmonicAngle)
    .value("UreyBradleyAngle", ost::mol::mm::UreyBradleyAngle)
    .value("PeriodicDihedral", ost::mol::mm::PeriodicDihedral)
    .value("PeriodicImproper", ost::mol::mm::PeriodicImproper)
    .value("HarmonicImproper", ost::mol::mm::HarmonicImproper)
    .value("CMap", ost::mol::mm::CMap)
    .value("LJ", ost::mol::mm::LJ)
    .value("LJPair", ost::mol::mm::LJPair)
    .value("GBSA", ost::mol::mm::GBSA)
    .value("DistanceConstraint", ost::mol::mm::DistanceConstraint)
    .value("Exclusion", ost::mol::mm::Exclusion)
    .value("HarmonicPositionRestraint", ost::mol::mm::HarmonicPositionRestraint)
    .value("HarmonicDistanceRestraint", ost::mol::mm::HarmonicDistanceRestraint)
  ;

  class_<ost::mol::mm::Interaction>("Interaction",init<ost::mol::mm::FuncType>())
    .def("SetTypes",&WrapSetTypes,(arg("types")))
    .def("SetNames",&WrapSetNames,(arg("names")))
    .def("SetParam",&WrapSetParam,(arg("parameters")))
    .def("GetTypes",&WrapGetTypes)
    .def("GetNames",&WrapGetNames)
    .def("GetParam",&WrapGetParam)
    .def("GetAtoms",&ost::mol::mm::Interaction::GetAtoms,(arg("residue")))
    .def("GetFuncType",&ost::mol::mm::Interaction::GetFuncType)
    .def("ReplaceAtom",&ost::mol::mm::Interaction::ReplaceAtom,(arg("name"),arg("new_name"),arg("new_type")))
    .def("MatchTypes",&ost::mol::mm::Interaction::MatchTypes,(arg("types")))
    .def("MatchNames",&ost::mol::mm::Interaction::MatchNames,(arg("names")))
    .def("HasName",&ost::mol::mm::Interaction::HasName,(arg("name")))
    .def("HasType",&ost::mol::mm::Interaction::HasType,(arg("type")))
    .def("IsParametrized",&ost::mol::mm::Interaction::IsParametrized)
    .def("HasTypeWildcard",&ost::mol::mm::Interaction::HasTypeWildcard)
    .def("HasNameWildcard",&ost::mol::mm::Interaction::HasNameWildcard)    
  ;

  boost::python::register_ptr_to_python<ost::mol::mm::InteractionPtr>();

  class_<std::vector<ost::mol::mm::InteractionPtr> >("InteractionList", init<>())
    .def(vector_indexing_suite<std::vector<ost::mol::mm::InteractionPtr>, true>())
  ;

}
