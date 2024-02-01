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
#include <boost/version.hpp>
#if BOOST_VERSION<103400
#include <boost/python/detail/api_placeholder.hpp>
#endif
using namespace boost::python;

#include <ost/export_helper/pair_to_tuple_conv.hh>

#include <ost/mol/entity_view.hh>
#include <ost/mol/query.hh>
#include <ost/mol/mol.hh>
#include <ost/mol/view_op.hh>
using namespace ost;
using namespace ost::mol;
#include "bounds.hh"

namespace {

template<class T>
std::vector<T> from_list(const object& seq)
{
  std::vector<T> nrvo;
  for (int i = 0; i < len(seq); ++i) {
    nrvo.push_back(extract<T>(seq[i]));
  }
  return nrvo;
}


typedef ChainView (EntityView::*StringMethod)(const String&) const;
typedef ChainView (EntityView::*StringMethod)(const String&) const;
typedef EntityView (EntityView::*QueryMethod)(const Query&, uint) const;
typedef EntityView (EntityView::*QSMethod)(const String&, uint) const;
typedef ChainView (EntityView::*HandleMethod)(const ChainHandle&) const;
typedef ChainView (EntityView::*HM)(const ChainHandle&, ViewAddFlags);
typedef ChainView (EntityView::*VM)(const ChainView&, ViewAddFlags);
VM vm=&EntityView::AddChain;
HM hm=&EntityView::AddChain;

HandleMethod find_chain_hnd=&EntityView::FindChain;
HandleMethod view_for_ch=&EntityView::ViewForHandle;
StringMethod find_chain_str=&EntityView::FindChain;
QSMethod select_string=&EntityView::Select;
QueryMethod  select_query=&EntityView::Select;

EntityView create_view(const object& seq)
{
  
  if(len(seq)==0) {
    return EntityView();
  }
  extract<AtomHandle> get_handle(seq[0]);
  if(get_handle.check()) {
    return CreateViewFromAtomList(from_list<AtomHandle>(seq));
  }
  extract<AtomView> get_view(seq[0]);
  if(get_view.check()) {
    return CreateViewFromAtomList(from_list<AtomView>(seq));
  }
  throw Error("expected sequence of atom handles or atom views");
  return EntityView();
}

ResidueView (EntityView::*add_res_a)(const ResidueHandle&, 
                                     ViewAddFlags)=&EntityView::AddResidue;
ResidueView (EntityView::*add_res_b)(const ResidueView&, 
                                     ViewAddFlags)=&EntityView::AddResidue;
AtomView (EntityView::*add_atom_a)(const AtomHandle&, 
                                     ViewAddFlags)=&EntityView::AddAtom;
AtomView (EntityView::*add_atom_b)(const AtomView&, 
                                     ViewAddFlags)=&EntityView::AddAtom;                                     
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(X_add_residue_overloads,
                                      EntityView::AddResidue, 1, 2)


BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(X_add_chain_overloads,
                                      EntityView::AddChain, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(X_add_atom_overloads,
                                      EntityView::AddAtom, 1, 2)

AtomView (EntityView::*find_atom_a)(const String&, const ResNum&,
                                    const String&) const=&EntityView::FindAtom;
AtomView (EntityView::*find_atom_b)(const AtomHandle&) const=&EntityView::FindAtom;
AtomView (EntityView::*view_for_ah)(const AtomHandle&) const=&EntityView::ViewForHandle;

ResidueView (EntityView::*find_res_a)(const ResidueHandle&) const=&EntityView::FindResidue;
ResidueView (EntityView::*find_res_b)(const String&,
                                      const ResNum&) const=&EntityView::FindResidue;

ResidueView (EntityView::*view_for_rh)(const ResidueHandle&) const=&EntityView::ViewForHandle;
}


void export_EntityView()
{
  enum_<ViewAddFlag::Type>("ViewAddFlag")
    .value("INCLUDE_ATOMS",ViewAddFlag::INCLUDE_ATOMS)
    .value("INCLUDE_RESIDUES",ViewAddFlag::INCLUDE_RESIDUES)
    .value("INCLUDE_CHAINS",ViewAddFlag::INCLUDE_CHAINS)
    .value("INCLUDE_ALL",ViewAddFlag::INCLUDE_ALL)
    .value("CHECK_DUPLICATES",ViewAddFlag::CHECK_DUPLICATES)
    .export_values()
  ;

  Real (EntityView::*get_angle1)(const AtomHandle&, const AtomHandle&, const AtomHandle&) const = &EntityView::GetAngle;
  Real (EntityView::*get_angle2)(const AtomView&, const AtomView&, const AtomView&) const = &EntityView::GetAngle;

  class_<EntityView, bases<EntityBase> >("EntityView", init<>())
    .def("Copy", &EntityView::Copy)
    .def("ExtendViewToResidues", &EntityView::ExtendViewToResidues)
    .def("ExtendViewToSurrounding", &EntityView::ExtendViewToSurrounding)
    .def("FindChain", find_chain_str)
    .def("FindChain", find_chain_hnd)
    .def("ViewForHandle", view_for_rh)
    .def("ViewForHandle", view_for_ch)
    .def("FindAtom", find_atom_a)
    .def("FindAtom", find_atom_b)
    .def("FindResidue", find_res_a)
    .def("FindResidue", find_res_b)

    .def("ViewForHandle", view_for_ah)
    .def("GetAtomCount", &EntityView::GetAtomCount)
    .def("GetGeometricStart", geom_start<EntityView>)
    .def("GetGeometricEnd", geom_end<EntityView>)
    .def("GetGeometricCenter", geom_center<EntityView>)
    .add_property("geometric_center", geom_center<EntityView>)
    .def("GetCenterOfMass", &EntityView::GetCenterOfMass)
    .def("GetCenterOfAtoms", &EntityView::GetCenterOfAtoms)
    .def("GetMass", &EntityView::GetMass)
    .def("GetResidueCount", &EntityView::GetResidueCount)
    .def("GetChainCount", &EntityView::GetChainCount)
    .def("GetBondCount", &EntityView::GetBondCount)
    .def("GetAngle", get_angle1)
    .def("GetAngle", get_angle2)
    .def("FindChain", find_chain_hnd)
    .def("RemoveChain", &EntityView::RemoveChain)
    .def("FindWithin", &EntityView::FindWithin)
    .def("RemoveResidue", &EntityView::RemoveResidue)
    .def("RemoveAtom", &EntityView::RemoveAtom)
    .def("CreateEmptyView", &EntityView::CreateEmptyView)
    .def("IsValid", &EntityView::IsValid)    
    .def("CreateFullView", &EntityView::CreateFullView)
    .def("AddAllInclusiveBonds", &EntityView::AddAllInclusiveBonds)
    .add_property("chain_count", &EntityView::GetChainCount)
    .add_property("residue_count", &EntityView::GetResidueCount)
    .add_property("atom_count", &EntityView::GetAtomCount)
    .add_property("bond_count", &EntityView::GetBondCount)
    .add_property("mass", &EntityView::GetMass)
    .add_property("center_of_mass", &EntityView::GetCenterOfMass)
    .add_property("center_of_atoms", &EntityView::GetCenterOfAtoms)
    .add_property("valid", &EntityView::IsValid)    
    .def("GetResidueList", &EntityView::GetResidueList)
    .def("GetAtomList", &EntityView::GetAtomList)
    .add_property("atoms", &EntityView::GetAtomList)
    .add_property("residues", &EntityView::GetResidueList)
    .def("AddBond", &EntityView::AddBond)
    .def("AddAtom", add_atom_a,
         X_add_atom_overloads(args("atom_handle", "view_add_flags")))
    .def("AddAtom", add_atom_b,
         X_add_atom_overloads(args("atom_view", "view_add_flags")))         
    .def("AddResidue", add_res_a,
         X_add_residue_overloads(args("residue_handle", "view_add_flags")))
    .def("AddResidue", add_res_b,
         X_add_residue_overloads(args("residue_view", "view_add_flags")))         
    .def("AddChain", hm,
         X_add_chain_overloads(args("chain_handle", "view_add_flags")))
    .def("AddChain", vm,
         X_add_chain_overloads(args("chain_view", "view_add_flags")))
    .def("Select", select_string, arg("flags")=0)
    .def("Select", select_query, arg("flags")=0)
    .def("GetBondList", &EntityView::GetBondList,
         return_value_policy<reference_existing_object>())
    .def("AddAllInclusiveBonds",&EntityView::AddAllInclusiveBonds)
    .add_property("handle", &EntityView::GetHandle)
    .def("GetHandle", &EntityView::GetHandle)
    .add_property("bonds",
                  make_function(&EntityView::GetBondList,
                                return_value_policy<reference_existing_object>()))
    .def("GetChainList", &EntityView::GetChainList,
         return_value_policy<copy_const_reference>())
    .add_property("chains",
                  make_function(&EntityView::GetChainList,
                                return_value_policy<copy_const_reference>()))
    .def(self==self)
    .def(self!=self)
    .def("Dump", &EntityView::Dump)
    .def("GetBounds", &EntityView::GetBounds)
    .add_property("bounds", &EntityView::GetBounds)    
  ;
  def("Union", &Union);
  def("Difference", &Difference);
  def("Intersection", &Intersection);

  def("CreateViewFromAtoms", create_view);
  def("CreateViewFromAtomList", create_view);
  
  def("CreateEntityFromView", &CreateEntityFromView,
      (arg("view"), arg("include_exlusive_atoms"), arg("handle")=EntityHandle()));
  class_<EntityViewList>("EntityViewList", init<>())
    .def(vector_indexing_suite<EntityViewList>())
  ;

  to_python_converter<std::pair<EntityView, EntityView>,
                      PairToTupleConverter<EntityView, EntityView> >();
}
