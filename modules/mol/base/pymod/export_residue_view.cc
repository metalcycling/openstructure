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
#include <ost/geom/export_helper/vector.hh>
using namespace ost;
using namespace ost::mol;
#include "bounds.hh"

namespace {
typedef AtomView (ResidueView::*StringMethod)(const String&) const;
typedef AtomView (ResidueView::*HandleMethod)(const AtomHandle&) const;
typedef AtomView (ResidueView::*ViewMethod)(const AtomView&, ViewAddFlags);
typedef AtomView (ResidueView::*HandleMethodNonConst)(const AtomHandle&, ViewAddFlags);

StringMethod string_find_atom=&ResidueView::FindAtom;
HandleMethod handle_find_atom=&ResidueView::FindAtom;
HandleMethod view_for_ah=&ResidueView::ViewForHandle;
HandleMethodNonConst add_atom_handle=&ResidueView::AddAtom;
ViewMethod add_atom_view=&ResidueView::AddAtom;
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(X_add_atom_overloads,
                                       ResidueView::AddAtom, 1, 2);
typedef EntityView (ResidueView::*QueryMethod2)(const Query&, uint) const;
typedef EntityView (ResidueView::*StringMethod2)(const String&, uint) const;
QueryMethod2 select_query=&ResidueView::Select;
StringMethod2 select_string=&ResidueView::Select;
}

void export_ResidueView()
{
  class_<ResidueViewList>("ResidueViewList", no_init)
    .def(vector_indexing_suite<ResidueViewList>())
    .def(geom::VectorAdditions<ResidueViewList>()) 
  ;

  class_<ResidueView, bases<ResidueBase> >("ResidueView", init<>())
    .def("GetChain",&ResidueView::GetChain)
    .add_property("chain", &ResidueView::GetChain)
    .def("GetEntity", &ResidueView::GetEntity)
    .add_property("entity", &ResidueView::GetEntity)
    .def("GetAtomList", &ResidueView::GetAtomList,
         return_value_policy<copy_const_reference>())
    .def("GetAtomCount", &ResidueView::GetAtomCount)
    .add_property("atoms",
                  make_function(&ResidueView::GetAtomList,
                                return_value_policy<copy_const_reference>()))
    .add_property("atom_count", &ResidueView::GetAtomCount)
    .def("FindAtom", string_find_atom, args("atom_name"))
    .def("AddAtom", add_atom_handle, X_add_atom_overloads(args("atom_handle", "flags")))
    .def("AddAtom", add_atom_view, X_add_atom_overloads(args("atom_view", "flags")))
    .def("FindAtom", handle_find_atom, args("atom_handle"))
    .def("ViewForHandle", view_for_ah, arg("atom_handle"))
    .def("IsAtomIncluded", &ResidueView::IsAtomIncluded, args("atom_handle"))
    .def("GetIndex", &ResidueView::GetIndex)  
    .add_property("index", &ResidueView::GetIndex)
    .def("RemoveAtoms", &ResidueView::RemoveAtom)
    .def("IsValid", &ResidueView::IsValid)
    .add_property("handle", &ResidueView::GetHandle)
    .def("GetHandle", &ResidueView::GetHandle)
    .def("Select", select_string, arg("flags")=0)
    .def("Select", select_query, arg("flags")=0)
    .def("RemoveAtom", &ResidueView::RemoveAtom)
    .def(self==self)
    .def(self!=self)    
    .def("GetMass", &ResidueView::GetMass)
    .def("GetCenterOfMass", &ResidueView::GetCenterOfMass)
    .def("GetCenterOfAtoms", &ResidueView::GetCenterOfAtoms)
    .def("GetGeometricCenter", geom_center<ResidueView>)
    .add_property("mass", &ResidueView::GetMass)
    .add_property("center_of_mass", &ResidueView::GetCenterOfMass)
    .add_property("center_of_atoms", &ResidueView::GetCenterOfAtoms)    
    .add_property("geometric_center", geom_center<ResidueView>)
    .add_property("valid", &ResidueView::IsValid)    
    .def("GetGeometricStart", geom_start<ResidueView>)
    .def("GetGeometricEnd", geom_end<ResidueView>) 
    .def("GetBounds", &ResidueView::GetBounds)
    .add_property("bounds", &ResidueView::GetBounds)
    .def("__hash__", &ResidueView::GetHashCode)
    .def("GetHashCode", &ResidueView::GetHashCode)
    .add_property("hash_code", &ResidueView::GetHashCode)
  ;


}
