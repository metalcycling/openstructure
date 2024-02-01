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

using namespace ost;
using namespace ost::mol;

#include <ost/export_helper/generic_property_def.hh>
#include <ost/geom/export_helper/vector.hh>

namespace {
  ChainHandle get_chain(AtomHandle& a)
  {
    return a.GetResidue().GetChain();
  }
}

void export_Atom()
{
  class_<AtomBase> atom_base("AtomBase", no_init);
  atom_base
    .def("GetPos", &AtomBase::GetPos,
         return_value_policy<copy_const_reference>())
    .def("GetOriginalPos", &AtomBase::GetOriginalPos,
         return_value_policy<copy_const_reference>())
    .def("GetName", &AtomBase::GetName,
         return_value_policy<copy_const_reference>())
    .def("GetQualifiedName", &AtomBase::GetQualifiedName)
    .add_property("qualified_name", &AtomBase::GetQualifiedName)
    .def(self_ns::str(self))
    .def("GetIndex", &AtomBase::GetIndex)
    .add_property("pos",
                  make_function(&AtomBase::GetPos,
                                return_value_policy<copy_const_reference>()))
    .add_property("original_pos",
                  make_function(&AtomBase::GetOriginalPos,
                                return_value_policy<copy_const_reference>()))
    .add_property("name",
                  make_function(&AtomBase::GetName,
                                return_value_policy<copy_const_reference>()),
                  &AtomBase::SetName)
    .add_property("index", &AtomBase::GetIndex)
    .def("GetRadius", &AtomBase::GetRadius)
    .def("SetRadius", &AtomBase::SetRadius)
    .def("GetElement", &AtomBase::GetElement, 
         return_value_policy<copy_const_reference>())
    .def("SetElement", &AtomBase::SetElement)
    .def("GetCharge", &AtomBase::GetCharge)
    .def("GetBFactor", &AtomBase::GetBFactor)
    .def("SetBFactor", &AtomBase::SetBFactor)
    .def("SetOccupancy", &AtomBase::SetOccupancy)
    .def("GetOccupancy", &AtomBase::GetOccupancy)
    .def("GetMass", &AtomBase::GetMass)
    .def("IsHetAtom", &AtomBase::IsHetAtom)
    .def("SetMass", &AtomBase::SetMass)
    .add_property("b_factor", &AtomBase::GetBFactor, &AtomBase::SetBFactor)
    .add_property("occupancy", &AtomBase::GetOccupancy, 
                  &AtomBase::SetOccupancy)
    .def("SetCharge", &AtomBase::SetCharge)
    .add_property("radius", &AtomBase::GetRadius, &AtomBase::SetRadius)
    .add_property("element", make_function(&AtomBase::GetElement, 
                  return_value_policy<copy_const_reference>()),
                  &AtomBase::SetElement)
    .add_property("is_hetatom", &AtomBase::IsHetAtom, &AtomBase::SetHetAtom)
    .add_property("charge", &AtomBase::GetCharge, &AtomBase::SetCharge)
    .add_property("mass", &AtomBase::GetMass, &AtomBase::SetMass)
    .add_property("valid", &AtomBase::IsValid)
    .def("IsValid", &AtomBase::IsValid)
  ;
  generic_prop_def<AtomBase>(atom_base);

  class_<AtomHandle, bases<AtomBase> >("AtomHandle", init<>())
    .def("GetResidue",&AtomHandle::GetResidue)
    .add_property("residue",&AtomHandle::GetResidue)
    .def("GetChain",get_chain)
    .add_property("chain",get_chain)
    .def("GetBondList", &AtomHandle::GetBondList)
    .def("GetBondCount", &AtomHandle::GetBondCount)
    .def("GetEntity", &AtomHandle::GetEntity)
    .add_property("bonds", &AtomHandle::GetBondList)
    .def("GetHandle", &AtomHandle::GetHandle)
    .add_property("handle", &AtomHandle::GetHandle)
    .add_property("entity", &AtomHandle::GetEntity)
    .def("GetBondPartners", &AtomHandle::GetBondPartners)
    .def("GetHashCode", &AtomHandle::GetHashCode)    
    .def("FindBondToAtom", &AtomHandle::FindBondToAtom, args("other_atom"))
    .def(self==self)
    .def(self!=self)
    .def("__hash__", &AtomHandle::GetHashCode)
    .add_property("hash_code", &AtomHandle::GetHashCode)
  ;

  class_<AtomHandleList>("AtomHandleList", no_init)
    .def(vector_indexing_suite<AtomHandleList>())
    .def(geom::VectorAdditions<AtomHandleList>())
  ;
}

