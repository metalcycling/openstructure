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

namespace {
  ChainView get_chain(AtomView& a)
  {
    return a.GetResidue().GetChain();
  }
}

void export_AtomView()
{

  class_<AtomView, bases<AtomBase> >("AtomView", init<>())
    .def("GetResidue",&AtomView::GetResidue)
    .add_property("residue",&AtomView::GetResidue)
    .def("GetEntity",&AtomView::GetEntity)
    .add_property("entity",&AtomView::GetEntity)
    .def("GetChain",get_chain)
    .add_property("chain",get_chain)
    .def(self==self)
    .def(self!=self)
    .def("RemoveBonds", &AtomView::RemoveBonds)
    .def("GetHandle", &AtomView::GetHandle)
    .add_property("handle", &AtomView::GetHandle)
    .def("GetBondCount", &AtomView::GetBondCount)
    .def("GetBondList", &AtomView::GetBondList)
    .add_property("bonds", &AtomView::GetBondList)
    .def("GetHashCode", &AtomView::GetHashCode)
    .def("__hash__", &AtomView::GetHashCode)
    .add_property("hash_code", &AtomView::GetHashCode)
    .def("IsValid", &AtomView::IsValid)
    .add_property("valid", &AtomView::IsValid)
    .def("GetBondPartners", &AtomView::GetBondPartners)
  ;
  class_<AtomViewList>("AtomViewList", init<>())
    .def(vector_indexing_suite<AtomViewList>())
    .def(geom::VectorAdditions<AtomViewList>())
  ;
}

