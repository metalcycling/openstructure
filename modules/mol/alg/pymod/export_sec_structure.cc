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

#include <ost/mol/alg/sec_struct.hh>
#include <ost/mol/alg/sec_structure_segments.hh>

using namespace boost::python;

namespace{

void AssignSecStructView(ost::mol::EntityView& view) {
  ost::mol::alg::AssignSecStruct(view);
}

void AssignSecStructHandle(ost::mol::EntityHandle& handle) {
  ost::mol::alg::AssignSecStruct(handle);
}

ost::mol::alg::SecStructureSegments 
ExtractHelicalSegments_handle(const ost::mol::ChainHandle& chain) {
  return ost::mol::alg::ExtractHelicalSegments(chain);
}
ost::mol::alg::SecStructureSegments 
ExtractHelicalSegments_view(const ost::mol::ChainView& chain) {
  return ost::mol::alg::ExtractHelicalSegments(chain);
}
ost::mol::alg::SecStructureSegments 
ExtractExtendedSegments_handle(const ost::mol::ChainHandle& chain) {
  return ost::mol::alg::ExtractExtendedSegments(chain);
}
ost::mol::alg::SecStructureSegments 
ExtractExtendedSegments_view(const ost::mol::ChainView& chain) {
  return ost::mol::alg::ExtractExtendedSegments(chain);
}
ost::mol::alg::SecStructureSegments 
ExtractSecStructureSegments_handle(const ost::mol::ChainHandle& chain) {
  return ost::mol::alg::ExtractSecStructureSegments(chain);
}
ost::mol::alg::SecStructureSegments 
ExtractSecStructureSegments_view(const ost::mol::ChainView& chain) {
  return ost::mol::alg::ExtractSecStructureSegments(chain);
}

} // ns

void export_sec_struct() {
  def("AssignSecStruct", &AssignSecStructView, (arg("ent")));
  def("AssignSecStruct", &AssignSecStructHandle, (arg("ent")));
}

void export_sec_struct_segments() {

  class_<ost::mol::alg::SecStructureSegment>("SecStructureSegment", init<int,int,ost::mol::SecStructure>())
    .def_readwrite("first", &ost::mol::alg::SecStructureSegment::first)
    .def_readwrite("last", &ost::mol::alg::SecStructureSegment::last)
    .def_readwrite("ss_type", &ost::mol::alg::SecStructureSegment::ss_type)
   ;

   class_<ost::mol::alg::SecStructureSegments>("SecStructureSegments", init<>())
     .def(vector_indexing_suite<ost::mol::alg::SecStructureSegments>())
   ;

   def("ExtractHelicalSegments", &ExtractHelicalSegments_handle, (arg("chain")));
   def("ExtractHelicalSegments", &ExtractHelicalSegments_view, (arg("chain")));
   def("ExtractExtendedSegments", &ExtractExtendedSegments_handle, (arg("chain")));
   def("ExtractExtendedSegments", &ExtractExtendedSegments_view, (arg("chain")));
   def("ExtractSecStructureSegments", &ExtractSecStructureSegments_handle, (arg("chain")));
   def("ExtractSecStructureSegments", &ExtractSecStructureSegments_view, (arg("chain")));
}

