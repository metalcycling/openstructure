//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2023 by the OpenStructure authors
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
using namespace boost::python;

#include <ost/mol/alg/biounit.hh>

namespace{

  PyObject* wrap_to_bytes(ost::mol::alg::BUInfo& bu) {
    String str = bu.ToString();
    return PyBytes_FromStringAndSize(str.c_str(), str.size());
  }
  
  ost::mol::alg::BUInfo wrap_from_bytes(boost::python::object obj) {
    String str(PyBytes_AsString(obj.ptr()), PyBytes_Size(obj.ptr()));
    return ost::mol::alg::BUInfo::FromString(str);
  }

  list wrap_get_au_chains(const ost::mol::alg::BUInfo& info) {
    list return_list;
    const std::vector<std::vector<String> >& au_chains = info.GetAUChains();
    for(uint i = 0; i < au_chains.size(); ++i) {
      list tmp;
      for(uint j = 0; j < au_chains[i].size(); ++j) {
        tmp.append(au_chains[i][j]);
      }
      return_list.append(tmp);
    }
    return return_list;
  }

  list wrap_get_transformations(const ost::mol::alg::BUInfo& info) {
    list return_list;
    const std::vector<std::vector<geom::Mat4> >& tfs = info.GetTransformations();
    for(uint i = 0; i < tfs.size(); ++i) {
      list tmp;
      for(uint j = 0; j < tfs[i].size(); ++j) {
        tmp.append(tfs[i][j]);
      }
      return_list.append(tmp);
    }
    return return_list;
  }

  ost::mol::EntityHandle wrap_CreateBU_one(const ost::mol::EntityHandle& asu,
                                           const ost::io::MMCifInfoBioUnit& bu) {
    return ost::mol::alg::CreateBU(asu, bu);
  }

  ost::mol::EntityHandle wrap_CreateBU_two(const ost::mol::EntityHandle& asu,
                                           const ost::mol::alg::BUInfo& bu) {
    return ost::mol::alg::CreateBU(asu, bu);
  }

} // anon ns

void export_biounit() {

  class_<ost::mol::alg::BUInfo>("BUInfo", init<const ost::io::MMCifInfoBioUnit&>())
    .def("FromBytes", &wrap_from_bytes).staticmethod("FromBytes")
    .def("ToBytes", &wrap_to_bytes)
    .def("GetAUChains", &wrap_get_au_chains)
    .def("GetTransformations", &wrap_get_transformations)
  ;

  def("CreateBU", &wrap_CreateBU_one);
  def("CreateBU", &wrap_CreateBU_two);
}
