//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2021 by the OpenStructure authors
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

#include <ost/io/mol/omf.hh>

using namespace ost;
using namespace ost::io;

namespace{

  PyObject* wrap_to_bytes(OMFPtr omf) {
    String str = omf->ToString();
    return PyBytes_FromStringAndSize(str.c_str(), str.size());
  }
  
  OMFPtr wrap_from_bytes(boost::python::object obj) {
    String str(PyBytes_AsString(obj.ptr()), PyBytes_Size(obj.ptr()));
    return OMF::FromString(str);
  }

}

void export_omf_io() {
  class_<OMF, OMFPtr>("OMF",no_init)
    .def("FromEntity", &OMF::FromEntity).staticmethod("FromEntity")
    .def("FromMMCIF", &OMF::FromMMCIF).staticmethod("FromMMCIF")
    .def("FromFile", &OMF::FromFile).staticmethod("FromFile")
    .def("FromBytes", &wrap_from_bytes).staticmethod("FromBytes")
    .def("ToFile", &OMF::ToFile)
    .def("ToBytes", &wrap_to_bytes)
    .def("GetAU", &OMF::GetAU)
    .def("GetAUChain", &OMF::GetAUChain)
    .def("GetBU", &OMF::GetBU)
  ;
}
