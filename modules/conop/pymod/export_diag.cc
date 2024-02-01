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
//
#include <boost/python.hpp>
#include <ost/conop/diag.hh>
using namespace ost;

using namespace boost::python;
using namespace ost::conop;


void export_diag() {

  enum_<DiagType>("DiagType")
    .value("DIAG_UNK_ATOM", DIAG_UNK_ATOM)
    .value("DIAG_UNK_RESIDUE", DIAG_UNK_RESIDUE)
    .value("DIAG_MISSING_ATOM", DIAG_MISSING_ATOM)
    .value("DIAG_NONSTD_RESIDUE", DIAG_NONSTD_RESIDUE)
    .export_values()
  ;

  class_<Diag>("Diag", init<DiagType,const char*>())
    .def("Format", &Diag::Format, arg("colored")=true)
  ;
  class_<Diagnostics, DiagnosticsPtr>("Diagnostics")
  ;
}
