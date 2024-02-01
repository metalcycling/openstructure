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

#include <ost/gui/input/spnav_input.hh>

#include "sip_handler.hh"

using namespace boost::python;
using namespace ost;
using namespace ost::gui;

namespace {

object spnav_get_instance()
{
  static object pyqt5_module=import("PyQt5.QtCore");
  static object sip_module=import("sip");
  SpnavInput* si = SpnavInput::Instance();
  if(si->isValid()) {
    return ost::gui::get_py_qobject<SpnavInput>(si);
  } else {
    return object();
  }
}

}

void export_Input()
{
  class_<SpnavInput, boost::noncopyable >("SpnavInput",no_init)
    .def("GetQThread", &spnav_get_instance).staticmethod("GetQThread")
    .def("GetQObject", &get_py_qobject<SpnavInput>)
    .add_property("qobject", &get_py_qobject<SpnavInput>)
    .def("Instance", &SpnavInput::Instance,
     return_value_policy<reference_existing_object>()).staticmethod("Instance")
  ;
}

