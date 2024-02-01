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

#include <ost/gui/python_shell/python_shell.hh>
#include <ost/gui/python_shell/python_shell_widget.hh>

#include "sip_handler.hh"

using namespace boost::python;
using namespace ost;
using namespace ost::gui;

namespace {

void pyshell_set_tab_width(PythonShell* py_shell, int width) {
  py_shell->PyShell()->SetTabWidth(width);
}

int pyshell_get_tab_width(PythonShell* py_shell) {
  return py_shell->PyShell()->GetTabWidth();
}

}

void export_PyShell()
{ 
  class_<PythonShell, boost::noncopyable>("PythonShell", no_init)
    .def("Show", &PythonShell::show)
    .def("Hide", &PythonShell::hide)
    .def("SetTabWidth", &pyshell_set_tab_width)
    .def("GetTabWidth", &pyshell_get_tab_width)
    .def("GetQObject",&get_py_qobject<PythonShell>)
    .add_property("qobject", &get_py_qobject<PythonShell>)
  ;
}
