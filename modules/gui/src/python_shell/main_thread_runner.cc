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
/*
  MainThreadRunner

  Author: Andreas Schenk
*/

#include <boost/python.hpp>
#include "main_thread_runner.hh"

namespace ost { namespace gui {
namespace bp = boost::python;

MainThreadRunner::MainThreadRunner():
  QObject()
{
}
void MainThreadRunner::Run(const QString& widget)
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  bp::object  main_module = bp::import("__main__");
  bp::object main_namespace = main_module.attr("__dict__");
  std::string cmd = "Proxy._widget_=QtGui." + widget.toStdString() + "()\n";
  bp::exec(cmd.c_str(), main_namespace, main_namespace);
  PyGILState_Release(gstate);
}

}}//ns
