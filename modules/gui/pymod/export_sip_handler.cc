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

#include <QWidget>

#include "sip_handler.hh"


using namespace boost::python;
using namespace ost;
using namespace ost::gui;

namespace {

unsigned long sip_handle(QObject* w)
{
  return reinterpret_cast<unsigned long>(w);
}

QWidget* bpqt_handle(unsigned long a)
{
  return reinterpret_cast<QWidget*>(a);
}

}

void export_SipHandler()
{
  class_<SipHandlerBase, boost::noncopyable>("SipHandlerBase", no_init)
    .def("GetSipHandle", &SipHandlerBase::GetSipHandle)
    .add_property("sip_handle", &SipHandlerBase::GetSipHandle)
  ;
  def("SipHandle", &sip_handle);
  def("BPQtHandle", &bpqt_handle,
      return_value_policy<return_opaque_pointer>());  
}

