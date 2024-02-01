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
#include <vector>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <ost/gui/messages/message_widget.hh>

#include "sip_handler.hh"

using namespace boost::python;
using namespace ost;
using namespace ost::gui;

namespace {

void log_message_a(MessageWidget* message_widget, const QString& message, QMessageBox::Icon icon = QMessageBox::Information)
{
  message_widget->LogMessage(message,icon);
}

void log_message_b(MessageWidget * message_widget, object py_object)
{
  if(py_object.ptr() != Py_None){
    if(PyObject_IsInstance(py_object.ptr(), (PyObject*)&PyUnicode_Type)){
      String message = extract < std::string > (py_object);
      message_widget->LogMessage(QString(message.c_str()));
    }
    else if(QStandardItem* item = get_cpp_qobject<QStandardItem>(py_object)){
      message_widget->LogMessage(item);
    }
  }
}

}

void export_MessageWidget()
{
  enum_<QMessageBox::Icon>("MessageIcon")
      .value("QUESTION", QMessageBox::Question)
      .value("INFORMATION", QMessageBox::Information)
      .value("WARNING", QMessageBox::Warning)
      .value("CRITICAL", QMessageBox::Critical)
      .export_values()
      ;

  class_<MessageWidget, boost::noncopyable>("MessageWidget", no_init)
     .def("Show", &MessageWidget::show)
     .def("Hide", &MessageWidget::hide)
     .def("Clear", &MessageWidget::Clear)
     .def("RemoveSelected", &MessageWidget::RemoveSelected)
     .def("LogMessage", &log_message_a)
     .def("LogMessage", &log_message_b)
     .def("GetQObject",&get_py_qobject<MessageWidget>)
     .add_property("qobject", &get_py_qobject<MessageWidget>)
   ;
}

