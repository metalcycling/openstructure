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
#include <string>
#include <boost/python/register_ptr_to_python.hpp>


#include <ost/gui/widget.hh>

#include "sip_handler.hh"

using namespace boost::python;
using namespace ost;
using namespace ost::gui;


struct WrappedWidget : public Widget
{
  WrappedWidget(PyObject *p, object py_object):
       Widget(NULL,NULL){
    if (QWidget* widget=get_cpp_qobject<QWidget>(py_object)) {
      this->SetInternalWidget(widget);
    }
  }
  virtual bool Restore(const QString& prefix){return true;}

  virtual bool Save(const QString& prefix){return true;}

  virtual ~WrappedWidget(){ }
private:
  QString   unique_id_;
};

void export_Widget()
{
  class_<Widget, WrappedWidget, boost::noncopyable>("Widget", init<object>())
    .def("Save", &WrappedWidget::Save)
    .def("Restore", &WrappedWidget::Restore)
    .def("SetDestroyOnClose", &WrappedWidget::SetDestroyOnClose)
    .def("DestroyOnClose", &WrappedWidget::DestroyOnClose)
    .def("GetQObject",&get_py_qobject<Widget>)
    .add_property("qobject", &get_py_qobject<Widget>)
  ;
}

