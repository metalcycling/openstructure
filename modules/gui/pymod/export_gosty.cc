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
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;

#include <ost/gui/gosty_app.hh>

#include <ost/gui/gl_win.hh>
#include <ost/gui/perspective.hh>
#include <ost/gui/python_shell/python_shell.hh>
#include <ost/gui/scene_win/scene_win.hh>
#include <ost/gui/tools/tool_options_win.hh>
#include <ost/gui/scene_win/scene_win.hh>
#include <ost/gui/sequence_viewer/sequence_viewer.hh>
#include <ost/gui/messages/message_widget.hh>

#include "transfer_ownership.hh"
#include "sip_handler.hh"

#include "data_viewer_proxy.hh"
using namespace ost::img::gui;

using namespace ost;
using namespace ost::gui;

namespace {

DataViewerProxyPtr app_create_data_viewer1(GostyApp* app, const ost::img::ImageHandle& d, const QString& name)
{
  return DataViewerProxyPtr(new DataViewerProxy(app->CreateDataViewer(d,name)));
}
DataViewerProxyPtr app_create_data_viewer2(GostyApp* app, const ost::img::ImageHandle& d)
{
  return DataViewerProxyPtr(new DataViewerProxy(app->CreateDataViewer(d)));
}
DataViewerProxyPtr app_create_data_viewer3(GostyApp* app, const ost::img::ImageHandle& d, const QString& name, int noparent)
{
  return DataViewerProxyPtr(new DataViewerProxy(app->CreateDataViewer(d,name,noparent)));
}
}

void app_add_widget_to_app_a(GostyApp* app, const QString& ident,
                               QWidget* object)
{
  if(QWidget* widget = dynamic_cast<QWidget*>(object)){
    TransferOwnership(widget);
    app->AddWidgetToApp(ident, reinterpret_cast<QWidget*>(widget));
  }
}

void app_add_widget_to_app_b(GostyApp* app, const QString& ident,
                               object py_object)
{
  if(QWidget* widget = get_cpp_qobject<QWidget>(py_object)){
    app->AddWidgetToApp(ident, widget);
  }
}

object get_widget(GostyApp* app, const QString& ident)
{
  return get_py_qobject<QWidget>(app->GetWidget(ident));
}

void export_Gosty()
{
  class_<GostyApp, boost::noncopyable>("GostyApp", no_init)
    .def("Instance", &GostyApp::Instance,
         return_value_policy<reference_existing_object>()).staticmethod("Instance")
    .def("SetAppTitle", &GostyApp::SetAppTitle)
    .def("StopScript",&GostyApp::StopScript)
    .def("GetPyShell", &GostyApp::GetPyShell,
        return_value_policy<reference_existing_object>())
    .add_property("py_shell", make_function(&GostyApp::GetPyShell,
        return_value_policy<reference_existing_object>()))
    .def("GetGLWin", &GostyApp::GetGLWin,
        return_value_policy<reference_existing_object>())
    .add_property("gl_win", make_function(&GostyApp::GetGLWin,
        return_value_policy<reference_existing_object>()))
    .def("GetSceneWin", &GostyApp::GetSceneWin,
        return_value_policy<reference_existing_object>())
    .add_property("scene_win", make_function(&GostyApp::GetSceneWin,
        return_value_policy<reference_existing_object>()))
    .def("GetSequenceViewer", &GostyApp::GetSequenceViewer,
        return_value_policy<reference_existing_object>())
    .add_property("seq_viewer", make_function(&GostyApp::GetSequenceViewer,
        return_value_policy<reference_existing_object>()))
    .def("GetToolOptionsWin", &GostyApp::GetToolOptionsWin,
        return_value_policy<reference_existing_object>())
    .def("GetWidget", &get_widget)
    .add_property("tool_options_win", make_function(&GostyApp::GetToolOptionsWin,
        return_value_policy<reference_existing_object>()))
    .def("GetMessageWidget", &GostyApp::GetMessageWidget,
        return_value_policy<reference_existing_object>())
    .add_property("message_widget", make_function(&GostyApp::GetMessageWidget,
        return_value_policy<reference_existing_object>()))
    .def("CreateDataViewer", &app_create_data_viewer1)
    .def("CreateDataViewer", &app_create_data_viewer2)
    .def("CreateDataViewer", &app_create_data_viewer3)
    .def("ProcessEvents", &GostyApp::ProcessEvents)           
    .add_property("perspective", 
                  make_function(&GostyApp::GetPerspective, 
                                return_value_policy<reference_existing_object>()))
    .def("AddWidgetToApp", &app_add_widget_to_app_a)
    .def("AddWidgetToApp", &app_add_widget_to_app_b)
    .def("GetPerspective", &GostyApp::GetPerspective, 
         return_value_policy<reference_existing_object>())
    .def("TryStereo",&GostyApp::TryStereo)
  ;
  register_ptr_to_python<GostyApp*>();
}
