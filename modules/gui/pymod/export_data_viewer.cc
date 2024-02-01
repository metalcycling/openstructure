//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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
  Authors: Ansgar Philippsen, Andreas Schenk
*/

#include <boost/python.hpp>
#include <ost/base.hh>
#include <ost/gui/data_viewer/data_viewer.hh>
#include "data_viewer_proxy.hh"
#include <ost/gui/data_viewer/overlay_manager.hh>

using namespace boost::python;
using namespace ost::img;
using namespace ost;
using namespace ost::gui;
using namespace ost::img::gui;

namespace {

void add_dock1(DataViewerProxy* p, QWidget* w, const String& name) 
{
  p->AddDockWidget(w,QString::fromStdString(name));
}

void add_dock2(DataViewerProxy* p, QWidget* w, const String& name, bool s)
{
  p->AddDockWidget(w,QString::fromStdString(name),s);
}

void add_dock3(DataViewerProxy* p, QWidget* w, const String& name, bool s, int area)
{
  p->AddDockWidget(w,QString::fromStdString(name),s,area);
}

} // anon ns

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(o_AddOverlay, AddOverlay, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(o_AddDockWidget, AddDockWidget, 2, 4)

void export_data_viewer()
{
  class_<DataViewerProxy, DataViewerProxyPtr >("DataViewer",no_init)
    .def("SetData",&DataViewerProxy::SetData)
    .def("SetName",&DataViewerProxy::SetName)
    .def("GetOverlayManager",&DataViewerProxy::GetOverlayManager)
    .def("GetNormalizer",&DataViewerProxy::GetNormalizer,return_value_policy<return_by_value>())
    .def("Renormalize",&DataViewerProxy::Renormalize)
    .def("AddOverlay",&DataViewerProxy::AddOverlay,o_AddOverlay())
    .def("ClearOverlays",&DataViewerProxy::ClearOverlays)
    .def("ClearActiveOverlay",&DataViewerProxy::ClearActiveOverlay)
    .def("GetSelection",&DataViewerProxy::GetSelection)
    .def("SetSelection",&DataViewerProxy::SetSelection)
    .def("UpdateView",&DataViewerProxy::UpdateView)
    .def("Recenter",&DataViewerProxy::Recenter)
    .def("AddDockWidget",&DataViewerProxy::AddDockWidget,o_AddDockWidget())
    .def("AddDockWidget",add_dock1)
    .def("AddDockWidget",add_dock2)
    .def("AddDockWidget",add_dock3)
    .def("RemoveDockWidget",&DataViewerProxy::RemoveDockWidget)
    .def("SetAntialiasing",&DataViewerProxy::SetAntialiasing)
    .def("Show",&DataViewerProxy::Show)
    .def("Hide", &DataViewerProxy::Hide)
    .def("SetSlab", &DataViewerProxy::SetSlab)
    .def("GetSlab", &DataViewerProxy::GetSlab)
    .add_property("slab", &DataViewerProxy::GetSlab, &DataViewerProxy::SetSlab)
    .def("SetZoomScale", &DataViewerProxy::SetZoomScale)
    .def("GetZoomScale", &DataViewerProxy::GetZoomScale)
    .add_property("zoomscale", &DataViewerProxy::GetZoomScale, &DataViewerProxy::SetZoomScale)
    .def("SetViewerMin", &DataViewerProxy::SetViewerMin)
    .def("GetViewerMin", &DataViewerProxy::GetViewerMin)
    .add_property("viewer_min", &DataViewerProxy::GetViewerMin, &DataViewerProxy::SetViewerMin)
    .def("SetViewerMax", &DataViewerProxy::SetViewerMax)
    .def("GetViewerMax", &DataViewerProxy::GetViewerMax)
    .add_property("viewer_max", &DataViewerProxy::GetViewerMax, &DataViewerProxy::SetViewerMax)
    .def("SetGamma", &DataViewerProxy::SetGamma)
    .def("GetGamma", &DataViewerProxy::GetGamma)
    .add_property("gamma", &DataViewerProxy::GetGamma, &DataViewerProxy::SetGamma)
    .def("SetInvert", &DataViewerProxy::SetInvert)
    .def("GetInvert", &DataViewerProxy::GetInvert)
    .add_property("invert", &DataViewerProxy::GetInvert, &DataViewerProxy::SetInvert)
    .def("SetOffset", &DataViewerProxy::SetOffset)
    .def("GetOffset", &DataViewerProxy::GetOffset)
    .add_property("offset", &DataViewerProxy::GetOffset, &DataViewerProxy::SetOffset)
    .def("GetButtonGroup", &DataViewerProxy::GetButtonGroup)
    .def("GetQObject",&DataViewerProxy::GetQObject)
    .add_property("qobject",&DataViewerProxy::GetQObject)

      /*
    .def("SetAutoNormalize",&DataViewerProxy::SetAutoNormalize)
    .def("ShowClickedPosition",&DataViewerProxy::ShowClickedPosition,o_ShowClickedPosition())
    .def("GetClickedPosition",&DataViewerProxy::GetClickedPosition)
    */    
    ;
}
