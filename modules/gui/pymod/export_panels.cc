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

#include <ost/gui/panels/panel_manager.hh>
#include <ost/gui/panels/panel_bar.hh>
#include "ost/gui/widget.hh"

#include "sip_handler.hh"



using namespace boost::python;
using namespace ost;
using namespace ost::gui;

namespace{
void panels_add_widget_a(PanelManager * panels, PanelPosition pos, const SipHandlerBase& sh, bool hidden=false)
{
  if(Widget* widget = reinterpret_cast<Widget*>(sh.GetSipHandle())){
    panels->AddWidget(pos,widget,hidden);
  }
}

void panels_add_widget_b(PanelManager * panels, PanelPosition pos, Widget* widget, bool hidden=false)
{
  panels->AddWidget(pos,widget,hidden);
}

void panels_add_widget_c(PanelManager * panels, PanelPosition pos, object py_object, bool hidden=false)
{
  if(Widget* widget = get_cpp_qobject<Widget>(py_object)){
    panels->AddWidget(pos,widget,hidden);
  }
}

void panels_remove_widget_a(PanelManager * panels, Widget* w)
{
  panels->RemoveWidget(w);
}

void panels_remove_widget_b(PanelManager * panels, const SipHandlerBase& sh)
{
  if(Widget* widget = reinterpret_cast<Widget*>(sh.GetSipHandle())){
    panels->RemoveWidget(widget);
  }
}

void panels_remove_widget_c(PanelManager * panels, object py_object)
{
  if(Widget* widget = get_cpp_qobject<Widget>(py_object)){
    panels->RemoveWidget(widget);
  }
}

void panels_add_widget_to_pool_a(PanelManager * panels, const QString& full_name, Widget* widget)
{
  panels->AddWidgetToPool(full_name,widget);
}

void panels_add_widget_to_pool_b(PanelManager * panels, const QString& full_name, const SipHandlerBase& sh)
{
  if(Widget* widget = reinterpret_cast<Widget*>(sh.GetSipHandle())){
    panels->AddWidgetToPool(full_name,widget);
  }
}

void panels_add_widget_to_pool_c(PanelManager * panels, const QString& full_name, object py_object)
{
  if(Widget* widget = get_cpp_qobject<Widget>(py_object)){
    panels->AddWidgetToPool(full_name,widget);
  }
}

void panels_add_widget_to_pool_d(PanelManager * panels, const QString& name, int limit=-1)
{
  panels->AddWidgetToPool(name,limit);
}

object panels_get_menu(PanelManager* panels)
{
  return get_py_qobject<QMenu>(panels->GetMenu());
}

void panel_bar_add_widget_a(PanelBar * pb, Widget* widget, bool hidden=false)
{
  pb->AddWidget(widget,hidden);
}

void panel_bar_add_widget_b(PanelBar * pb, const SipHandlerBase& sh, bool hidden=false)
{
  if(Widget* widget = reinterpret_cast<Widget*>(sh.GetSipHandle())){
    pb->AddWidget(widget,hidden);
  }
}

void panel_bar_add_widget_c(PanelBar * pb, object py_object, bool hidden=false)
{
  if(Widget* widget = get_cpp_qobject<Widget>(py_object)){
    pb->AddWidget(widget,hidden);
  }
}

}

void export_Panels()
{
  enum_<PanelPosition>("PanelPosition")
      .value("LEFT_PANEL", LEFT_PANEL)
      .value("RIGHT_PANEL", RIGHT_PANEL)
      .value("BOTTOM_PANEL", BOTTOM_PANEL)
      .value("NONE", NONE)
      .export_values()
      ;

  class_<PanelManager, boost::noncopyable>("PanelManager", no_init)
    .def("AddWidget", &panels_add_widget_a, arg("hidden")=false)
    .def("AddWidget", &panels_add_widget_b, arg("hidden")=false)
    .def("AddWidget", &panels_add_widget_c, arg("hidden")=false)
    .def("AddWidgetByName", &PanelManager::AddWidgetByName)
    .def("RemoveWidget", &panels_remove_widget_a)
    .def("RemoveWidget", &panels_remove_widget_b)
    .def("RemoveWidget", &panels_remove_widget_c)
    .def("GetMenu", &panels_get_menu)
    .def("AddWidgetToPool", &panels_add_widget_to_pool_a)
    .def("AddWidgetToPool", &panels_add_widget_to_pool_b)
    .def("AddWidgetToPool", &panels_add_widget_to_pool_c)
    .def("AddWidgetToPool", &panels_add_widget_to_pool_d)
    .def("Save", &PanelManager::Save)
    .def("Restore", &PanelManager::Restore)
    .add_property("menu", &panels_get_menu)
    .def("GetQObject",&get_py_qobject<PanelManager>)
    .add_property("qobject", &get_py_qobject<PanelManager>)
  ;

  class_<PanelBar, boost::noncopyable>("PanelBar", no_init)
    .def("AddWidget", &panel_bar_add_widget_a, arg("hidden")=false)
    .def("AddWidget", &panel_bar_add_widget_b, arg("hidden")=false)
    .def("AddWidget", &panel_bar_add_widget_c, arg("hidden")=false)
    .def("Restore", &PanelBar::Restore)
    .def("Save", &PanelBar::Save)
    .def("Hide", &PanelBar::hide)
    .def("Show", &PanelBar::show)
    .def("GetQObject",&get_py_qobject<PanelBar>)
    .add_property("qobject", &get_py_qobject<PanelBar>)
  ;

}

