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

#include <QMdiArea>
#include <QMenuBar>

#include <ost/gui/perspective.hh>
#include <ost/gui/main_area.hh>
#include <ost/gui/panels/panel_manager.hh>
#include <ost/gui/panels/button_bar.hh>
#include <ost/gui/panels/panel_bar.hh>

#include "sip_handler.hh"

using namespace boost::python;
using namespace ost;
using namespace ost::gui;

object app_get_menu_bar(Perspective* pers)
{
  return get_py_qobject<QMenuBar>(pers->GetMenuBar());
}

object app_get_menu(Perspective* pers, const QString& name)
{
  return get_py_qobject<QMenu>(pers->GetMenu(name));
}

void export_Perspective()
{
  class_<Perspective, boost::noncopyable>("Perspective", no_init)
    .def("GetPanels", &Perspective::GetPanels,
         return_value_policy<reference_existing_object>())
    .def("GetMainArea", &Perspective::GetMainArea,
         return_value_policy<reference_existing_object>())
    .def("GetMenuBar", &app_get_menu_bar)
    .def("GetMenu", &app_get_menu)
    .def("HideAllBars", &Perspective::HideAllBars)
    .def("StatusMessage", &Perspective::StatusMessage)
    .def("Restore", &Perspective::Restore)
    .def("Save", &Perspective::Save)
    .def("ChangeVerbositySlider", &Perspective::ChangeVerbositySlider)
    .add_property("main_area", make_function(&Perspective::GetMainArea,
                  return_value_policy<reference_existing_object>()))
    .add_property("panels", make_function(&Perspective::GetPanels,
                  return_value_policy<reference_existing_object>()))
    .add_property("menubar", make_function(&app_get_menu_bar))
  ;
}

