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

#include <ost/gui/scene_win/scene_win.hh>
#include <ost/gui/scene_win/context_menu.hh>

#include "sip_handler.hh"

#include <QAction>

using namespace boost::python;
using namespace ost;
using namespace ost::gui;

namespace {

void context_menu_add_action(ContextMenu* cm, object action, int types)
{
  QAction* act = get_cpp_qobject<QAction>(action);
  if(act){
    ContextActionTypes type = ContextActionTypes(types);
    cm->AddAction(act,type);
  }
}

}

void export_SceneWin()
{
  enum_<ContextActionType>("ContextActionType")
      .value("GFX_OBJECT", GFX_OBJECT)
      .value("ENTITY", ENTITY)
      .value("ENTITY_VIEW", ENTITY_VIEW)
      .value("CUSTOM_VIEW", CUSTOM_VIEW)
      .value("NOT_VISIBLE", NOT_VISIBLE)
      .value("NOT_HIDDEN", NOT_HIDDEN)
      .value("NOT_SCENE", NOT_SCENE)
      .value("VIEWS_SAME_OBJECT", VIEWS_SAME_OBJECT)
      .value("SINGLE", SINGLE)
      .value("MULTI", MULTI)
      .value("MAP", MAP)
      .export_values()
      ;

  class_<ContextMenu, boost::noncopyable>("ContextMenu", no_init)
    .def("AddAction", &context_menu_add_action)
    .def("GetQObject",&get_py_qobject<ContextMenu>)
    .add_property("qobject", &get_py_qobject<ContextMenu>)
  ;

  class_<SceneWin, boost::noncopyable>("SceneWin", no_init)
    .def("Show", &SceneWin::show)
    .def("Hide", &SceneWin::hide)
    .def("GetContextMenu", &SceneWin::GetContextMenu, return_value_policy<reference_existing_object>())
    .def("GetQObject",&get_py_qobject<SceneWin>)
    .add_property("qobject", &get_py_qobject<SceneWin>)
  ;

}

