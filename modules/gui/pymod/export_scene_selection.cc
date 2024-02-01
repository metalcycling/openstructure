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
using namespace boost::python;

#include <ost/gui/scene_selection.hh>
#include <ost/gfx/gfx_node.hh>
#include <ost/gfx/entity.hh>
#include <ost/gfx/gfx_object.hh>
using namespace ost;
using namespace ost::gui;

void export_SceneSelection()
{
  class_<SceneSelection, boost::noncopyable>("SceneSelection",no_init)
    .def("Instance", &SceneSelection::Instance,
     return_value_policy<reference_existing_object>()).staticmethod("Instance")
    .def("GetActiveNodeCount",&SceneSelection::GetActiveNodeCount)
    .def("GetActiveNode",&SceneSelection::GetActiveNode)
    .def("GetViewEntity", &SceneSelection::GetViewEntity)
    .def("GetActiveViewCount", &SceneSelection::GetActiveViewCount)
    .def("GetActiveView", &SceneSelection::GetActiveView)
    .def("GetViewUnion", &SceneSelection::GetViewUnion)
    ;
}
