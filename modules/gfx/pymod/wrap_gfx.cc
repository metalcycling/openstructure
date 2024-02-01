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

#include <ost/gfx/module_config.hh>
#if OST_SHADER_SUPPORT_ENABLED
#include <ost/gfx/shader.hh>
#endif
#include <ost/info/info.hh>
#include <ost/gfx/prim_list.hh>
#include <ost/gfx/gfx_test_object.hh>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

extern void export_Scene();
extern void export_GfxNode();
extern void export_GfxObj();
extern void export_Entity();
extern void export_Surface();
extern void export_primlist();
extern void export_primitives();
extern void export_color();
extern void export_gradient();
extern void export_bitmap();
extern void export_Exporter();
extern void export_Map();
extern void export_SymmetryNode();
extern void export_SceneObserver();

extern void export_RenderOptions();

extern void export_ColorOps();

extern void export_GLWinBase();

using namespace ost;
using namespace ost::gfx;

BOOST_PYTHON_MODULE(_ost_gfx)
{
  export_Scene();
  export_GfxNode();
  export_GfxObj();
  export_Entity();
  export_Surface();
  export_RenderOptions();
  export_Map();
  export_SymmetryNode();
  export_SceneObserver();
  export_ColorOps();
  export_GLWinBase();
  export_primitives();
  export_primlist();
  export_color();
  export_gradient();
  export_bitmap();

  enum_<RenderMode::Type>("RenderMode")
    .value("SIMPLE",RenderMode::SIMPLE)
    .value("DEBUG",RenderMode::DEBUG)
    .value("DEBUG2",RenderMode::DEBUG2)
    .value("TEST",RenderMode::TEST)
    .value("TEST2",RenderMode::TEST2)
    .value("CUSTOM",RenderMode::CUSTOM)
    .value("CPK",RenderMode::CPK)
    .value("LINE_TRACE",RenderMode::LINE_TRACE)
    .value("TRACE",RenderMode::TRACE)    
    .value("SLINE",RenderMode::SLINE)
    .value("HSC",RenderMode::HSC)
    .value("TUBE",RenderMode::TUBE)
    .value("FILL",RenderMode::FILL)
    .value("OUTLINE",RenderMode::OUTLINE)
    .export_values()
    ;
  enum_<InputCommand>("InputCommand")
    .value("INPUT_COMMAND_NONE", INPUT_COMMAND_NONE)
    .value("INPUT_COMMAND_ROTX", INPUT_COMMAND_ROTX)
    .value("INPUT_COMMAND_ROTY", INPUT_COMMAND_ROTY)
    .value("INPUT_COMMAND_ROTZ", INPUT_COMMAND_ROTZ)
    .value("INPUT_COMMAND_TRANSX", INPUT_COMMAND_TRANSX)
    .value("INPUT_COMMAND_TRANSY", INPUT_COMMAND_TRANSY)
    .value("INPUT_COMMAND_TRANSZ", INPUT_COMMAND_TRANSZ)
    .value("INPUT_COMMAND_SLABN", INPUT_COMMAND_SLABN)
    .value("INPUT_COMMAND_SLABF", INPUT_COMMAND_SLABF)
    .value("INPUT_COMMAND_AUTOSLAB", INPUT_COMMAND_AUTOSLAB)
    .value("INPUT_COMMAND_TOGGLE_FOG", INPUT_COMMAND_TOGGLE_FOG)
    .value("INPUT_COMMAND_CUSTOM1", INPUT_COMMAND_CUSTOM1)
    .value("INPUT_COMMAND_CUSTOM2", INPUT_COMMAND_CUSTOM2)
    .value("INPUT_COMMAND_CUSTOM3", INPUT_COMMAND_CUSTOM3)
    .value("INPUT_COMMAND_REBUILD", INPUT_COMMAND_REBUILD)
    .export_values()
    ;
  enum_<InputDevice>("InputDevice")
    .value("INPUT_DEVICE_MOUSE", INPUT_DEVICE_MOUSE)
    .export_values()
  ;
  enum_<TransformTarget>("TransformTarget")
   .value("TRANSFORM_VIEW", TRANSFORM_VIEW)
   .value("TRANSFORM_RIGID", TRANSFORM_RIGID)
   .value("TRANSFORM_TORSION", TRANSFORM_TORSION)
   .value("TRANSFORM_ANGLE", TRANSFORM_ANGLE)
   .value("TRANSFORM_ROTAMER", TRANSFORM_ROTAMER)
   .export_values()
  ;
  class_<InputEvent>("InputEvent", init<InputDevice, InputCommand, float>())
   .def(init<InputDevice,InputCommand,int,int,float>())
  ;

  class_<GfxTestObj, bases<GfxObj>, boost::noncopyable>("GfxTestObj", init<>());
  
#if OST_SHADER_SUPPORT_ENABLED
  class_<Shader, boost::noncopyable>("Shader", no_init)
    .def("Instance",&Shader::Instance,
         return_value_policy<reference_existing_object>()).staticmethod("Instance")
    .def("PushProgram",&Shader::PushProgram)
    .def("PopProgram",&Shader::PopProgram)
    .def("Activate",&Shader::Activate)
    ;
#endif

  export_primitives();
  export_Exporter();
}



