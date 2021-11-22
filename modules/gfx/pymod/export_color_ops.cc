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

#include <ost/config.hh>

#include <ost/gfx/color_ops/color_op.hh>
#include <ost/gfx/color_ops/basic_gradient_color_op.hh>
#include <ost/gfx/color_ops/by_element_color_op.hh>
#include <ost/gfx/color_ops/by_chain_color_op.hh>
#include <ost/gfx/color_ops/uniform_color_op.hh>
#include <ost/gfx/color_ops/gradient_color_op.hh>
#include <ost/gfx/color_ops/gradient_level_color_op.hh>
#include <ost/gfx/color_ops/entity_view_color_op.hh>
#include <ost/gfx/color_ops/map_handle_color_op.hh>


using namespace boost::python;

using namespace ost;
using namespace ost::gfx;

void export_ColorOps()
{
  class_<ColorOp>("ColorOp", init<const String& >())
    .def(init<const ColorOp&>())
    .def(init<const mol::QueryViewWrapper&, int>())
    .def("CanApplyTo",&ColorOp::CanApplyTo)
    .def("GetName",&ColorOp::GetName,return_value_policy<copy_const_reference>())
    .def("SetSelection",&ColorOp::SetSelection)
    .def("GetSelection",&ColorOp::GetSelection)
    .def("SetSelectionFlags",&ColorOp::SetSelectionFlags)
    .def("GetSelectionFlags",&ColorOp::GetSelectionFlags)
    .def("SetMask",&ColorOp::SetMask)
    .def("GetMask",&ColorOp::GetMask)
    .def("ToInfo",&ColorOp::ToInfo)
    .def("FromInfo",&ColorOp::FromInfo)
    .staticmethod("FromInfo")
  ;

  class_<BasicGradientColorOp, bases<ColorOp> >("BasicGradientColorOp", init<const String& , const gfx::Gradient& , mol::Prop::Level>())
    .def(init<const BasicGradientColorOp&>())
    .def("SetGradient",&BasicGradientColorOp::SetGradient)
    .def("GetGradient",&BasicGradientColorOp::GetGradient)
    .def("SetLevel",&BasicGradientColorOp::SetLevel)
    .def("GetLevel",&BasicGradientColorOp::GetLevel)
    .def("FromInfo",&BasicGradientColorOp::FromInfo)
    .staticmethod("FromInfo")
  ;

  class_<ByElementColorOp, bases<ColorOp> >("ByElementColorOp", init<>())
    .def(init<const ByElementColorOp&>())
    .def(init<const String&>())
    .def(init<const mol::QueryViewWrapper&>())
    .def(init<const String&, int>())
    .def(init<const mol::QueryViewWrapper&, int>())
    .def("FromInfo",&ByElementColorOp::FromInfo)
    .staticmethod("FromInfo")
  ;

  class_<ByChainColorOp, bases<ColorOp> >("ByChainColorOp", init<>())
    .def(init<const ByChainColorOp&>())
    .def(init<const String&>())
    .def(init<const mol::QueryViewWrapper&>())
    .def(init<const String&, int>())
    .def(init<const mol::QueryViewWrapper&, int>())
    .def("GetChainCount",&ByChainColorOp::GetChainCount)
    .def("SetChainCount",&ByChainColorOp::SetChainCount)
    .def("FromInfo",&ByChainColorOp::FromInfo)
    .staticmethod("FromInfo")
  ;

  class_<UniformColorOp, bases<ColorOp> >("UniformColorOp", init<>())
    .def(init<const UniformColorOp&>())
    .def(init<const String&, const gfx::Color&>())
    .def(init<const String&, int, const gfx::Color&>())
    .def(init<const mol::QueryViewWrapper&, const gfx::Color&>())
    .def(init<const mol::QueryViewWrapper&, int, const gfx::Color&>())
    .def("SetColor",&UniformColorOp::SetColor)
    .def("GetColor",&UniformColorOp::GetColor)
    .def("FromInfo",&UniformColorOp::FromInfo)
    .staticmethod("FromInfo")
  ;

  class_<GradientColorOp, bases<ColorOp> >("GradientColorOp", init<>())
    .def(init<const GradientColorOp&>())
    .def(init<const String&, const String&, const gfx::Gradient&, float, float>())
    .def(init<const String&, int, const String&, const gfx::Gradient&, float, float>())
    .def(init<const mol::QueryViewWrapper&, const String&, const gfx::Gradient&, float, float>())
    .def(init<const mol::QueryViewWrapper&, int, const String&, const gfx::Gradient&, float, float>())
    .def("SetProperty",&GradientColorOp::SetProperty)
    .def("GetProperty",&GradientColorOp::GetProperty)
    .def("SetGradient",&GradientColorOp::SetGradient)
    .def("GetGradient",&GradientColorOp::GetGradient)
    .def("GetCalculateMinMax",&GradientColorOp::GetCalculateMinMax)
    .def("SetMaxV",&GradientColorOp::SetMaxV)
    .def("GetMaxV",&GradientColorOp::GetMaxV)
    .def("SetMinV",&GradientColorOp::SetMinV)
    .def("GetMinV",&GradientColorOp::GetMinV)
    .def("FromInfo",&GradientColorOp::FromInfo)
    .staticmethod("FromInfo")
  ;

  class_<GradientLevelColorOp, bases<GradientColorOp> >("GradientLevelColorOp", init<>())
    .def(init<const GradientLevelColorOp&>())
    .def(init<const String&, const String&, const gfx::Gradient&, float, float, mol::Prop::Level>())
    .def(init<const String&, int, const String&, const gfx::Gradient&, float, float, mol::Prop::Level>())
    .def(init<const String&, const String&, const gfx::Gradient&, mol::Prop::Level>())
    .def(init<const String&, int, const String&, const gfx::Gradient&, mol::Prop::Level>())
    .def(init<const mol::QueryViewWrapper&, const String&, const gfx::Gradient&, float, float, mol::Prop::Level>())
    .def(init<const mol::QueryViewWrapper&, int, const String&, const gfx::Gradient&, float, float, mol::Prop::Level>())
    .def(init<const mol::QueryViewWrapper&, const String&, const gfx::Gradient&, mol::Prop::Level>())
    .def(init<const mol::QueryViewWrapper&, int, const String&, const gfx::Gradient&, mol::Prop::Level>())
    .def("SetLevel",&GradientLevelColorOp::SetLevel)
    .def("GetLevel",&GradientLevelColorOp::GetLevel)
    .def("FromInfo",&GradientLevelColorOp::FromInfo)
    .staticmethod("FromInfo")
  ;


  class_<EntityViewColorOp, bases<GradientColorOp> >("EntityViewColorOp", init<>())
    .def(init<const EntityViewColorOp&>())
    .def(init<const String&, const gfx::Gradient&, float, float, const mol::EntityView&>())
    .def(init<int, const String&, const gfx::Gradient&, float, float, const mol::EntityView&>())
    .def("SetEntityView",&EntityViewColorOp::SetEntityView)
    .def("GetEntityView",&EntityViewColorOp::GetEntityView,return_value_policy<copy_const_reference>())
    .def("FromInfo",&EntityViewColorOp::FromInfo)
    .staticmethod("FromInfo")
  ;

  class_<MapHandleColorOp, bases<GradientColorOp> >("MapHandleColorOp", init<>())
    .def(init<const MapHandleColorOp&>())
    .def(init<const String&, const String&, const gfx::Gradient&, float, float, const img::MapHandle&>())
    .def(init<const String&, int, const String&, const gfx::Gradient&, float, float, const img::MapHandle&>())
    .def(init<const mol::QueryViewWrapper&, const String&, const gfx::Gradient&, float, float, const img::MapHandle&>())
    .def(init<const mol::QueryViewWrapper&, int, const String&, const gfx::Gradient&, float, float, const img::MapHandle&>())
    .def("SetMapHandle",&MapHandleColorOp::SetMapHandle)
    .def("GetMapHandle",&MapHandleColorOp::GetMapHandle,return_value_policy<copy_const_reference>())
    .def("FromInfo",&MapHandleColorOp::FromInfo)
    .staticmethod("FromInfo")
  ;
}



