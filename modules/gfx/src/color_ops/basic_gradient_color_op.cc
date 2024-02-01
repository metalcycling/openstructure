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

#include "basic_gradient_color_op.hh"

#include <ost/gfx/surface.hh>
#include <ost/dyn_cast.hh>


namespace ost { namespace gfx {

BasicGradientColorOp::BasicGradientColorOp() : ColorOp(), gradient_(){
  this->Init();
}

BasicGradientColorOp::BasicGradientColorOp(const String& selection, const gfx::Gradient& gradient, mol::Prop::Level level) :
  ColorOp(selection), gradient_(gradient), level_(level){
  this->Init();
}

void BasicGradientColorOp::Init()
{
  this->SetName("Basic gradient");
}

bool BasicGradientColorOp::CanApplyTo(const GfxObjP& obj) const{
  if(dynamic_cast<Surface*>(obj.get()))
    return true;
  return false;
}

void BasicGradientColorOp::ApplyTo(GfxObjP& objP) const{
  GfxObj* obj = objP.get();
  if(Surface* sur=dynamic_cast<Surface*>(obj)){
    sur->Apply(*this, false);
  }
}

void BasicGradientColorOp::SetGradient(const gfx::Gradient& gradient){
  gradient_ = gradient;
}

gfx::Gradient BasicGradientColorOp::GetGradient() const{
  return gradient_;
}

void BasicGradientColorOp::SetLevel(mol::Prop::Level level){
  level_ = level;
}

mol::Prop::Level BasicGradientColorOp::GetLevel() const{
  return level_;
}

void BasicGradientColorOp::ToInfo(info::InfoGroup& group) const{
  info::InfoGroup super_group = group.CreateGroup("ColorOp");
  ColorOp::ToInfo(super_group);
  info::InfoGroup gradient_group = group.CreateGroup("Gradient");
  gradient_.GradientToInfo(gradient_group);
  std::ostringstream ss;
  ss << int(level_);
  group.SetTextData(ss.str());
}

gfx::BasicGradientColorOp BasicGradientColorOp::FromInfo(info::InfoGroup& group){
  info::InfoGroup super_group = group.GetGroup("ColorOp");
  gfx::ColorOp op = ColorOp::FromInfo(super_group);
  info::InfoGroup gradient_group = group.GetGroup("Gradient");
  gfx::Gradient gradient = gfx::Gradient::GradientFromInfo(gradient_group);
  std::istringstream ss(group.GetTextData());
  int l;
  ss >> l;
  mol::Prop::Level level = mol::Prop::Level(l);
  String selection = op.GetSelection();
  return gfx::BasicGradientColorOp(selection, gradient, level);
}

}}

