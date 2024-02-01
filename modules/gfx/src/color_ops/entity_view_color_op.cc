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


#include <ost/gfx/entity.hh>
#include <ost/gfx/surface.hh>
#include <ost/dyn_cast.hh>

#include <ost/gfx/color_ops/entity_view_color_op.hh>

namespace ost { namespace gfx {

EntityViewColorOp::EntityViewColorOp() : GradientColorOp(), ev_(){
  this->Init();
}
EntityViewColorOp::EntityViewColorOp(int mask, const String& property, const gfx::Gradient& gradient, float minv, float maxv, const mol::EntityView& ev) :
  GradientColorOp("", mask, property, gradient, minv, maxv), ev_(ev){
  this->Init();
}

EntityViewColorOp::EntityViewColorOp(const String& property, const gfx::Gradient& gradient, float minv, float maxv, const mol::EntityView& ev) :
  GradientColorOp("", property, gradient, minv, maxv), ev_(ev){
  this->Init();
}

void EntityViewColorOp::Init()
{
  this->SetName("EntityView gradient");
}

bool EntityViewColorOp::CanApplyTo(const GfxObjP& obj) const{
  if(dynamic_cast<Entity*>(obj.get())){
    return true;
  }
  else if(dynamic_cast<Surface*>(obj.get())){
    return true;
  }
  return false;
}

void EntityViewColorOp::ApplyTo(GfxObjP& objP) const{
  GfxObj* obj = objP.get();
  if(Entity* ent=dynamic_cast<Entity*>(obj)){
    ent->Apply(*this, false);
  }
  else if(Surface* sur=dynamic_cast<Surface*>(obj)){
	sur->Apply(*this, false);
  }
}

void EntityViewColorOp::SetEntityView(const mol::EntityView& ev){
  ev_ = ev;
}

void EntityViewColorOp::SetSelection(const String& selection){
  //Do nothing
}

const mol::EntityView& EntityViewColorOp::GetEntityView() const{
  return ev_;
}

void EntityViewColorOp::ToInfo(info::InfoGroup& group) const{
  info::InfoGroup super_group = group.CreateGroup("GradientColorOp");
  GradientColorOp::ToInfo(super_group);
  //TODO persist entity view
  std::ostringstream ss;
  ss << "";
  group.SetTextData(ss.str());
}

gfx::EntityViewColorOp EntityViewColorOp::FromInfo(info::InfoGroup& group){
  info::InfoGroup super_group = group.GetGroup("GradientColorOp");
  gfx::GradientColorOp gop = GradientColorOp::FromInfo(super_group);
  std::istringstream ss(group.GetTextData());
  gfx::Gradient gradient = gop.GetGradient();
  int mask = gop.GetMask();
  String property = gop.GetProperty();
  float minv = gop.GetMinV();
  float maxv = gop.GetMaxV();
  //TODO load entity view
  return gfx::EntityViewColorOp(mask, property, gradient, minv, maxv, mol::EntityView());
}

}}

