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

#include <ost/gfx/color_ops/map_handle_color_op.hh>

namespace ost { namespace gfx {

MapHandleColorOp::MapHandleColorOp() : GradientColorOp(), mh_(){
  this->Init();
}

MapHandleColorOp::MapHandleColorOp(const String& selection, const String& property, const gfx::Gradient& gradient, float minv, float maxv, const img::MapHandle& mh) :
  GradientColorOp(selection, property, gradient, minv, maxv), mh_(mh){
  this->Init();
}

MapHandleColorOp::MapHandleColorOp(const String& selection, int mask, const String& property, const gfx::Gradient& gradient, float minv, float maxv, const img::MapHandle& mh) :
  GradientColorOp(selection, mask, property, gradient, minv, maxv), mh_(mh){
  this->Init();
}

MapHandleColorOp::MapHandleColorOp(const mol::QueryViewWrapper& query_view, const String& property, const gfx::Gradient& gradient, float minv, float maxv, const img::MapHandle& mh) :
  GradientColorOp(query_view, property, gradient, minv, maxv), mh_(mh){
  this->Init();
}

MapHandleColorOp::MapHandleColorOp(const mol::QueryViewWrapper& query_view, int mask, const String& property, const gfx::Gradient& gradient, float minv, float maxv, const img::MapHandle& mh) :
  GradientColorOp(query_view, mask, property, gradient, minv, maxv), mh_(mh){
  this->Init();
}

void MapHandleColorOp::Init()
{
  this->SetName("MapHandle gradient");
}

bool MapHandleColorOp::CanApplyTo(const GfxObjP& obj) const{
  if(dynamic_cast<Entity*>(obj.get())){
    return true;
  }
  else if(dynamic_cast<Surface*>(obj.get())){
    return true;
  }
  return false;
}

void MapHandleColorOp::ApplyTo(GfxObjP& objP) const{
  GfxObj* obj = objP.get();
  if(Entity* ent=dynamic_cast<Entity*>(obj)){
    ent->Apply(*this, false);
  }
  else if(Surface* sur=dynamic_cast<Surface*>(obj)){
    sur->Apply(*this, false);
  }
}

void MapHandleColorOp::SetMapHandle(const img::MapHandle& mh){
  mh_ = mh;
}

const img::MapHandle& MapHandleColorOp::GetMapHandle() const{
  return mh_;
}

void MapHandleColorOp::ToInfo(info::InfoGroup& group) const{
  info::InfoGroup super_group = group.CreateGroup("GradientColorOp");
  GradientColorOp::ToInfo(super_group);
  //TODO persist map handle
  std::ostringstream ss;
  ss << "";
  group.SetTextData(ss.str());
}

gfx::MapHandleColorOp MapHandleColorOp::FromInfo(info::InfoGroup& group){
  info::InfoGroup super_group = group.GetGroup("GradientColorOp");
  gfx::GradientColorOp gop = GradientColorOp::FromInfo(super_group);
  std::istringstream ss(group.GetTextData());
  mol::QueryViewWrapper wrapper(gop.GetSelection(),gop.GetSelectionFlags());
  gfx::Gradient gradient = gop.GetGradient();
  int mask = gop.GetMask();
  String property = gop.GetProperty();
  float minv = gop.GetMinV();
  float maxv = gop.GetMaxV();
  //TODO load map handle
  return gfx::MapHandleColorOp(wrapper, mask, property, gradient, minv, maxv, img::MapHandle());
}

}}
