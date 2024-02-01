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

#include <ost/dyn_cast.hh>

#include <ost/gfx/entity.hh>
#include <ost/gfx/surface.hh>

#include <ost/gfx/color_ops/uniform_color_op.hh>

namespace ost { namespace gfx {

UniformColorOp::UniformColorOp() : ColorOp(), color_(){
  this->Init();
}

UniformColorOp::UniformColorOp(const String& selection, 
                               const gfx::Color& color): 
  ColorOp(selection), color_(color){
  this->Init();
}

UniformColorOp::UniformColorOp(const String& selection, int mask,
                               const gfx::Color& color):
  ColorOp(selection,mask), color_(color){
  this->Init();
}

UniformColorOp::UniformColorOp(const mol::QueryViewWrapper& query_view,
                               const gfx::Color& color):
  ColorOp(query_view), color_(color){
  this->Init();
}

UniformColorOp::UniformColorOp(const mol::QueryViewWrapper& query_view,
                               int mask, const gfx::Color& color):
  ColorOp(query_view,mask), color_(color){
  this->Init();
}

void UniformColorOp::Init()
{
  this->SetName("Uniform color");
}

bool UniformColorOp::CanApplyTo(const GfxObjP& obj) const
{
  return dynamic_cast<Entity*>(obj.get()) || dynamic_cast<Surface*>(obj.get());
}

void UniformColorOp::ApplyTo(GfxObjP& objP) const
{
  GfxObj* obj = objP.get();
  if(Entity* ent=dynamic_cast<Entity*>(obj)){
    ent->Apply(*this, false);
  }
  else if(Surface* sur=dynamic_cast<Surface*>(obj)){
    sur->Apply(*this, false);
  }
}

void UniformColorOp::SetColor(const gfx::Color& color)
{
  color_ = color;
}

gfx::Color UniformColorOp::GetColor() const
{
  return color_;
}

void UniformColorOp::ToInfo(info::InfoGroup& group) const
{
  info::InfoGroup super_group = group.CreateGroup("ColorOp");
  ColorOp::ToInfo(super_group);
  std::ostringstream ss;
  ss << color_.Red() << "\t" << color_.Green() << "\t" << color_.Blue() 
     << "\t" << color_.Alpha();
  group.SetTextData(ss.str());
}

gfx::UniformColorOp UniformColorOp::FromInfo(info::InfoGroup& group)
{
  info::InfoGroup super_group = group.GetGroup("ColorOp");
  gfx::ColorOp op = ColorOp::FromInfo(super_group);
  std::istringstream ss(group.GetTextData());

  float r, g, b, a;
  ss >> r >> g >> b >> a;
  Color c = Color(r,g,b,a);
  mol::QueryViewWrapper wrapper(op.GetSelection(),op.GetSelectionFlags());
  int mask = op.GetMask();
  return gfx::UniformColorOp(wrapper,mask,c);
}

}}

