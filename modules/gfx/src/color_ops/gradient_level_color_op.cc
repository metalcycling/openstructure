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
#include <ost/gfx/color_ops/gradient_level_color_op.hh>

namespace ost { namespace gfx {

GradientLevelColorOp::GradientLevelColorOp() : GradientColorOp(), level_()
{
  this->Init();
}

GradientLevelColorOp::GradientLevelColorOp(const String& selection, 
                                            const String& property,
                                            const gfx::Gradient& gradient, 
                                            float minv, 
                                            float maxv, mol::Prop::Level level):
  GradientColorOp(selection, property, gradient, minv, maxv), level_(level)
{
  this->Init();
}

GradientLevelColorOp::GradientLevelColorOp(const String& selection, int mask,
                                            const String& property,
                                            const gfx::Gradient& gradient,
                                            float minv,
                                            float maxv, mol::Prop::Level level):
  GradientColorOp(selection, mask, property, gradient, minv, maxv), level_(level)
{
  this->Init();
}

GradientLevelColorOp::GradientLevelColorOp(const String& selection,
                                            const String& property,
                                            const gfx::Gradient& gradient,
                                            mol::Prop::Level level):
  GradientColorOp(selection, property, gradient), level_(level)
{
  this->Init();
}

GradientLevelColorOp::GradientLevelColorOp(const String& selection, int mask,
                                            const String& property,
                                            const gfx::Gradient& gradient,
                                            mol::Prop::Level level):
  GradientColorOp(selection, mask, property, gradient), level_(level)
{
  this->Init();
}

GradientLevelColorOp::GradientLevelColorOp(const mol::QueryViewWrapper& query_view,
                                            const String& property,
                                            const gfx::Gradient& gradient,
                                            float minv,
                                            float maxv, mol::Prop::Level level):
  GradientColorOp(query_view, property, gradient, minv, maxv), level_(level)
{
  this->Init();
}

GradientLevelColorOp::GradientLevelColorOp(const mol::QueryViewWrapper& query_view, int mask,
                                            const String& property,
                                            const gfx::Gradient& gradient,
                                            float minv,
                                            float maxv, mol::Prop::Level level):
  GradientColorOp(query_view, mask, property, gradient, minv, maxv), level_(level)
{
  this->Init();
}

GradientLevelColorOp::GradientLevelColorOp(const mol::QueryViewWrapper& query_view,
                                            const String& property,
                                            const gfx::Gradient& gradient,
                                            mol::Prop::Level level):
  GradientColorOp(query_view, property, gradient), level_(level)
{
  this->Init();
}

GradientLevelColorOp::GradientLevelColorOp(const mol::QueryViewWrapper& query_view, int mask,
                                            const String& property,
                                            const gfx::Gradient& gradient,
                                            mol::Prop::Level level):
  GradientColorOp(query_view, mask, property, gradient), level_(level)
{
  this->Init();
}

void GradientLevelColorOp::Init()
{
  this->SetName("Gradient level");
}

bool GradientLevelColorOp::CanApplyTo(const GfxObjP& obj) const{
  if(dynamic_cast<Entity*>(obj.get()))
    return true;
  else if(dynamic_cast<Surface*>(obj.get())){
    return true;
  }
  return false;
}

void GradientLevelColorOp::ApplyTo(GfxObjP& objP) const{
  GfxObj* obj = objP.get();
  if(Entity* ent=dynamic_cast<Entity*>(obj)){
    ent->Apply(*this, false);
  }
  else if(Surface* sur=dynamic_cast<Surface*>(obj)){
    sur->Apply(*this, false);
  }
}

void GradientLevelColorOp::SetLevel(mol::Prop::Level level)
{
  level_ = level;
}

mol::Prop::Level GradientLevelColorOp::GetLevel() const
{
  return level_;
}

void GradientLevelColorOp::ToInfo(info::InfoGroup& group) const
{
  info::InfoGroup super_group = group.CreateGroup("GradientColorOp");
  GradientColorOp::ToInfo(super_group);
  std::ostringstream ss;
  ss << int(level_);
  group.SetTextData(ss.str());
}

gfx::GradientLevelColorOp GradientLevelColorOp::FromInfo(info::InfoGroup& group)
{
  info::InfoGroup super_group = group.GetGroup("GradientColorOp");
  gfx::GradientColorOp gop = GradientColorOp::FromInfo(super_group);
  std::istringstream ss(group.GetTextData());
  int l;
  ss >> l;
  mol::Prop::Level level = mol::Prop::Level(l);
  mol::QueryViewWrapper wrapper(gop.GetSelection(),gop.GetSelectionFlags());
  String property = gop.GetProperty();
  int mask = gop.GetMask();
  gfx::Gradient gradient = gop.GetGradient();
  float minv = gop.GetMinV();
  float maxv = gop.GetMaxV();
  if(gop.GetCalculateMinMax()){
    return gfx::GradientLevelColorOp(wrapper, mask, property, gradient, level);
  }
  else{
    return gfx::GradientLevelColorOp(wrapper, mask, property, gradient, minv, maxv, level);
  }
}

}}

