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
#include <ost/dyn_cast.hh>

#include <ost/gfx/color_ops/gradient_color_op.hh>

namespace ost { namespace gfx {

GradientColorOp::GradientColorOp() : ColorOp(), property_(), gradient_(){
  this->Init();
}

GradientColorOp::GradientColorOp(const String& selection, const String& property,
                                 const gfx::Gradient& gradient, 
                                 float minv, float maxv):
  ColorOp(selection), property_(property), gradient_(gradient), clamp_(true), calculate_(false),
  minv_(minv), maxv_(maxv)
{
  this->Init();
}

GradientColorOp::GradientColorOp(const String& selection, int mask, const String& property,
                                 const gfx::Gradient& gradient,
                                 float minv, float maxv):
  ColorOp(selection,mask), property_(property), gradient_(gradient), clamp_(true), calculate_(false),
  minv_(minv), maxv_(maxv)
{
  this->Init();
}

GradientColorOp::GradientColorOp(const String& selection, const String& property,
                                 const gfx::Gradient& gradient):
  ColorOp(selection), property_(property), gradient_(gradient), clamp_(true), calculate_(true)
{
  this->Init();
}

GradientColorOp::GradientColorOp(const String& selection, int mask, const String& property,
                                 const gfx::Gradient& gradient):
  ColorOp(selection,mask), property_(property), gradient_(gradient), clamp_(true), calculate_(true)
{
  this->Init();
}

GradientColorOp::GradientColorOp(const mol::QueryViewWrapper& query_view, const String& property,
                                 const gfx::Gradient& gradient,
                                 float minv, float maxv):
  ColorOp(query_view), property_(property), gradient_(gradient), clamp_(true), calculate_(false),
  minv_(minv), maxv_(maxv)
{
  this->Init();
}

GradientColorOp::GradientColorOp(const mol::QueryViewWrapper& query_view, int mask, const String& property,
                                 const gfx::Gradient& gradient,
                                 float minv, float maxv):
  ColorOp(query_view,mask), property_(property), gradient_(gradient), clamp_(true), calculate_(false),
  minv_(minv), maxv_(maxv)
{
  this->Init();
}

GradientColorOp::GradientColorOp(const mol::QueryViewWrapper& query_view, const String& property,
                                 const gfx::Gradient& gradient):
  ColorOp(query_view), property_(property), gradient_(gradient), clamp_(true), calculate_(true)
{
  this->Init();
}

GradientColorOp::GradientColorOp(const mol::QueryViewWrapper& query_view, int mask, const String& property,
                                 const gfx::Gradient& gradient):
  ColorOp(query_view,mask), property_(property), gradient_(gradient), clamp_(true), calculate_(true)
{
  this->Init();
}

void GradientColorOp::Init()
{
  this->SetName("Gradient");
}

void GradientColorOp::SetProperty(const String& property)
{
  property_ = property;
}

String GradientColorOp::GetProperty() const
{
  return property_;
}

void GradientColorOp::SetGradient(const gfx::Gradient& gradient)
{
  gradient_ = gradient;
}

gfx::Gradient GradientColorOp::GetGradient() const
{
  return gradient_;
}

void GradientColorOp::SetMaxV(float maxv) const
{
  maxv_ = maxv;
}

float GradientColorOp::GetMaxV() const
{
  return maxv_;
}

bool GradientColorOp::GetCalculateMinMax() const
{
  return calculate_;
}

void GradientColorOp::SetMinV(float minv) const
{
  minv_ = minv;
}
float GradientColorOp::GetMinV() const
{
  return minv_;
}

void GradientColorOp::ToInfo(info::InfoGroup& group) const
{
  info::InfoGroup super_group = group.CreateGroup("ColorOp");
  ColorOp::ToInfo(super_group);
  info::InfoGroup gradient_group = group.CreateGroup("Gradient");
  gradient_.GradientToInfo(gradient_group);
  std::ostringstream ss;
  ss << calculate_ << "\t" << minv_ << "\t" << maxv_ << "\t" << property_;
  group.SetTextData(ss.str());
}

gfx::GradientColorOp GradientColorOp::FromInfo(info::InfoGroup& group)
{
  info::InfoGroup super_group = group.GetGroup("ColorOp");
  gfx::ColorOp op = ColorOp::FromInfo(super_group);
  int mask = op.GetMask();
  info::InfoGroup gradient_group = group.GetGroup("Gradient");
  gfx::Gradient gradient = gfx::Gradient::GradientFromInfo(gradient_group);
  std::istringstream ss(group.GetTextData());
  bool calculate;
  float minv, maxv;
  ss >> calculate >> minv >> maxv;
  String property;
  String part;
  while(ss.good()){
    ss >> part;
    property.append(part);
    if(ss.good())property.append(" ");
    part.clear();
  }
  mol::QueryViewWrapper wrapper(op.GetSelection(),op.GetSelectionFlags());
  if(calculate) {
    return gfx::GradientColorOp(wrapper,mask,property,gradient);
  }
  else {
    return gfx::GradientColorOp(wrapper,mask,property,gradient, minv, maxv);
  }
}

}}

