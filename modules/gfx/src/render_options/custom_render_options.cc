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
/*
  Author: Stefan Scheuber
*/

#include "custom_render_options.hh"

namespace ost { namespace gfx {

CustomRenderOptions::CustomRenderOptions(): 
  sphere_detail_(4),
  arc_detail_(4),
  sphere_rad_(0.25),
  bond_rad_(0.25)
{}

RenderMode::Type CustomRenderOptions::GetRenderMode(){
 return RenderMode::CUSTOM;
}

bool CustomRenderOptions::CanApplyRenderOptions(RenderOptionsPtr render_options){
  return render_options.get()->GetRenderMode()==RenderMode::CUSTOM;
}

void CustomRenderOptions::ApplyRenderOptions(RenderOptionsPtr render_options){
  CustomRenderOptionsPtr options = boost::static_pointer_cast<CustomRenderOptions>(render_options);

  sphere_detail_=options->GetSphereDetail();
  arc_detail_=options->GetArcDetail();
  sphere_rad_=options->GetSphereRad();
  bond_rad_=options->GetBondRad();
  this->NotifyStateChange();
}

void CustomRenderOptions::SetSphereDetail(uint sphere_detail){
  if(sphere_detail_ != sphere_detail){
    sphere_detail_ = sphere_detail;
    this->NotifyStateChange();
  }
}

uint CustomRenderOptions::GetSphereDetail(){
  return sphere_detail_;
}

void CustomRenderOptions::SetArcDetail(uint arc_detail){
  if(arc_detail_ != arc_detail){
    arc_detail_ = arc_detail;
    this->NotifyStateChange();
  }
}

uint CustomRenderOptions::GetArcDetail(){
  return arc_detail_;
}

void CustomRenderOptions::SetSphereRad(float r)
{
  if(r!=sphere_rad_) {
    sphere_rad_=r;
    this->NotifyStateChange();
  }
}

float CustomRenderOptions::GetSphereRad() const
{
  return sphere_rad_;
}

void CustomRenderOptions::SetBondRad(float r)
{
  if(r!=bond_rad_) {
    bond_rad_=r;
    this->NotifyStateChange();
  }
}

float CustomRenderOptions::GetBondRad() const
{
  return bond_rad_;
}

}} // ns
