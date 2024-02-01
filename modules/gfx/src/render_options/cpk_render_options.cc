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

#include "cpk_render_options.hh"

namespace ost { namespace gfx {

  CPKRenderOptions::CPKRenderOptions(): sphere_detail_(4),cpk_mode_(0),rad_mult_(1.0) {
#if OST_SHADER_SUPPORT_ENABLED
  cpk_mode_=1;
#else
  cpk_mode_=0;
#endif
}

RenderMode::Type CPKRenderOptions::GetRenderMode(){
 return RenderMode::CPK;
}

bool CPKRenderOptions::CanApplyRenderOptions(RenderOptionsPtr render_options){
  return render_options.get()->GetRenderMode()==RenderMode::CPK;
}

void CPKRenderOptions::ApplyRenderOptions(RenderOptionsPtr render_options){
  CPKRenderOptionsPtr options = boost::static_pointer_cast<CPKRenderOptions>(render_options);

  sphere_detail_=options->GetSphereDetail();
  cpk_mode_=options->GetSphereMode();
  this->NotifyStateChange();
}

void CPKRenderOptions::SetSphereDetail(uint sphere_detail){
  if(sphere_detail_ != sphere_detail){
    sphere_detail_ = sphere_detail;
    this->NotifyStateChange();
  }
}

uint CPKRenderOptions::GetSphereDetail(){
  return sphere_detail_;
}

void CPKRenderOptions::SetSphereMode(uint mode){
  if(cpk_mode_!=mode) {
    cpk_mode_=mode;
    this->NotifyStateChange();
  }
}

uint CPKRenderOptions::GetSphereMode(){
  return cpk_mode_;
}

void CPKRenderOptions::SetRadiusMult(float m)
{
  rad_mult_=std::max(0.0f,m);
  this->NotifyStateChange();
}

float CPKRenderOptions::GetRadiusMult() const
{
  return rad_mult_;
}

}} // ns
