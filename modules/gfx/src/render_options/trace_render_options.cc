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
  Author: Marco Biasini
*/

#include "trace_render_options.hh"

namespace ost { namespace gfx {

TraceRenderOptions::TraceRenderOptions(): 
  spline_detail_(6),
  arc_detail_(4),
  smooth_factor_(0.0),
  tube_radius_(0.4)
{}

RenderMode::Type TraceRenderOptions::GetRenderMode()
{
 return RenderMode::TRACE;
}

bool TraceRenderOptions::CanApplyRenderOptions(RenderOptionsPtr render_options)
{
  return render_options.get()->GetRenderMode()==RenderMode::TRACE;
}

void TraceRenderOptions::ApplyRenderOptions(RenderOptionsPtr render_options)
{
  TraceRenderOptionsPtr options = boost::static_pointer_cast<TraceRenderOptions>(render_options);
  
  arc_detail_=options->GetArcDetail();
  smooth_factor_=options->GetNormalSmoothFactor();
  tube_radius_=options->GetTubeRadius();
  this->NotifyStateChange();
}

void TraceRenderOptions::SetArcDetail(uint arc_detail){
  if(arc_detail_ != arc_detail){
    arc_detail_ = arc_detail;
    this->NotifyStateChange();
  }
}

uint TraceRenderOptions::GetArcDetail() const{
  return arc_detail_;
}

void TraceRenderOptions::SetNormalSmoothFactor(float smooth_factor){
  if(smooth_factor_ != smooth_factor){
    smooth_factor_ = smooth_factor;
    this->NotifyStateChange();
  }
}

float TraceRenderOptions::GetNormalSmoothFactor() const{
  return smooth_factor_;
}

void TraceRenderOptions::SetTubeRadius(float tube_radius){
  if(tube_radius_ != tube_radius){
    tube_radius_= tube_radius>0.0 ? tube_radius : tube_radius_;
    this->NotifyStateChange();
  }
}

float TraceRenderOptions::GetTubeRadius() const{
  return tube_radius_;
}

}}
