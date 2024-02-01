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
#include "sline_render_options.hh"

namespace ost { namespace gfx {

SlineRenderOptions::SlineRenderOptions():
  spline_detail_(6)
{}

RenderMode::Type SlineRenderOptions::GetRenderMode(){
 return RenderMode::SLINE;
}

bool SlineRenderOptions::CanApplyRenderOptions(RenderOptionsPtr render_options){
  return render_options.get()->GetRenderMode()==RenderMode::SLINE;
}

void SlineRenderOptions::ApplyRenderOptions(RenderOptionsPtr render_options){
  SlineRenderOptionsPtr options = boost::static_pointer_cast<SlineRenderOptions>(render_options);

  LineRenderOptions::ApplyRenderOptions(render_options);
  spline_detail_=options->GetSplineDetail();
}

void SlineRenderOptions::SetSplineDetail(uint spline_detail){
  if(spline_detail_!=spline_detail){
    spline_detail_=spline_detail;
    this->NotifyStateChange();
  }
}

uint SlineRenderOptions::GetSplineDetail(){
  return spline_detail_;
}

}} // ns
