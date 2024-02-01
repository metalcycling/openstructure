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
#include "line_render_options.hh"


namespace ost { namespace gfx {

LineRenderOptions::LineRenderOptions(): line_width_(2.0), aa_lines_(false){}

void LineRenderOptions::ApplyRenderOptions(RenderOptionsPtr render_options){
  LineRenderOptionsPtr options = boost::static_pointer_cast<LineRenderOptions>(render_options);

  line_width_=options->GetLineWidth();
  aa_lines_=options->GetAALines();
  this->NotifyStateChange();
}

void LineRenderOptions::SetLineWidth(float line_width){
  if(line_width_!=line_width) {
    line_width_ = line_width;
    this->NotifyStateChange();
  }
}

float LineRenderOptions::GetLineWidth(){
  return line_width_;
}

void LineRenderOptions::SetAALines(bool aa_lines){
  if(aa_lines_!=aa_lines) {
    aa_lines_=aa_lines;
    this->NotifyStateChange();
  }
}

bool LineRenderOptions::GetAALines(){
  return aa_lines_;
}

}} // ns
