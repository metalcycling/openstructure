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
#include "gradient.hh"
#include <iostream>
#include <vector>
#include <algorithm>

#include <ost/log.hh>

#include <ost/info/info.hh>
#include <ost/info/info_fw.hh>

#include <ost/gfx/gradient_manager.hh>

namespace ost { namespace gfx {


Gradient::Gradient():
  stops_(),
  hsv_mode_(false)
{}

Gradient::Gradient(const String& name)
{
  Gradient gradient = GradientManager::Instance().GetGradient(name);
  hsv_mode_ = gradient.hsv_mode_;
  // why doesn't this work:
  //  stops_=gradient.stops_
  // or even better
  //  *this = gradient
  StopList stops = gradient.GetStops();
  for(unsigned int i = 0; i < stops.size(); i++){
    Stop stop = stops[i];
    this->SetColorAt(stop.t, stop.color);
  }
}

void Gradient::SetColorAt(float t, const Color& c)
{
  std::vector<Stop>::iterator i;
  for (i=stops_.begin(); i!=stops_.end(); ++i) {
    if (i->t>t) {
      break;
    }
  }
  stops_.insert(i, Stop(t, c));
}

Color Gradient::GetColorAt(float t) const
{
  if (stops_.empty()) {
    return Color();
  } else if(stops_.size()==1) {
    return stops_[0].color;
  }
  
  uint c=0;
  while (c<stops_.size() && t>=stops_[c].t) {
    ++c;
  }
  if (c==0) {
    return stops_.front().color;
  } else if (c==stops_.size()) {
    return stops_.back().color;
  }
  float d = stops_[c].t-stops_[c-1].t;
  float tt = d>0.001  ? (t-stops_[c-1].t)/d : 0.0;
  if(hsv_mode_) {
    geom::Vec4 v=stops_[c-1].color.GetHSVA()*(1.0-tt)+stops_[c].color.GetHSVA()*tt;
    return HSVA(v[0],v[1],v[2],v[3]);
  } else {
    return stops_[c-1].color*(1.0-tt)+stops_[c].color*tt;
  }
}

void Gradient::GradientToInfo(info::InfoGroup& group) const
{
  std::ostringstream ss;

  ss << stops_.size() << "\n";

  for( StopList::size_type i = 0; i < stops_.size(); i++ ) {
    ss << stops_[i].t << "\t" << stops_[i].color.Red() << "\t" << stops_[i].color.Green() << "\t" << stops_[i].color.Blue() << "\t" << stops_[i].color.Alpha() << "\n";
  }
  group.SetTextData(ss.str());
  group.SetAttribute("hsv_mode",hsv_mode_? "1" : "0");
}

Gradient::StopList Gradient::GetStops() const{
	return stops_;
}

Gradient Gradient::GradientFromInfo(const info::InfoGroup& group)
{
  std::istringstream ss(group.GetTextData());
  bool hsv_color=false;
  if(group.HasAttribute("hsv_color")) {
    String hsv_color_s = group.GetAttribute("hsv_color");
    std::transform(hsv_color_s.begin(), hsv_color_s.end(), hsv_color_s.begin(), ::tolower);
    if(hsv_color_s=="0" || hsv_color_s=="false" || hsv_color_s=="no" || hsv_color_s=="off") {
      hsv_color=false;
    } else {
      hsv_color=true;
    }
  }
  Gradient grad;

  StopList::size_type size = 0;
  ss >> size;
  float t, r, g, b, a;
  for( StopList::size_type i = 0; i < size; i++ ) {
    ss >> t >> r >> g >> b >> a;
    Color c = hsv_color ? HSVA(r,g,b,a) : RGBA(r,g,b,a);
    grad.SetColorAt(t, c);
  }

  if(group.HasAttribute("hsv_mode")) {
    String hsv_mode_s = group.GetAttribute("hsv_mode");
    std::transform(hsv_mode_s.begin(), hsv_mode_s.end(), hsv_mode_s.begin(), ::tolower);
    if(hsv_mode_s=="0" || hsv_mode_s=="false" || hsv_mode_s=="no" || hsv_mode_s=="off") {
      grad.SetHSVMode(false);
    } else {
      grad.SetHSVMode(true);
    }
  }
  return grad;
}

}}
