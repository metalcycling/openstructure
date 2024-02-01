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

#include "cartoon_render_options.hh"

namespace ost { namespace gfx {

CartoonRenderOptions::CartoonRenderOptions(bool force_tube):
  force_tube_(force_tube),
  spline_detail_(6),
  poly_mode_(2),
  arc_detail_(4),
  smooth_factor_(0.0),
  tube_radius_(0.4),
  tube_ratio_(1.0),
  tube_profile_(0),
  helix_width_(1.1),
  helix_thickness_(0.2),
  helix_ecc_(0.3),
  helix_profile_(1),
  helix_mode_(0),
  strand_width_(1.2),
  strand_thickness_(0.2),
  strand_ecc_(0.3),
  strand_profile_(1),
  strand_mode_(0),
  color_blend_mode_(0),
  twist_hack_(true)
{}

RenderMode::Type CartoonRenderOptions::GetRenderMode(){
  if(force_tube_)
    return RenderMode::TUBE;
  return RenderMode::HSC;
}

bool CartoonRenderOptions::CanApplyRenderOptions(RenderOptionsPtr render_options){
  if(force_tube_)
    return render_options.get()->GetRenderMode()==RenderMode::TUBE;
  return render_options.get()->GetRenderMode()==RenderMode::HSC;
}

void CartoonRenderOptions::ApplyRenderOptions(RenderOptionsPtr render_options){
  CartoonRenderOptionsPtr options = boost::static_pointer_cast<CartoonRenderOptions>(render_options);

  spline_detail_=options->GetSplineDetail();
  poly_mode_=options->GetPolyMode();
  arc_detail_=options->GetArcDetail();
  smooth_factor_=options->GetNormalSmoothFactor();
  tube_radius_=options->GetTubeRadius();
  tube_ratio_=options->GetTubeRatio();
  tube_profile_=options->GetTubeProfileType();
  helix_width_=options->GetHelixWidth();
  helix_thickness_=options->GetHelixThickness();
  helix_ecc_=options->GetHelixEcc();
  helix_profile_=options->GetHelixProfileType();
  helix_mode_=options->GetHelixMode();
  strand_width_=options->GetStrandWidth();
  strand_thickness_=options->GetStrandThickness();
  strand_ecc_=options->GetStrandEcc();
  strand_profile_=options->GetStrandProfileType();
  strand_mode_=options->GetStrandMode();
  color_blend_mode_=options->GetColorBlendMode();
  twist_hack_=options->GetTwistHack();
  this->NotifyStateChange();
}

void CartoonRenderOptions::SetSplineDetail(uint spline_detail){
  if(spline_detail_ != spline_detail){
    spline_detail_ = spline_detail;
    this->NotifyStateChange();
  }
}

uint CartoonRenderOptions::GetSplineDetail() const{
  return spline_detail_;
}

void CartoonRenderOptions::SetPolyMode(uint poly_mode){
  if(poly_mode_ != poly_mode){
    poly_mode_ = poly_mode;
    this->NotifyStateChange();
  }
}

uint CartoonRenderOptions::GetPolyMode() const{
  return poly_mode_;
}

void CartoonRenderOptions::SetArcDetail(uint arc_detail){
  if(arc_detail_ != arc_detail){
    arc_detail_ = arc_detail;
    this->NotifyStateChange();
  }
}

uint CartoonRenderOptions::GetArcDetail() const{
  return arc_detail_;
}

void CartoonRenderOptions::SetNormalSmoothFactor(float smooth_factor){
  if(smooth_factor_ != smooth_factor){
    smooth_factor_ = smooth_factor;
    this->NotifyStateChange();
  }
}

float CartoonRenderOptions::GetNormalSmoothFactor() const{
  return smooth_factor_;
}

void CartoonRenderOptions::SetTubeRadius(float tube_radius){
  if(tube_radius_ != tube_radius){
    tube_radius_= tube_radius>0.0 ? tube_radius : tube_radius_;
    this->NotifyStateChange();
  }
}

float CartoonRenderOptions::GetTubeRadius() const{
  return tube_radius_;
}

void CartoonRenderOptions::SetTubeRatio(float tube_ratio){
  if(tube_ratio_ != tube_ratio){
    tube_ratio_= tube_ratio>0.0 ? tube_ratio : tube_ratio_;
    this->NotifyStateChange();
  }
}

float CartoonRenderOptions::GetTubeRatio() const{
  return tube_ratio_;
}

unsigned int CartoonRenderOptions::GetTubeProfileType() const
{
  return tube_profile_;
}

void CartoonRenderOptions::SetTubeProfileType(unsigned int t)
{
  tube_profile_=t;
  this->NotifyStateChange();
}

void CartoonRenderOptions::SetHelixWidth(float helix_width){
  if(helix_width_ != helix_width){
    helix_width_= helix_width>0.0 ? helix_width : helix_width_;
    this->NotifyStateChange();
  }
}

float CartoonRenderOptions::GetHelixWidth() const{
  return helix_width_;
}

void CartoonRenderOptions::SetHelixThickness(float helix_thickness){
  if(helix_thickness_ != helix_thickness){
    helix_thickness_= helix_thickness>0.0 ? helix_thickness : helix_thickness_;
    this->NotifyStateChange();
  }
}

float CartoonRenderOptions::GetHelixThickness() const{
  return helix_thickness_;
}

void CartoonRenderOptions::SetHelixEcc(float helix_ecc){
  if(helix_ecc_ != helix_ecc){
    helix_ecc_= helix_ecc>0.0 ? helix_ecc : helix_ecc_;
    this->NotifyStateChange();
  }
}

float CartoonRenderOptions::GetHelixEcc() const{
  return helix_ecc_;
}

unsigned int CartoonRenderOptions::GetHelixProfileType() const
{
  return helix_profile_;
}

void CartoonRenderOptions::SetHelixProfileType(unsigned int t)
{
  helix_profile_=t;
  this->NotifyStateChange();
}

unsigned int CartoonRenderOptions::GetHelixMode() const
{
  return helix_mode_;
}

void CartoonRenderOptions::SetHelixMode(unsigned int m)
{
  helix_mode_=m;
  this->NotifyStateChange();
}

void CartoonRenderOptions::SetStrandWidth(float strand_width){
  if(strand_width_ != strand_width){
    strand_width_= strand_width>0.0 ? strand_width : strand_width_;
    this->NotifyStateChange();
  }
}

float CartoonRenderOptions::GetStrandWidth() const{
  return strand_width_;
}

void CartoonRenderOptions::SetStrandThickness(float strand_thickness){
  if(strand_thickness_ != strand_thickness){
    strand_thickness_= strand_thickness>0.0 ? strand_thickness : strand_thickness_;
    this->NotifyStateChange();
  }
}

float CartoonRenderOptions::GetStrandThickness() const{
  return strand_thickness_;
}

void CartoonRenderOptions::SetStrandEcc(float strand_ecc){
  if(strand_ecc_ != strand_ecc){
    strand_ecc_= strand_ecc>0.0 ? strand_ecc : strand_ecc_;
    this->NotifyStateChange();
  }
}

float CartoonRenderOptions::GetStrandEcc() const{
  return strand_ecc_;
}

unsigned int CartoonRenderOptions::GetStrandProfileType() const
{
  return strand_profile_;
}

void CartoonRenderOptions::SetStrandProfileType(unsigned int t)
{
  strand_profile_=t;
  this->NotifyStateChange();
}

unsigned int CartoonRenderOptions::GetStrandMode() const
{
  return strand_mode_;
}

void CartoonRenderOptions::SetStrandMode(unsigned int m)
{
  strand_mode_=m;
  this->NotifyStateChange();
}

void CartoonRenderOptions::SetColorBlendMode(unsigned int m)
{
  color_blend_mode_=m;
  this->NotifyStateChange();
}
unsigned int CartoonRenderOptions::GetColorBlendMode() const
{
  return color_blend_mode_;
}

float CartoonRenderOptions::GetMaxRad() const{
  float max_rad=std::max(float(3.0),tube_radius_*tube_ratio_);
  max_rad=std::max(max_rad,tube_radius_);
  max_rad=std::max(max_rad,helix_width_);
  max_rad=std::max(max_rad,helix_thickness_);
  max_rad=std::max(max_rad,strand_width_);
  max_rad=std::max(max_rad,strand_thickness_);
  return max_rad;
}

void CartoonRenderOptions::SetTwistHack(bool f) 
{
  twist_hack_=f;
  this->NotifyStateChange();
}

}} // ns

