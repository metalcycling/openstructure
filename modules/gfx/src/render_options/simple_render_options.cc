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

#include "simple_render_options.hh"

namespace ost { namespace gfx {

SimpleRenderOptions::SimpleRenderOptions(): 
  blur_flag_(false), blur_factors_(1.0, 0.8),
  bond_order_flag_(false), bond_order_distance_(0.06)
{ }

RenderMode::Type SimpleRenderOptions::GetRenderMode()
{
 return RenderMode::SIMPLE;
}

bool SimpleRenderOptions::CanApplyRenderOptions(RenderOptionsPtr render_options){
  return render_options.get()->GetRenderMode()==RenderMode::SIMPLE;
}

bool SimpleRenderOptions::GetBlurFlag() const
{
  return blur_flag_;
}

bool SimpleRenderOptions::GetBondOrderFlag() const
{
  return bond_order_flag_;
}

Real SimpleRenderOptions::GetBondOrderDistance() const
{
  return bond_order_distance_;
}

void SimpleRenderOptions::SetBlurFlag(bool flag)
{
  blur_flag_=flag;
  this->NotifyStateChange();
}

void SimpleRenderOptions::SetBondOrderFlag(bool flag)
{
  bond_order_flag_=flag;
  this->NotifyStateChange();
}

void SimpleRenderOptions::SetBondOrderDistance(Real bod)
{
  bond_order_distance_=bod;
  this->NotifyStateChange();
}

const std::pair<Real, Real>& SimpleRenderOptions::GetBlurFactors() const
{
  return blur_factors_;
}

void SimpleRenderOptions::SetBlurFactors(Real bf1, Real bf2)
{
  blur_factors_=std::make_pair(bf1, bf2);
  this->NotifyStateChange();
}

}} // ns

