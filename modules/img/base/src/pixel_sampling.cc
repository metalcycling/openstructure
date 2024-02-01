//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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
  pixel sampling

  Author: Ansgar Philippsen
*/

#include "pixel_sampling.hh"

#include "extent.hh"
#include "point.hh"
#include "vecmat.hh"

#include <ost/log.hh>

namespace ost { namespace img {

PixelSampling::PixelSampling():
  dom_(SPATIAL),
  spat_scale_(Vec3(1.0,1.0,1.0)),
  i_spat_scale_(Vec3(1.0,1.0,1.0)),
  freq_scale_(Vec3(1.0,1.0,1.0)),
  i_freq_scale_(Vec3(1.0,1.0,1.0)),
  fac_(Vec3(1.0,1.0,1.0))
{
  SetPixelSampling(Vec3(1.0,1.0,1.0));
}

PixelSampling::PixelSampling(const Vec3& sampling, DataDomain d, const Extent& e):
  dom_(d),
  spat_scale_(Vec3(1.0,1.0,1.0)),
  i_spat_scale_(Vec3(1.0,1.0,1.0)),
  freq_scale_(Vec3(1.0,1.0,1.0)),
  i_freq_scale_(Vec3(1.0,1.0,1.0)),
  fac_(Vec3(e.GetSize().GetWidth(),
	    e.GetSize().GetHeight(),
	    e.GetSize().GetDepth()))
{
  if(dom_==HALF_FREQUENCY) dom_=FREQUENCY;
  //SetPixelSampling(sampling);
  set_spat_scale(sampling);
}

const Vec3& PixelSampling::GetPixelSampling() const
{
  return (dom_==SPATIAL ? spat_scale_ : freq_scale_);
}

void PixelSampling::SetPixelSampling(const Vec3& d)
{
  dom_==SPATIAL ? set_spat_scale(d) : set_freq_scale(d);
}
  
const Vec3& PixelSampling::GetSpatialSampling() const
{
  return spat_scale_;
}

const Vec3& PixelSampling::GetFrequencySampling() const
{
  return freq_scale_;
}

void PixelSampling::SetSpatialSampling(const Vec3& s)
{
  set_spat_scale(s);
}

void PixelSampling::SetFrequencySampling(const Vec3& s)
{
  set_freq_scale(s);
}

void PixelSampling::SetDomain(DataDomain d)
{
  if(d==HALF_FREQUENCY) d=FREQUENCY;
  dom_=d;
}

void PixelSampling::SetExtent(const Extent& e)
{
  fac_=Vec3(e.GetSize().GetWidth(),
	    e.GetSize().GetHeight(),
	    e.GetSize().GetDepth());
}

Vec3 PixelSampling::Point2Coord(const Point& p) const
{
  return (dom_==SPATIAL ? CompMultiply(spat_scale_,p.ToVec3()) : CompMultiply(freq_scale_,p.ToVec3()));
}

Vec3 PixelSampling::Coord2Point(const Vec3& c) const
{
  return (dom_==SPATIAL ? CompMultiply(c,i_spat_scale_) : CompMultiply(c,i_freq_scale_));
}

Vec3 PixelSampling::Vec2Coord(const Vec3& p) const
{
  return (dom_==SPATIAL ? CompMultiply(spat_scale_,p) : CompMultiply(freq_scale_,p));
}


void PixelSampling::set_spat_scale(const Vec3& dim)
{
  if(dim[0]<=0.0 || dim[1]<0.0 || dim[2]<0.0) {
    LOG_ERROR("invalid spatial pixel sampling, must be >0");
    return;
  }
  spat_scale_=dim;
  i_spat_scale_=1.0/spat_scale_;
  freq_scale_=Vec3(1.0/(fac_[0]*dim[0]),1.0/(fac_[1]*dim[1]),1.0/(fac_[2]*dim[2]));
  i_freq_scale_ = 1.0/freq_scale_;
}

void PixelSampling::set_freq_scale(const Vec3& dim)
{
  freq_scale_ = dim;
  i_freq_scale_ = 1.0/dim;
  spat_scale_=Vec3(1.0/(fac_[0]*dim[0]),1.0/(fac_[1]*dim[1]),1.0/(fac_[2]*dim[2]));
  i_spat_scale_ = 1.0/spat_scale_;
}

}} // namespace
