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
  pixel sampling encapsulation

  Author: Ansgar Philippsen
*/

#ifndef IMG_PIXEL_SAMPLING_H
#define IMG_PIXEL_SAMPLING_H

#include <ost/img/module_config.hh>
#include <ost/img/data_types.hh>
#include "vecmat.hh"

#include <ost/message.hh>

namespace ost { namespace img {

class Extent;
class Point;

class DLLEXPORT InvalidSampling: public Error {
public:
  InvalidSampling():
    Error("invalid sampling value used, must be >0")
  {}
};

//! Helper class to handle pixel sampling
/*
  This class encapsulates the logic dealing with
  all aspects of domain-dependent pixel sampling.
  It is used by the derived Data classes to
  implement the sampling functionality.
*/
class DLLEXPORT_OST_IMG_BASE PixelSampling {
public:
  PixelSampling();
  // initialize with spatial sampling always!
  PixelSampling(const Vec3& spat_sampling, DataDomain d, const Extent& e);

  const Vec3& GetPixelSampling() const;
  void SetPixelSampling(const Vec3&);

  const Vec3& GetSpatialSampling() const;
  void SetSpatialSampling(const Vec3& d);

  const Vec3& GetFrequencySampling() const;
  void SetFrequencySampling(const Vec3& d);

  void SetDomain(DataDomain d);
  void SetExtent(const Extent& e);

  Vec3 Point2Coord(const Point& p) const;
  Vec3 Coord2Point(const Vec3& c) const;

  Vec3 Vec2Coord(const Vec3& p) const;

private:
  DataDomain dom_;
  Vec3 spat_scale_, i_spat_scale_;
  Vec3 freq_scale_, i_freq_scale_;
  Vec3 fac_;

  void set_spat_scale(const Vec3& d);
  void set_freq_scale(const Vec3& d);
};


}} // namespace

#endif
