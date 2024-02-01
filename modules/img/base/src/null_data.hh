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
  convenience null data class

  Author: Ansgar Philippsen
*/

#ifndef IMG_NULL_DATA_H
#define IMG_NULL_DATA_H

#include "data.hh"
#include "pixel_sampling.hh"
#include "image_factory.hh"
#include "image_handle.hh"

namespace ost { namespace img {

//! null data singleton
/*!
  Provides a full implementation of the Data interface, 
  but without providing any real functionality; all
  writes are ignored, and all reads return zero or
  otherwise null defaults.
*/
class DLLEXPORT_OST_IMG_BASE NullData: public Data
{
public:
  virtual DataType GetType() const {return REAL;}

  virtual DataDomain GetDomain() const {return SPATIAL;}

  virtual Extent GetExtent() const {return Extent();}

  virtual Real GetReal(const Point &p) const {return 0.0;}

  virtual Complex GetComplex(const Point &p) const {return 0.0;}

  virtual Real GetIntpolReal(const Vec3 &v) const {return 0.0;}
  virtual Real GetIntpolReal(const Vec2 &v) const {return 0.0;}
  virtual Real GetIntpolReal(const Real &d) const {return 0.0;}

  virtual Complex GetIntpolComplex(const Vec3 &v) const {return 0.0;}
  virtual Complex GetIntpolComplex(const Vec2 &v) const {return 0.0;}
  virtual Complex GetIntpolComplex(const Real &d) const {return 0.0;}

  virtual void Apply(NonModAlgorithm& a) const {}
  virtual void ApplyIP(NonModAlgorithm& a) const {}

  virtual void SetSpatialOrigin(const Point& o) {}

  virtual Point GetSpatialOrigin() const {return Point();}

  virtual void Attach(DataObserver *o) const {}

  virtual void Detach(DataObserver *o) const {}

  virtual void Notify() const {}

  virtual ImageHandle Crop(const Extent& e) const {
    return CreateImage();
  }

  static Data& Instance() {
    static NullData instance;
    return instance;
  }

protected:
  virtual PixelSampling& Sampling() {
    return sampling_;
  }
  virtual const PixelSampling& Sampling() const {
    return sampling_;
  }

private:
  PixelSampling sampling_;
};

}} // namespace

#endif
