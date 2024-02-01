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
  image state base

  Author: Ansgar Philippsen
*/

#ifndef IMAGE_STATE_BASE_H
#define IMAGE_STATE_BASE_H

#include "image_state_base_fw.hh"

#include <ost/base.hh>
#include <ost/img/module_config.hh>
#include <ost/img/data.hh>

#include "image_state_visitor_fw.hh"
#include "type_fw.hh"

namespace ost { namespace img { 

class PixelSampling;

namespace image_state {

/*
  Abstract base class for the image state. A polymorphic pointer
  to this is managed by instantiations of the Image class.

  At this level, the SpatialOrigin and SpatialSampling are
  stored, but the interpretation is left to the 'user' of
  the image state, either the Image class or any ImageStateVisitor
*/
class DLLEXPORT_OST_IMG_BASE ImageStateBase {
public:

  virtual ~ImageStateBase() {}

  // provide deep copy
  virtual ImageStateBasePtr Clone(bool cc=true) const = 0;

  virtual size_t MemSize() const = 0;

  ///////////////////////////////////////
  // abstract virtual methods
  // to be implemented in ImageState

  // Reset extent, only re-allocating memory
  /* 
     still needed ? its probably better to offer this
     on the lower level of ImageStateImpl
   */
  //virtual void ResetExtent(const Extent& e) = 0;

  //! Retrieve data type: REAL or COMPLEX
  virtual DataType GetType() const = 0;

  //! Retrieve domain
  virtual DataDomain GetDomain() const = 0;

  //! Set spatial origin
  virtual void SetSpatialOrigin(const Point& o) = 0;

  //! Retrieve spatial origin
  virtual Point GetSpatialOrigin() const = 0;


  //! Retrieve Extent
  /*!
    Given by size, spatial origin and data domain, therefore the extent
    is a read-only property. The extent that is returned here corresponds
    to the physical extent, ie for half-frequency data only one half
  */
  virtual Extent GetExtent() const = 0;

  //! Retrieve logical extent
  /*!
    For half-frequency data, returns the full, virtual extent. For 
    spatial and frequency domain equivalent to GetExtent().
  */
  virtual Extent GetLogicalExtent() const = 0;

  //! retrieve value at Point as Real
  /*!
    This call is potentially slow, but boundary checked
  */
  virtual Real GetReal(const Point &p) const = 0;

  //! set real value at Point
  /*!
    This call is potentially slow, but boundary checked
  */
  virtual void SetReal(const Point &p, const Real& r) = 0;

  // retrive value at Point as Complex
  /*
    This call is potentially slow, but boundary checked
  */
  virtual Complex GetComplex(const Point &p) const = 0;

  //! set complex value at Point
  /*!
    This call is potentially slow, but boundary checked
  */
  virtual void SetComplex(const Point &p, const Complex& c) = 0;

  virtual Real GetIntpolReal(const Vec3 &v) const = 0;
  virtual Real GetIntpolReal(const Vec2 &v) const = 0;
  virtual Real GetIntpolReal(const Real &d) const = 0;

  virtual Complex GetIntpolComplex(const Vec3 &v) const = 0;
  virtual Complex GetIntpolComplex(const Vec2 &v) const = 0;
  virtual Complex GetIntpolComplex(const Real &d) const = 0;

  virtual PixelSampling& GetSampling() = 0;

  virtual const PixelSampling& GetSampling() const = 0;

  virtual void SetSampling(const PixelSampling& s) = 0;

  virtual Vec3 GetAbsoluteOrigin() const = 0;

  virtual void SetAbsoluteOrigin(const Vec3& c) = 0;

    // Index to coordinate conversion
  /*
    Returns the physical dimensions of a given Point based on
    the sampling of the current domain and the absolute coordinates
    of the origin of the image
  */
  virtual Vec3 IndexToCoord(const Point &p) const = 0;

  // Convert coordinates back to (fractional) indices
  /*
    Returns the location within the Data that corresponds to
    the given physical coordinate, bases on the sampling of
    the current domain and the absolute coordinates
    of the origin of the image
  */
  virtual Vec3 CoordToIndex(const Vec3& c) const = 0;

  /*
    Returns the physical dimensions of a given Vec which contains
    fractional Point coordinated. The returned value is based on
    the sampling of the current domain and the absolute coordinates
    of the origin of the image
  */
  virtual Vec3 FractionalIndexToCoord(const Vec3 &p) const =0;

  // visitor interface

  // apply const visitor, ie to extract info
  virtual void ApplyIP(ImageStateNonModVisitorBase& v) const = 0;
  // apply const visitor (same as ApplyIP(ImageStateNonModVisitorBase&)
  virtual void Apply(ImageStateNonModVisitorBase& v) const = 0;

  // apply in-place visitor, usually modifies this
  virtual void ApplyIP(ImageStateModIPVisitorBase& v) = 0;
  // apply in-place visitor to clone of this, return result
  virtual ImageStateBasePtr Apply(ImageStateModIPVisitorBase& v) const = 0;

  // apply in-place const visitor, usually modifies this
  virtual void ApplyIP(const ImageStateConstModIPVisitorBase& v) = 0;
  // apply in-place const visitor to clone of this, return result
  virtual ImageStateBasePtr Apply(const ImageStateConstModIPVisitorBase& v) const = 0;

  // apply out-of-place visitor and return result
  // NOTE: it is not possible to apply an out-of-place visitor in-place!
  // that is handled by the calling API (ie ImageHandle)
  virtual ImageStateBasePtr Apply(ImageStateModOPVisitorBase& v) const = 0;

  // apply out-of-place const visitor and return result
  virtual ImageStateBasePtr Apply(const ImageStateConstModOPVisitorBase& v) const = 0;

  // apply morph visitor, with the expectation that returning state has modified this
  virtual ImageStateBasePtr Apply(ImageStateMorphVisitorBase& v) = 0;

  // operators
  virtual ImageStateBase& operator+=(Real v) = 0;
  virtual ImageStateBase& operator+=(const Complex& v) = 0;
  virtual ImageStateBase& operator-=(Real v) = 0;
  virtual ImageStateBase& operator-=(const Complex& v) = 0;
  virtual ImageStateBase& operator*=(Real v) = 0;
  virtual ImageStateBase& operator*=(const Complex& v) = 0;
  virtual ImageStateBase& operator/=(Real v) = 0;
  virtual ImageStateBase& operator/=(const Complex& v) = 0;

  ImageStateBase& operator+=(const ImageStateBase& b);
  ImageStateBase& operator-=(const ImageStateBase& b);
  ImageStateBase& operator*=(const ImageStateBase& b);
  ImageStateBase& operator/=(const ImageStateBase& b);

protected:
  ImageStateBase() {}
  ImageStateBase(const ImageStateBase& s) {}

private:
  // assignement is forbidden
  ImageStateBase& operator=(const ImageStateBase& s);
};

}}} // namespaces

#endif
