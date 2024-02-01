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
  abstract base class for generic data

  Author: Ansgar Philippsen
*/

#ifndef IMG_DATA_H
#define IMG_DATA_H

#include <cmath>

#include <boost/shared_ptr.hpp>

#include <ost/img/module_config.hh>
#include <ost/base.hh>
#include "point.hh"
#include "extent.hh"
#include "data_types.hh"
#include "extent_iterator.hh"
#include "vecmat.hh"


namespace ost { namespace img {


//using namespace info;

// forward declaration
class NonModAlgorithm;
class ModIPAlgorithm;
class DataObserver;
class PixelSampling;
class ImageHandle;

//! Abstract base class for read-only data
/*!
  This class represents the abstract concept of a 
  read-only (const) data source which allows:

   - retrieving values based on triple integer indices
   - retrieving type, domain, and extent
   - applying NonModAlgorithms
   - (de-)attaching DataObservers
   - retrieving the spatial origin
   - retrieving pixel sampling

  In addition, a meta-info (InfoHandle) is associated
  with this data
*/
class DLLEXPORT_OST_IMG_BASE ConstData {
public:
  /*! @name Initialization and Deconstruction
   */
  //@{
  ConstData();
  virtual ~ConstData();
  //@}

  /*! @name Properties
  */
  //@{
  //! Return type of data: \c REAL or \c COMPLEX
  virtual DataType GetType() const = 0;

  //! Return domain: \c SPATIAL, \c FREQUENCY or \c HALF_FREQUENCY
  virtual DataDomain GetDomain() const = 0;

  //! Returns true if data is currently in spatial domain 
  /*!
    Convenience wrapper around GetDomain()
  */
  bool IsSpatial() const {return GetDomain()==SPATIAL;}

  //! Returns true if data is currently in frequency domain
  /*!
    Convenience wrapper around GetDomain()
  */
  bool IsFrequency() const {
    return GetDomain()==FREQUENCY || GetDomain()==HALF_FREQUENCY;
  }

  //! Returns true if data is currently of real type 
  /*!
    Convenience wrapper around GetType()
  */
  bool IsReal() const {return GetType()==REAL;}

  //! Returns true if data is currently of complex type
  /*!
    Convenience wrapper around GetType()
  */
  bool IsComplex() const {return GetType()==COMPLEX;}

  //! Return (logical) extent of data
  virtual Extent GetExtent() const = 0;

  //! Retrieve spatial origin
  virtual Point GetSpatialOrigin() const = 0;

  //! Return Size
  /*!
    Convenience method
  */
  Size GetSize() const {return GetExtent().GetSize();}

  //@}

  /*! @name Value retrieval
    The pure virtual methods GetReal() and GetComplex() form the basis for the
    interpolation and normalization methods.
  */
  //@{
  //! Get explicit real value
  virtual Real GetReal(const Point &p) const = 0;

  //! Get explicit complex value
  virtual Complex GetComplex(const Point &p) const = 0;

  //! Get interpolated real value from 3D,2D or 1D vec
  virtual Real GetIntpolReal(const Vec3 &v) const = 0;
  virtual Real GetIntpolReal(const Vec2 &v) const = 0;
  virtual Real GetIntpolReal(const Real &d) const = 0;

  //! Get interpolated complex value from 3D,2D or 1D vec
  virtual Complex GetIntpolComplex(const Vec3 &v) const = 0;
  virtual Complex GetIntpolComplex(const Vec2 &v) const = 0;
  virtual Complex GetIntpolComplex(const Real &d) const = 0;


  /*! @name Algorithm interface
    On the level of Data, only NonModAlgorithms may be applied
  */
  //@{
  //! Algorithm interface: NonModAlgorithm
  /*
    Only a NonModAlgorithm can be applied to a Data instance,
    the Mod IP/OP variants are part of the ImageHandle interface
  */
  virtual void Apply(NonModAlgorithm& a) const = 0;
  virtual void ApplyIP(NonModAlgorithm& a) const=0;
  //@}

  
  /*! @name Observer interface
   */
  //@{
  //! DataObserver interface: attach observer
  virtual void Attach(DataObserver *o) const = 0;

  //! DataObserver interface: detach observer
  virtual void Detach(DataObserver *o) const = 0;

  //! Notify each observer that state has changed
  virtual void Notify() const = 0;
  //@}

  /*! @name Misc
    should really be located somewhere else...
   */
  //@{
  virtual Real OverallDifference(const ConstData& d) const;
  virtual Real NormDifference(const ConstData& d) const;
  //@}


  /*! @name Sampling interface
    The pixel sampling depends on the current domain (\c SPATIAL or 
    \c FREQUENCY ) and may change, e.g. after a fourier-transform.
    Domain-independent access to the sampling is provided
    by the Set/GetSpatialSampling methods.
    The IndexToCoord() and CoordToIndex() methods can be used for conversion
    between pixel values and coordinates; they honour the current domain.
  */
  //@{

  //! return pixel sampling of the current domain
  Vec3 GetPixelSampling() const;

  //! return spatial pixel dimensions
  Vec3 GetSpatialSampling() const;

  //! return frequency pixel dimensions
  Vec3 GetFrequencySampling() const;

  // Index to coordinate conversion
  /*
    Returns the physical dimensions of a given Point based on
    the sampling of the current domain.
  */
  //Vec3 IndexToCoord(const Point &p) const;

  // Convert coordinates back to (fractional) indices
  /*
    Returns the location within the Data that corresponds to
    the given physical coordinate.
  */
  //Vec3 CoordToIndex(const Vec3& c) const;
  //@}



 protected:
  ConstData(const ConstData& d);
  ConstData& operator=(const ConstData& d);

  /*! @name Sampling implementation
    A reference to a PixelSampling instance must be provided by subclasses
    in order for the sampling interface to work.
  */
  //@{
  virtual PixelSampling& Sampling() = 0;
  virtual const PixelSampling& Sampling() const = 0;
  //@}

  //! implements interpolation
  void CalcIntpolWeights(Vec3 v, Real w[8], Point p[8]) const;

private:

};


//! Abstract base class for data
/*!
  The non-const Data class is expanded slightly in functionality, 
  as it offers some modification routines, namely

   - setting the spatial origin
   - modifying the pixel sampling

  It servers mainly as the base class for ImageHandle.
*/
class DLLEXPORT_OST_IMG_BASE Data: public ConstData {
public:
  /*! @name Initialization and Deconstruction
   */
  //@{
  Data();
  virtual ~Data();
  //@}

  
  /*! @name Properties
  */
  //@{
  //! Set the spatial origin
  virtual void SetSpatialOrigin(const Point& o) = 0;
  //@}

  /*! @name Sampling interface
  */
  //@{

  void SetPixelSampling(const Vec3&);
  //! Set uniform pixel dimensions of the current domain
  void SetPixelSampling(Real d) {SetPixelSampling(Vec3(d,d,d));}
  //! set spatial pixel dimensions (non-uniform)
  void SetSpatialSampling(const Vec3&);
  //! set spatial pixel dimensions (uniform)
  void SetSpatialSampling(Real d) {SetSpatialSampling(Vec3(d,d,d));}

  virtual void Apply(NonModAlgorithm& a) const=0;
  virtual void ApplyIP(NonModAlgorithm& a) const=0;

 protected:
  Data(const Data& d);
  Data& operator=(const Data& d);
};



}} // namespace img

#endif
