//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2010 by the OpenStructure authors
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
  Authors: Ansgar Philippsen, Andreas Schenk
*/

#ifndef IMAGE_STATE_IMPL_HH
#define IMAGE_STATE_IMPL_HH

#include <iostream>
#include <boost/shared_ptr.hpp>

#include <ost/base.hh>
#include <ost/iplt/value_util.hh>
#include <ost/iplt/pixel_sampling.hh>

#include "image_state_base.hh"
#include "image_state_visitor_fw.hh"
#include "value_holder.hh"

namespace ost { namespace iplt {

class Size;

namespace image_state {

/*
  Provides concrete implementation of abstract 
  ImageStateBase interface. The two template
  parameters allow a combination of value
  and domain to be specified, which is
  utilized in the image state algorithms
  as well as the binary dispatch code.

  Internally, the ValueHolder and Domain interplay
  to provide the functionality

*/
// \internal
template <typename T, class D>
class TEMPLATE_EXPORT ImageStateImpl: public ImageStateBase
{
public:
  typedef ImageStateImpl<T,D> StateType;
  typedef boost::shared_ptr<ImageStateImpl<T,D> > SharedPtrType;

  ImageStateImpl();

  ImageStateImpl(const ImageStateImpl<T,D>& s);

  ImageStateImpl<T,D>& operator=(const ImageStateImpl<T,D>& s);

  // requires logical extent, not physical one!
  ImageStateImpl(const Extent& e, const PixelSampling& s);

  // also takes absolute center, requires logical extent, not physical one!
  ImageStateImpl(const Extent& e, const PixelSampling& s, const Vec3& c);

  /* 
     special ctor to utilize pre-allocated memory: the data area allocated
     by ValueHolder will be transfered to this image state
  */
  template <typename U>
  ImageStateImpl(const Extent& e, ValueHolder<U>& d, const PixelSampling& s):
    domain_(e),
    data_(domain_.GetExtent().GetSize(),d),
    sampling_(s)
  {
    assert(d.GetPhysicalSize()*sizeof(U)==data_.GetPhysicalSize()*sizeof(T));
  }

  virtual ~ImageStateImpl();

  SharedPtrType CloneState(bool cc=true) const;

  // image state base interface
  virtual ImageStateBasePtr Clone(bool cc=true) const;
  
  virtual long MemSize() const;

  virtual DataType GetType() const;

  virtual DataDomain GetDomain() const;

  virtual void SetSpatialOrigin(const Point& o);

  virtual void CenterSpatialOrigin();

  virtual Point GetSpatialOrigin() const;

  virtual Extent GetExtent() const;

  virtual Extent GetLogicalExtent() const;

  virtual Real GetReal(const Point &p) const;

  virtual void SetReal(const Point &p, const Real& r);

  virtual Complex GetComplex(const Point &p) const;

  virtual void SetComplex(const Point &p, const Complex& c);

  virtual Real GetIntpolReal(const Vec3 &v) const;
  virtual Real GetIntpolReal(const Vec2 &v) const;
  virtual Real GetIntpolReal(const Real &d) const;

  virtual Complex GetIntpolComplex(const Vec3 &v) const;
  virtual Complex GetIntpolComplex(const Vec2 &v) const;
  virtual Complex GetIntpolComplex(const Real &d) const;

  virtual PixelSampling& GetSampling();

  virtual const PixelSampling& GetSampling() const;

  virtual void SetSampling(const PixelSampling& s);

  virtual Vec3 GetAbsoluteOrigin() const;

  virtual void SetAbsoluteOrigin(const Vec3& c);

  virtual Vec3 IndexToCoord(const Point &p) const;

  virtual Vec3 CoordToIndex(const Vec3& c) const;

  virtual Vec3 FractionalIndexToCoord(const Vec3 &p) const;

  void AdjustPhaseOrigin(const Point& p);

  // value interface
  //////////////////////////////////
  // non-virtual member functions

  // retrieve value at specified point
  /*
    First converts point to index using domain policy,
    then uses ValueHolderType::Value member function to
    retrieve actual value. No boundary check is performed
    here, this is the responsibility of the caller!
  */
  T& Value(const Point& p);

  // retrieve ro value at specified point
  const T& Value(const Point& p) const;

  //! retrieve boundary checked value
  T GetCheckedValue(const Point& p) const;

  //! set boundary checked value
  void SetCheckedValue(const Point& p, const T& v);

  // calculate interpolated value at given point
  /*
    out-of-bounds values are automatically set to zero
  */
  T CalcIntpolValue(const Vec3& v) const;
  T CalcIntpolValue(const Vec2& v) const;
  T CalcIntpolValue(const Real& v) const;

  // pass-through to value holder
  T& Value(const Index& i);

  // pass-through to value holder
  const T& Value(const Index& i) const;

  // pass-through to value holder
  T& Value(unsigned int i);

  // pass-through to value holder
  const T& Value(unsigned int i) const;

  //! direct access to value holder
  ValueHolder<T>& Data();
  //! direct access to value holder
  const ValueHolder<T>& Data() const;

  // other
  Size GetSize() const;

  // visitor interface

  virtual void ApplyIP(ImageStateNonModVisitorBase& v) const;

  virtual void Apply(ImageStateNonModVisitorBase& v) const;

  virtual void ApplyIP(ImageStateModIPVisitorBase& v);

  virtual void ApplyIP(const ImageStateConstModIPVisitorBase& v);

  virtual ImageStateBasePtr Apply(ImageStateModIPVisitorBase& v) const;

  virtual ImageStateBasePtr Apply(const ImageStateConstModIPVisitorBase& v) const;

  virtual ImageStateBasePtr Apply(ImageStateModOPVisitorBase& v) const;

  virtual ImageStateBasePtr Apply(const ImageStateConstModOPVisitorBase& v) const;

  virtual ImageStateBasePtr Apply(ImageStateMorphVisitorBase& v);

  // operators
  virtual ImageStateBase& operator+=(Real v);
  virtual ImageStateBase& operator+=(const Complex& v);
  virtual ImageStateBase& operator-=(Real v);
  virtual ImageStateBase& operator-=(const Complex& v);
  virtual ImageStateBase& operator*=(Real v);
  virtual ImageStateBase& operator*=(const Complex& v);
  virtual ImageStateBase& operator/=(Real v);
  virtual ImageStateBase& operator/=(const Complex& v);

  virtual void operator+=(const Function& b);
  virtual void operator-=(const Function& b);
  virtual void operator*=(const Function& b);
  virtual void operator/=(const Function& b);

protected:
  // helper methods for operators above
  template <typename U> void add_scalar(const U& v);
  template <typename U> void mul_scalar(const U& v);

private:
  D domain_;
  ValueHolder<T> data_;
  PixelSampling sampling_;
  geom::Vec3 absolute_origin_;
};


}}} // ns

#endif
