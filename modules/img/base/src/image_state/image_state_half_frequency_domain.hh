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
  Author: Ansgar Philippsen
*/

#ifndef IMAGE_STATE_HALF_FREQUENCY_DOMAIN_HH
#define IMAGE_STATE_HALF_FREQUENCY_DOMAIN_HH

#include <ost/base.hh>
#include <ost/img/point.hh>
#include <ost/img/size.hh>
#include <ost/img/extent.hh>
#include <ost/img/pixel_sampling.hh>

#include "value_holder.hh"
#include "index.hh"

#define CONJPOINT(a)  ((a)[2]==0 ?((a)[1]==0?(a)[0]>0:(a)[1]>0):(a)[2]>0) ? (a) : -(a)

namespace ost { namespace img { namespace image_state {

namespace {

int half_pos(int d) 
{
  return (d&0x1) ? (d-1)/2 : d/2;
}

int half_neg(int d) 
{
  return (d&0x1) ? -((d-1)/2) : -(d/2-1);
}

} // anon ns

/*
  1D half-frequency domain

  N  |   logical range  |  physical range
 ----------------------------------------
 even  -(N/2-1) : N/2       0:N/2
 odd   -(N-1)/2 : (N-1)/2   0:N/2+1
 

  2D half-frequency domain

  N is height 
  M is width

  N  |   logical range  |  physical range
 ----------------------------------------
 even  -(N/2-1) : N/2       0:N/2
 odd   -(N-1)/2 : (N-1)/2   0:N/2+1

  M  |   logical range  |  physical range
 ----------------------------------------
 even  -(M/2-1) : M/2       0:M-1
 odd   -(M-1)/2 : (M-1)/2   0:M-1

 the negative indices for the width are wrapped around, so
 that -1 == M-1, -2 == M-2 etc.


  3D half-frequency domain

  N is depth
  M is height
  L is width

  N  |   logical range  |  physical range
 ----------------------------------------
 even  -(N/2-1) : N/2       0:N/2
 odd   -(N-1)/2 : (N-1)/2   0:N/2+1

 M/L  |   logical range  |  physical range
 ----------------------------------------
 even  -(M/2-1) : M/2       0:M-1
 odd   -(M-1)/2 : (M-1)/2   0:M-1

 the negative indices for the width and height are wrapped around, so
 that -1 == M-1, -2 == M-2 etc.

*/

class DLLEXPORT_OST_IMG_BASE HalfFrequencyDomain {
public:
  HalfFrequencyDomain(const Extent& e):
    logical_extent_(ExtractLogicalExtent(e)),
    physical_extent_(ExtractPhysicalExtent(e)),
    dim_(e.GetDim()),
    spat_ori_()
  {
    assert(dim_>=1 && dim_<=3);
  }

  DataDomain GetDomain() const {return HALF_FREQUENCY;}

  void SetSpatialOrigin(const Point& o) {spat_ori_=o;}

  Point GetSpatialOrigin() const {return spat_ori_;}

  // this is the physical extent
  Extent GetExtent() const {return physical_extent_;}

  Extent GetLogicalExtent() const { return logical_extent_;}

  Extent GetPhysicalExtent() const { return physical_extent_;}

  template <typename V>
  Real GetReal(const Point &p, const ValueHolder<V>& data) const {
    if(logical_extent_.Contains(p)) {
      return Val2Val<V,Real>(data.Value(Point2Index(CONJPOINT(p))));
    }
    return Real();
  }

  template <typename V>
  void SetReal(const Point &p, const Real& r, ValueHolder<V>& data) {
    if(logical_extent_.Contains(p)) {
      Point cp=CONJPOINT(p);
      if(cp==p){
	data.Value(Point2Index(cp))=Val2Val<Real,V>(r);
      }
    }
  }

  template <typename V>
  Complex GetComplex(const Point &p, const ValueHolder<V>& data) const {
    if(logical_extent_.Contains(p)) {
      Point cp=CONJPOINT(p);
      if(cp==p){
	return Val2Val<V,Complex>(data.Value(Point2Index(cp)));
      }else{
	return conj(Val2Val<V,Complex>(data.Value(Point2Index(cp))));
      }
    }
    return Complex();
  }

  template <typename V>
  void SetComplex(const Point &p, const Complex& c, ValueHolder<V>& data) {
    if(logical_extent_.Contains(p)) {
      Point cp=CONJPOINT(p);
      if(cp==p){
	data.Value(Point2Index(cp))=Val2Val<Complex,V>(c);
      }
    }
  }

  Index Point2Index(const Point& p) const {
    const Size& size=physical_extent_.GetSize();
    return Index(p[0]>=0?p[0]:size[0]+p[0],
                 p[1]>=0?p[1]:size[1]+p[1],
                 p[2]);
  }

  Extent ExtractLogicalExtent(const Extent& e) const {
#if 1
      return Extent(Point(half_neg(e.GetWidth()),half_neg(e.GetHeight()),half_neg(e.GetDepth())),
		    Point(half_pos(e.GetWidth()),half_pos(e.GetHeight()),half_pos(e.GetDepth())));
#else

    if(e.GetDim()==1) {
      return Extent(Point(half_neg(e.GetWidth())),Point(half_pos(e.GetWidth())));
    } else if(e.GetDim()==2) {
      return Extent(Point(half_neg(e.GetWidth()),half_neg(e.GetHeight())),
		    Point(half_pos(e.GetWidth()),half_pos(e.GetHeight())));
    } else {
      return Extent(Point(half_neg(e.GetWidth()),half_neg(e.GetHeight()),half_neg(e.GetDepth())),
		    Point(half_pos(e.GetWidth()),half_pos(e.GetHeight()),half_pos(e.GetDepth())));
    }
#endif
  }

  Extent ExtractPhysicalExtent(const Extent& e) const {
    if(e.GetDim()==1) {
      return Extent(Point(0),Point(half_pos(e.GetWidth())));
    } else if(e.GetDim()==2) {
      return Extent(Point(half_neg(e.GetWidth()),0,0),
		    Point(half_pos(e.GetWidth()),half_pos(e.GetHeight()),0));
    } else {
      return Extent(Point(half_neg(e.GetWidth()),half_neg(e.GetHeight()),0),
		    Point(half_pos(e.GetWidth()),half_pos(e.GetHeight()),half_pos(e.GetDepth())));
    }
  }

private:
  Extent logical_extent_;
  Extent physical_extent_;
  int dim_;
  Point spat_ori_;

};



}}} // ns

#endif
