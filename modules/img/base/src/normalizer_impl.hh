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
  normalizer interface

  Authors: Ansgar Philippsen, Andreas Schenk
*/

#ifndef IMG_NORMALIZER_IMPL_H
#define IMG_NORMALIZER_IMPL_H

#include <boost/shared_ptr.hpp>

#include <ost/img/image_state.hh>
#include <ost/img/value_util.hh>
#include "normalizer_fw.hh"
#include <ost/img/module_config.hh>
namespace ost { namespace img {

//! Abstract normalizer base class
// \internal
/*!
  A normalizer allows the values in a data set to be
  normalized with a specific formula, without encoding
  this functionality within the data class itself.
*/
class DLLEXPORT_OST_IMG_BASE NormalizeFnc {
public:
  virtual ~NormalizeFnc() {}
  template <typename T, class D>
  void VisitState(ImageStateImpl<T,D>& s) const {

    for(T* ptr = s.Data().GetData(); ptr<s.Data().GetEnd(); ++ptr) {
      *ptr = static_cast<T>(Convert(*ptr));
    }
  }

  // normalizer abstract interface
  virtual Real BackConvert(Real v) const = 0;
  virtual Complex BackConvert(Complex v) const = 0;
  virtual Real Convert(Real v) const = 0;
  virtual Complex Convert(Complex v) const = 0;

  static String GetAlgorithmName() {return "Normalizer";}
};

/// \internal
typedef ImageStateConstModIPAlgorithm<NormalizeFnc> NormalizerBase;

/// \internal
class NormalizerImpl: public NormalizerBase{
public:
  virtual NormalizerImpl* Clone() const =0;
};

/// \internal
class DLLEXPORT_OST_IMG_BASE NoOpNormalizer: public NormalizerImpl {
public:
  // normalizer abstract interface
  virtual Real BackConvert(Real v) const  {return v;} 
  virtual Complex BackConvert(Complex v) const  {return v;}
  virtual Real Convert(Real v) const {return v;}
  virtual Complex Convert(Complex v) const {return v;}
  virtual NormalizerImpl* Clone() const {return new NoOpNormalizer((*this));}
};

/// \internal
class RangeHandler
{
public:
  RangeHandler():min_(),max_(){}
  RangeHandler(Real min,Real max):min_(min),max_(max){}
  Real GetMinimum() const{return min_;}
  Real GetMaximum() const{return max_;}
protected:
  Real Clip(Real val) const{return ClipMinMax(val,min_,max_);}  
  Complex Clip(Complex val) const{return ClipMinMax(val,min_,max_);}  
  Real ClipMinMax(Real val,Real minval,Real maxval) const{return std::max<Real>(std::min<Real>(val,maxval),minval);}  
  Complex ClipMinMax(Complex val,Real minval,Real maxval) const {
      return (val == Complex(0.0)? 1.0 : val/std::abs(val))*std::max<Real>(std::min<Real>(std::abs(val),maxval),minval);  
  }
private:
 Real min_;
 Real max_;
};

}} // namespace

#endif
