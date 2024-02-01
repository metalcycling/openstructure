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
  Author: Andreas Schenk
*/

#include <ost/img/image_state.hh>
#include <ost/img/alg/module_config.hh>

namespace ost { namespace img { namespace alg {

namespace {

template<typename T>
bool bigger(T& v, Real rv, const Complex& cv);

template<>
bool bigger<Real>(Real& v, Real rv, const Complex& cv) {return v>=rv;}

template<>
bool bigger<Complex>(Complex& v, Real rv, const Complex& cv) {return std::abs(v)>=std::abs(cv);}

template<typename T>
bool bigger(T& v, Real rv, const Complex& cv)
{
  return static_cast<Real>(v)>=rv;
}

struct ThresholdFnc {
  explicit ThresholdFnc(): r_fv_(0.0) { set_rval(0.0); }
  explicit ThresholdFnc(Real f) {set_rval(f);}
  explicit ThresholdFnc(const Complex& f) {set_cval(f);}

  void SetThreshold(Real f) {set_rval(f);}
  void SetThreshold(const Complex& f) {set_cval(f);}
  
  template <typename T, class D>
  void VisitState(ImageStateImpl<T,D>& isi) const {
    static const T zero=T(); // default ctor for zero value
    static const T one=T(1); // default ctor for one value
    for(T* p = isi.Data().GetData(); p<isi.Data().GetEnd();++p) {
      (*p)= bigger<T>(*p,r_fv_,c_fv_)? one : zero ;
    }
  }

  static String GetAlgorithmName() {return "Threshold";}

private:
  Real r_fv_;
  Complex c_fv_;

  void set_rval(Real f) {r_fv_=f; c_fv_=Val2Val<Real,Complex>(f);}
  void set_cval(const Complex& f) {c_fv_=f; r_fv_=Val2Val<Complex,Real>(f);}
};

} // anon ns

typedef ImageStateConstModIPAlgorithm<ThresholdFnc> Threshold;

}

OST_IMG_ALG_EXPLICIT_INST_DECL(class, ImageStateModIPAlgorithm<alg::ThresholdFnc>) 


}} // ns


