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
  "overloaded" transcendentals for images; instead of providing the
  typical func(image) interface, which would limit the use to out-of-place
  application, a full fletched img::alg interface is given for each:

  cos, exp, log, log10, pow, sin, sqrt, tan

  Author: Ansgar Philippsen
*/

#ifndef IMG_ALG_TRANS_HH
#define IMG_ALG_TRANS_HH

#include <ost/img/image.hh>
#include <ost/img/image_state.hh>
#include <ost/img/alg/module_config.hh>


#define IMG_ALG_TRANSCENDENTALS_BLOCK(FF,NN,SS)                                \
struct FF {					                                                   \
  FF() {}                                                                      \
  ~FF() {}                                                                     \
  template <typename T, class D>                                               \
  void VisitState(ImageStateImpl<T,D>& isi) const {	                           \
    const T* end = isi.Data().GetEnd();                                        \
    for(T* it = isi.Data().GetData(); it!=end; ++it) {                         \
      (*it) = SS (*it);                                                        \
    }                                                                          \
  }                                                                            \
  template <class D>                                                           \
  void VisitState(ImageStateImpl<Word,D>& isi) const {                         \
    const Word* end = isi.Data().GetEnd();                                     \
    for(Word* it = isi.Data().GetData(); it!=end; ++it) {                      \
      (*it) = static_cast<Word>(SS(static_cast<Real>(*it)));                   \
    }                                                                          \
  }                                                                            \
  static String GetAlgorithmName() {return "";}                                \
};                                                                             \
typedef ImageStateConstModIPAlgorithm<FF> NN;       

namespace ost { namespace img { namespace alg {

IMG_ALG_TRANSCENDENTALS_BLOCK(CosFnc,Cos,std::cos)
IMG_ALG_TRANSCENDENTALS_BLOCK(ExpFnc,Exp,std::exp)
IMG_ALG_TRANSCENDENTALS_BLOCK(LogFnc,Log,std::log)
IMG_ALG_TRANSCENDENTALS_BLOCK(Log10Fnc,Log10,std::log)
IMG_ALG_TRANSCENDENTALS_BLOCK(SinFnc,Sin,std::sin)
IMG_ALG_TRANSCENDENTALS_BLOCK(SqrtFnc,Sqrt,std::sqrt)
IMG_ALG_TRANSCENDENTALS_BLOCK(TanFnc,Tan,std::tan)


struct PowFnc {
  PowFnc():exp_(1.0) {}
  PowFnc(Real exp): exp_(exp) {}

  template <typename T, class D>
  void VisitState(ImageStateImpl<T,D>& isi) const {
    const T* end = isi.Data().GetEnd();
    for(T* it = isi.Data().GetData(); it!=end; ++it) {
      (*it) = std::pow(*it,exp_);
    }
  }
  static String GetAlgorithmName() {return "";}
 private:
  Real exp_;
};

typedef ImageStateConstModIPAlgorithm<PowFnc> Pow;           


}}} // ns

#undef IMG_ALG_TRANSCENDENTALS_BLOCK


#endif
