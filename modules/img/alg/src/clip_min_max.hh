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

#ifndef CLIP_MIN_MAX_HH_
#define CLIP_MIN_MAX_HH_

namespace ost { namespace img { namespace alg {

namespace {

template<typename T>
bool cbigger(T& v, Real rv, const Complex& cv);

template<>
bool cbigger<Real>(Real& v, Real rv, const Complex& cv) {return v>=rv;}

template<>
bool cbigger<Complex>(Complex& v, Real rv, const Complex& cv) {return std::abs(v)>=std::abs(cv);}

template<typename T>
bool cbigger(T& v, Real rv, const Complex& cv)
{
  return static_cast<Real>(v)>=rv;
}
  template<typename T>
  T get_t(Real rv, const Complex& cv);
  template<>
  Real get_t<Real>(Real rv, const Complex& cv){return rv;}
  template<>
  Complex get_t<Complex>(Real rv, const Complex& cv){return cv;}
  template<typename T>
   T get_t(Real rv, const Complex& cv){return static_cast<T>(rv);}

struct ClipMinMaxFnc {
  explicit ClipMinMaxFnc() { set_rl(0.0); set_ru(1.0); }
  explicit ClipMinMaxFnc(Real l,Real u) {set_rl(l);set_ru(u);}
  explicit ClipMinMaxFnc(const Complex& l, const Complex& u) {set_cl(l);set_cu(u);}

  void SetThresholds(Real l,Real u) {set_rl(l);set_ru(u);}
  void SetThresholds(const Complex& l, const Complex& u) {set_cl(l);set_cu(u);}
  
  template <typename T, class D>
  void VisitState(ImageStateImpl<T,D>& isi) const {
    for(T* p = isi.Data().GetData(); p<isi.Data().GetEnd();++p) {
      (*p)= cbigger<T>(*p,r_tl_,c_tl_)? (cbigger<T>(*p,r_tu_,c_tu_)? get_t<T>(r_tu_,c_tu_) :(*p)) :get_t<T>(r_tl_,c_tl_) ;
    }
  }

  static String GetAlgorithmName() {return "ClipMinMax";}

private:
  Real r_tl_;
  Real r_tu_;
  Complex c_tl_;
  Complex c_tu_;

  void set_rl(Real f) {r_tl_=f; c_tl_=Val2Val<Real,Complex>(f);}
  void set_ru(Real f) {r_tu_=f; c_tu_=Val2Val<Real,Complex>(f);}
  void set_cl(const Complex& f) {c_tl_=f; r_tl_=Val2Val<Complex,Real>(f);}
  void set_cu(const Complex& f) {c_tu_=f; r_tu_=Val2Val<Complex,Real>(f);}

};

} // anon ns

typedef ImageStateConstModIPAlgorithm<ClipMinMaxFnc> ClipMinMax;

}

//no expl inst header only impl

  }} // ns



#endif /*SET_MIN_MAX_HH_*/
