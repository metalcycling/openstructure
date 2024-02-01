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
#include "density_slice.hh"

namespace ost { namespace img { namespace alg {

template<typename T>
bool between(const T& v, Real rl, const Complex& cl, Real ru, const Complex& cu);

template<>
bool between<Real>(const Real& v, Real rl, const Complex& cl, Real ru, const Complex& cu) {return v>=rl && v<=ru ;}

template<>
bool between<Complex>(const Complex& v, Real rl, const Complex& cl, Real ru, const Complex& cu) {return std::abs(v)>=std::abs(cl) && std::abs(v)<=std::abs(cu);}

template<typename T>
bool between(const T& v, Real rl, const Complex& cl, Real ru, const Complex& cu)
{
  return between<Real>(static_cast<Real>(v),rl,cl,ru,cu);
}

SliceFnc::SliceFnc() {set_rl(0.0); set_ru(1.0);}
SliceFnc::SliceFnc(Real l,Real u) {set_rl(l);set_ru(u);}
SliceFnc::SliceFnc(const Complex& l, const Complex& u) {set_cl(l);set_cu(u);}

void SliceFnc::SetThresholds(Real l,Real u) {set_rl(l);set_ru(u);}
void SliceFnc::SetThresholds(const Complex& l, const Complex& u) {set_cl(l);set_cu(u);}
  
template <typename T, class D>
void SliceFnc::VisitState(ImageStateImpl<T,D>& isi) const 
{
  static const T zero=T(); // default ctor for zero value
  static const T one=T(1); // default ctor for one value
  for(T* p = isi.Data().GetData(); p<isi.Data().GetEnd();++p) {
    (*p)= between<T>(*p,r_tl_,c_tl_,r_tu_,c_tu_)? one : zero ;
  }
}

}
template class TEMPLATE_DEF_EXPORT image_state::ImageStateConstModIPAlgorithm<alg::SliceFnc>;
}} // ns
