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

#ifndef IMG_ALG_FILL_HH
#define IMG_ALG_FILL_HH

#include <ost/img/image.hh>
#include <ost/img/image_state.hh>
#include <ost/img/alg/module_config.hh>

namespace ost { namespace img { namespace alg {

namespace {

template<typename T>
void set_val(T& v, Real rv, const Complex& cv);

/* While the functions below are employed by the OST unit tests, GCC does not
   notice this and complains about "unused functions".
*/
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
template<>
void set_val<Real>(Real& v, Real rv, const Complex& cv) {v=rv;}

template<>
void set_val<Complex>(Complex& v, Real rv, const Complex& cv) {v=cv;}
#pragma GCC diagnostic pop

template<typename T>
void set_val(T& v, Real rv, const Complex& cv)
{
  v=static_cast<T>(rv);
}

struct FillFnc {
  FillFnc() {}
   
  explicit FillFnc(Real f) {set_rval(f);}
  explicit FillFnc(const Complex& f) {set_cval(f);}

  void SetFillValue(Real f) {set_rval(f);}
  void SetFillValue(const Complex& f) {set_cval(f);}
  
  template <typename T, class D>
  void VisitState(ImageStateImpl<T,D>& isi) const {
    for(T* p = isi.Data().GetData(); p<isi.Data().GetEnd();++p) {
      /*
	set_val construct necessary due to missing 
	template function specialization...
      */
      set_val<T>(*p,r_fv_,c_fv_);
    }
  }

  static String GetAlgorithmName() {return "Fill";}

private:
  Real r_fv_;
  Complex c_fv_;

  void set_rval(Real f) {r_fv_=f; c_fv_=Val2Val<Real,Complex>(f);}
  void set_cval(const Complex& f) {c_fv_=f; r_fv_=Val2Val<Complex,Real>(f);}
};

} // anon ns

typedef ImageStateConstModIPAlgorithm<FillFnc> Fill;

}

// completely inlined definition does not need template instantiation

}} // ns



#endif
