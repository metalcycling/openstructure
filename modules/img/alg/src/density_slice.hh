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

#ifndef IMG_ALG_DENSITY_SLICE_HH
#define IMG_ALG_DENSITY_SLICE_HH

#include <ost/img/image_state.hh>
#include <ost/img/image_state/image_state_algorithm.hh>
#include <ost/img/alg/module_config.hh>


namespace ost { namespace img { namespace alg {

struct DLLEXPORT_IMG_ALG SliceFnc {
  SliceFnc();
  explicit SliceFnc(Real l,Real u);
  explicit SliceFnc(const Complex& l, const Complex& u);

  void SetThresholds(Real l,Real u); 
  void SetThresholds(const Complex& l, const Complex& u);
  
  template <typename T, class D>
  void VisitState(ImageStateImpl<T,D>& isi) const;

  static String GetAlgorithmName() {return "DensitySlice";}

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

typedef ImageStateConstModIPAlgorithm<SliceFnc> DensitySlice;

}

OST_IMG_ALG_EXPLICIT_INST_DECL(class,ImageStateConstModIPAlgorithm<alg::SliceFnc>)

}} // ns

#endif
