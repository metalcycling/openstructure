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

#ifndef IMG_ALG_FFT_HH
#define IMG_ALG_FFT_HH

#include <ost/message.hh>

#include <ost/img/alg/module_config.hh>
#include <ost/img/image_state.hh>
#include <ost/img/image_state/image_state_algorithm.hh>
#include <ost/img/value_util.hh>

namespace ost { namespace img { namespace alg {

struct DLLEXPORT FFTException: public Error {
  FFTException(const String& m): Error(m) {}
};

struct DLLEXPORT_IMG_ALG FFTFnc {
  FFTFnc();
  FFTFnc(bool);
  ~FFTFnc();
  template <typename T, class D>
  ImageStateBasePtr VisitState(const ImageStateImpl<T,D>& s) const;

  static String GetAlgorithmName() {return "FFT";}
private:
  bool ori_flag_;
};

typedef image_state::ImageStateConstModOPAlgorithm<FFTFnc> FFT;

}


}} // ns

OST_IMG_ALG_EXPLICIT_INST_DECL(class, ImageStateConstModOPAlgorithm<ost::img::alg::FFTFnc>)


#endif

