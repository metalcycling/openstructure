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

#include <ost/message.hh>

#include "fft.hh"
#include "auto_correlate.hh"

namespace ost { namespace img { namespace alg {

AutoCorrelateFnc::AutoCorrelateFnc()
{}

namespace {

template <typename T>
T abs_square(T val);

template <>
Complex abs_square<Complex>(Complex val) {return std::norm(val);}

template <typename T>
T abs_square(T val) {return val*val;}

}

template <typename V, class D>
ImageStateBasePtr AutoCorrelateFnc::VisitState(const ImageStateImpl<V,D>& isi) const
{
  ImageStateImpl<V,D> tmp(isi);
  
  for(ExtentIterator it(isi.GetExtent());!it.AtEnd();++it) {
    tmp.Value(it)=abs_square<V>(isi.Value(it));
  }

  ImageStateBasePtr fft_state_bp=FFT().VisitState(tmp);
  return fft_state_bp;
}

}

namespace image_state {
template class TEMPLATE_DEF_EXPORT ImageStateConstModOPAlgorithm<alg::AutoCorrelateFnc>;
}}} //
