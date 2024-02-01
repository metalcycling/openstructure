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

#ifndef IMG_PASTE_IMPL_HH
#define IMG_PASTE_IMPL_HH

#include "image_state/image_state_algorithm.hh"

namespace ost { namespace img { namespace detail {

/// \internal
class DLLEXPORT_OST_IMG_BASE PasteFnc {
public:
  PasteFnc();
  PasteFnc(const ImageHandle& target);

  template<typename V, class D>
  void VisitState(const image_state::ImageStateImpl<V,D>& isi);
  
  static String GetAlgorithmName() {return "Paste";}

private:
  ImageHandle target_;
};

/// \internal
typedef image_state::ImageStateNonModAlgorithm<PasteFnc> PasteImpl;

}

OST_IMG_BASE_EXPLICIT_INST_DECL(class,image_state::ImageStateNonModAlgorithm<detail::PasteFnc>)

}} // ns

#endif
