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
  typedefs for all 'official' image states

  Author: Ansgar Philippsen
*/

#ifndef IMG_IMAGE_STATE_DEF_H
#define IMG_IMAGE_STATE_DEF_H

#include <ost/base.hh>

#include "image_state_impl.hh"
#include "image_state_spatial_domain.hh"
#include "image_state_frequency_domain.hh"
#include "image_state_half_frequency_domain.hh"

namespace ost { namespace img { namespace image_state {

typedef ImageStateImpl<Real, SpatialDomain> RealSpatialImageState;
typedef ImageStateImpl<Real, FrequencyDomain> RealFrequencyImageState;
typedef ImageStateImpl<Complex, SpatialDomain> ComplexSpatialImageState;
typedef ImageStateImpl<Complex, FrequencyDomain> ComplexFrequencyImageState;
typedef ImageStateImpl<Complex, HalfFrequencyDomain> ComplexHalfFrequencyImageState;
typedef ImageStateImpl<Word, SpatialDomain> WordSpatialImageState;


} // namespace image_state

#if defined(OST_MODULE_OST_IMG)
#  if defined(_MSC_VER)
#    define OST_DECL_IMAGE_STATE(c, v, d)
#  else
#    define OST_DECL_IMAGE_STATE(c, v, d) extern template c TEMPLATE_EXPORT image_state::ImageStateImpl<v, image_state::d>;
#  endif
#else
#  define OST_DECL_IMAGE_STATE(c, v, d) extern template c DLLIMPORT image_state::ImageStateImpl<v, image_state::d>;
#endif

OST_DECL_IMAGE_STATE(class, Real, SpatialDomain)
OST_DECL_IMAGE_STATE(class, Real, FrequencyDomain)
OST_DECL_IMAGE_STATE(class, Complex, SpatialDomain)
OST_DECL_IMAGE_STATE(class, Complex, FrequencyDomain)
OST_DECL_IMAGE_STATE(class, Complex, HalfFrequencyDomain)
OST_DECL_IMAGE_STATE(class, Word, SpatialDomain)

#undef OST_DECL_IMAGE_STATE

}} // namespaces

#endif
