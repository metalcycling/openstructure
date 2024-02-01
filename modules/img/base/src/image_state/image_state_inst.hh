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
  explicit instantiation of image state implementation

  Author: Ansgar Philippsen
*/

#ifndef IMG_IMAGE_STATE_INST_H
#define IMG_IMAGE_STATE_INST_H

/*
  At this point the various template instantiations
  in src code files are included
*/
#include "image_state_impl.cc"
#include "image_state_algorithm.cc"
#include "value_holder.cc"

#include "image_state_impl.hh"
#include "image_state_spatial_domain.hh"
#include "image_state_frequency_domain.hh"
#include "image_state_half_frequency_domain.hh"

namespace ost { namespace img { namespace image_state {

// explicit template instantiation so that all functionality is available
template class TEMPLATE_DEF_EXPORT ImageStateImpl<Real,SpatialDomain>;
template class TEMPLATE_DEF_EXPORT ImageStateImpl<Complex,SpatialDomain>;
template class TEMPLATE_DEF_EXPORT ImageStateImpl<Real,FrequencyDomain>;
template class TEMPLATE_DEF_EXPORT ImageStateImpl<Complex,FrequencyDomain>;
template class TEMPLATE_DEF_EXPORT ImageStateImpl<Complex,HalfFrequencyDomain>;
template class TEMPLATE_DEF_EXPORT ImageStateImpl<Word,SpatialDomain>;


template class TEMPLATE_DEF_EXPORT ValueHolder<Real>;
template class TEMPLATE_DEF_EXPORT ValueHolder<Complex>;
template class TEMPLATE_DEF_EXPORT ValueHolder<Word>;

}}} // namespaces

#endif
