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
  top level include for all relevant image state classes

  Author: Ansgar Philippsen
*/

#ifndef IMG_IMAGE_STATE_H
#define IMG_IMAGE_STATE_H


#include "image_state/image_state_def.hh"
#include "image_state/image_state_factory.hh"
#include "image_state/index.hh"
#include "image_state/image_state_visitor.hh"
#include "image_state/image_state_algorithm.hh"
#include "image_state/image_state_base.hh"

namespace ost { namespace img {

// move essential names into img namespace

using image_state::ImageStateImpl;

using image_state::Index;
using image_state::IndexIterator;

using image_state::CreateState;
using image_state::ImageStateBase;

using image_state::RealSpatialImageState;
using image_state::ComplexSpatialImageState;
using image_state::RealFrequencyImageState;
using image_state::ComplexFrequencyImageState;
using image_state::ComplexHalfFrequencyImageState;

using image_state::SpatialDomain;
using image_state::FrequencyDomain;
using image_state::HalfFrequencyDomain;

using image_state::InvalidImageStateVisitor;
using image_state::ImageStateNonModVisitor;
using image_state::ImageStateModIPVisitor;
using image_state::ImageStateConstModIPVisitor;
using image_state::ImageStateModOPVisitor;
using image_state::ImageStateConstModOPVisitor;

using image_state::ImageStateNonModAlgorithm;
using image_state::ImageStateModIPAlgorithm;
using image_state::ImageStateConstModIPAlgorithm;
using image_state::ImageStateModOPAlgorithm;
using image_state::ImageStateConstModOPAlgorithm;
using image_state::ImageStateMorphAlgorithm;

using image_state::ImageStateBasePtr;

}}// namespace

#endif
