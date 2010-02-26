//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2010 by the OpenStructure authors
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
  convenience fw for all relevant image state classes

  Author: Ansgar Philippsen 
*/

#ifndef IPLT_IMAGE_STATE_FW_H
#define IPLT_IMAGE_STATE_FW_H

#include "image_state/image_state_factory.hh"
#include "image_state/image_state_visitor_fw.hh"
#include "image_state/image_state_base_fw.hh"

namespace ost { namespace iplt {

namespace image_state {

struct Index;
class IndexIterator;
class ImageStateBase;

}

// move essential names into iplt namespace

using image_state::Index;
using image_state::IndexIterator;

using image_state::CreateState;
using image_state::ImageStateBase;

using image_state::ImageStateBasePtr;

}} // namespace

#endif
