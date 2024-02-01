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
  Authors: Andreas Schenk, Ansgar Philippsen
*/

#ifndef IMG_MASK_HH_
#define IMG_MASK_HH_

#include <ost/message.hh> 
#include "polygon_mask.hh" 
#include "extent_mask.hh" 
#include "circle_mask.hh" 
#include "spherical_mask.hh" 
#include "mask_op.hh" 

namespace ost { namespace img {


MaskPtr DLLEXPORT_OST_IMG_BASE Mask(const Polygon2& p);
MaskPtr DLLEXPORT_OST_IMG_BASE Mask(const Extent& e);
MaskPtr DLLEXPORT_OST_IMG_BASE Mask(const Circle2& c);
MaskPtr DLLEXPORT_OST_IMG_BASE Mask(const Sphere& s);


}}//ns

#endif /*IPLT_MASK_HH_*/
