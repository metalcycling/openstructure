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
  Author: Ansgar Philippsen
*/

#ifndef MASK_VISITOR_HH_
#define MASK_VISITOR_HH_

#include "circle_mask.hh"
#include "composite_mask.hh"
#include "extent_mask.hh"
#include "inverted_mask.hh"
#include "polygon_mask.hh"
#include "spherical_mask.hh"

namespace ost { namespace img {

class DLLEXPORT_OST_IMG_BASE MaskVisitor
{
public:
  virtual ~MaskVisitor() {}
  virtual void VisitCircleMask(CircleMask& m) {}
  virtual void VisitCompositeMask(CompositeMask& m) {}
  virtual void VisitCompositeMaskFinish(CompositeMask& m) {}
  virtual void VisitExtentMask(ExtentMask& m) {}
  virtual void VisitInvertedMask(InvertedMask& m) {}
  virtual void VisitInvertedMaskFinish(InvertedMask& m) {}
  virtual void VisitPolygonMask(PolygonMask& m) {}
  virtual void VisitSphericalMask(SphericalMask& m) {}
};

}} // ns

#endif

