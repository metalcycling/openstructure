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

#ifndef EXTENT_MASK_HH_
#define EXTENT_MASK_HH_

#include "extent.hh"
#include "mask_base.hh"

namespace ost {
namespace img
{

class DLLEXPORT_OST_IMG_BASE ExtentMask : public MaskBase, public Extent
{
public:
  ExtentMask(const Extent& e);
  virtual MaskPtr Clone();
  virtual bool IsInside(const Vec2& v);
  virtual void Shift(const Vec2& v);
  virtual void Expand(Real d);
  virtual void Scale(Real d);

  // traverse mask hierarchy with a visitor
  virtual void Apply(MaskVisitor& v);
};


}} // ns

#endif /*EXTENT_MASK_HH_*/
