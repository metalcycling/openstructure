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

#include "extent_mask.hh"
#include "mask_visitor.hh"

namespace ost {
namespace img
{
ExtentMask::ExtentMask(const Extent& e):
  MaskBase(),
  Extent(e)
{
  
}

MaskPtr ExtentMask::Clone()
{
  return MaskPtr(new ExtentMask(*this));
}
bool ExtentMask::IsInside(const Vec2& v)
{
  return Contains(Point(static_cast<int>(floor(v[0])),static_cast<int>(floor(v[1]))));
}

void ExtentMask::Shift(const Vec2& v)
{
  Point p(v);
  Extent::Shift(p);
}
void ExtentMask::Expand(Real d)
{
  int id=static_cast<int>(d);
  *this=Extent(GetStart()-Point(id,id,id),GetEnd()+Point(id,id,id));
}

void ExtentMask::Scale(Real d)
{
  *this=Extent(Point(GetStart().ToVec2()*d),Point(GetEnd().ToVec2()*d));
}

void ExtentMask::Apply(MaskVisitor& v)
{
  v.VisitExtentMask(*this);
}

}}//ns
