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

#include "circle_mask.hh"
#include "mask_visitor.hh"

namespace ost { namespace img {

CircleMask::CircleMask():
  Circle2()
{
}
CircleMask::CircleMask(const Circle2& c):
  Circle2(c)
{
}
CircleMask::CircleMask(const Vec2& v, Real r):
  Circle2(v,r)
{
}

CircleMask::~CircleMask()
{
}

MaskPtr CircleMask::Clone()
{
  return MaskPtr(new CircleMask(*this));

}
bool CircleMask::IsInside(const Vec2& v)
{
  return IsInCircle(*this,Vec2(v));
}
void CircleMask::Shift(const Vec2& v)
{
  center_+=Vec2(v);
}
void CircleMask::Expand(Real d)
{
  radius_+=d;
}
void CircleMask::Scale(Real d)
{
  center_*=d;
  radius_*=d;  
}

void CircleMask::Apply(MaskVisitor& v)
{
  v.VisitCircleMask(*this);
}


}}//ns
