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

#include "spherical_mask.hh"
#include "mask_visitor.hh"

namespace ost {

namespace img
{

SphericalMask::SphericalMask():
  MaskBase(),
  Sphere()
{
}

SphericalMask::~SphericalMask()
{
}



SphericalMask::SphericalMask(const Sphere& s):
  MaskBase(),
  Sphere(s)
{
}

SphericalMask::SphericalMask(const Vec2& v, Real d):
  MaskBase(),
  Sphere(v,d)
{
}

MaskPtr SphericalMask::Clone()
{
  return MaskPtr(new SphericalMask(*this));
}
bool SphericalMask::IsInside(const Vec2& v)
{
  return IsInSphere(*this,v);
}
void SphericalMask::Shift(const Vec2& v)
{
  SetOrigin(GetOrigin()+v);
}
void SphericalMask::Expand(Real d)
{
  SetRadius(GetRadius()+d);
}
void SphericalMask::Scale(Real d)
{
  SetOrigin(GetOrigin()*d);
  SetRadius(GetRadius()*d);
}

void SphericalMask::Apply(MaskVisitor& v)
{
  v.VisitSphericalMask(*this);
}

}}//ns
