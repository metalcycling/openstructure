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
  Author: Andreas Schenk
*/
#include "inverted_mask.hh"
#include "mask_visitor.hh"

namespace ost {

namespace img
{

InvertedMask::InvertedMask(const MaskPtr& lhs):
  MaskBase(),
  lhs_(lhs)
{
}
InvertedMask::InvertedMask(const InvertedMask& m):
  MaskBase(),
  lhs_(m.lhs_)
{
}

InvertedMask::~InvertedMask()
{
}
MaskPtr InvertedMask::Clone()
{
  return MaskPtr(new InvertedMask(*this));
}
bool InvertedMask::IsInside(const Vec2& v)
{
  return ! lhs_->IsInside(v);
}

void InvertedMask::Apply(MaskVisitor& v)
{
  v.VisitInvertedMask(*this);
  lhs_->Apply(v);
  v.VisitInvertedMaskFinish(*this);
}


void InvertedMask::Shift(const Vec2& v)
{
  lhs_->Shift(v);
}
void InvertedMask::Expand(Real d)
{
  lhs_->Expand(d);
}
void InvertedMask::Scale(Real d)
{
  lhs_->Scale(d);
}

}}//ns
