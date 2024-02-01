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

#include "composite_mask.hh"
#include "mask_visitor.hh"

namespace ost { namespace img {


bool AndOp::operator()(const MaskPtr& lhs, const MaskPtr& rhs, const Vec2& v) const
{
  return lhs->IsInside(v) & rhs->IsInside(v);
}

bool OrOp::operator()(const MaskPtr& lhs, const MaskPtr& rhs, const Vec2& v) const
{
  return lhs->IsInside(v) | rhs->IsInside(v);
}

bool XorOp::operator()(const MaskPtr& lhs, const MaskPtr& rhs, const Vec2& v) const
{
  return lhs->IsInside(v) ^ rhs->IsInside(v);
}





CompositeMask::CompositeMask(const MaskPtr& lhs, const MaskPtr& rhs, const String& op):
  MaskBase(),
  lhs_(lhs),
  rhs_(rhs),
  op_()
{
  op_=create_op_(op);
}

CompositeMask::CompositeMask(const CompositeMask& m):
  MaskBase(),
  lhs_(m.lhs_),
  rhs_(m.rhs_),
  op_(m.op_->Clone())
{
}
CompositeMask::~CompositeMask()
{
  delete op_;
}
MaskPtr CompositeMask::Clone()
{
  return MaskPtr(new CompositeMask(*this));
}
bool CompositeMask::IsInside(const Vec2& v)
{
  return (*op_)(lhs_,rhs_,v);
}


void CompositeMask::Shift(const Vec2& v)
{
  lhs_->Shift(v);
  rhs_->Shift(v);
}
void CompositeMask::Expand(Real d)
{
  lhs_->Expand(d);
  rhs_->Expand(d);
}
void CompositeMask::Scale(Real d)
{
  lhs_->Scale(d);
  rhs_->Scale(d);
}

void CompositeMask::Apply(MaskVisitor& v)
{
  v.VisitCompositeMask(*this);
  lhs_->Apply(v);
  rhs_->Apply(v);
  v.VisitCompositeMaskFinish(*this);
}
OpBase* CompositeMask::create_op_(const String& name)
{
  if(name=="AND") return new AndOp;
  else if(name=="OR") return new OrOp;
  else if(name=="XOR") return new XorOp;
  else throw("Unknown Operator");
  return new AndOp;
}


}}//ns
