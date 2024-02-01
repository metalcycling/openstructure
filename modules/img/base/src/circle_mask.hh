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

#include <ost/img/vecmat.hh>
#include "mask_base.hh"

#ifndef CIRCLE_MASK_HH_
#define CIRCLE_MASK_HH_



namespace ost { namespace img {

class CircleMask: public MaskBase, public Circle2
{
public:
  CircleMask();
  CircleMask(const Circle2& c);
  CircleMask(const Vec2& v, Real d);
  virtual ~CircleMask();
	
  virtual MaskPtr Clone();

  virtual bool IsInside(const Vec2& v);
  virtual void Shift(const Vec2& v);
  virtual void Expand(Real d);
  virtual void Scale(Real d);

  // traverse mask hierarchy with a visitor
  virtual void Apply(MaskVisitor& v);

};

}} // ns

#endif /*CIRCLE_MASK_HH_*/
