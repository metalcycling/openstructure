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
#ifndef MASK_BASE_HH_
#define MASK_BASE_HH_

#include <ost/img/vecmat.hh>
#include "mask_base_fw.hh"
namespace ost {
namespace img
{

class MaskVisitor;

class DLLEXPORT_OST_IMG_BASE MaskBase
{
public:
  MaskBase();
  virtual ~MaskBase();
  virtual MaskPtr Clone()=0;
  virtual bool IsInside(const Vec2& v)=0;
  bool IsOutside(const Vec2& v){return !IsInside(v);}
  virtual void Shift(const Vec2& v)=0;
  virtual void Expand(Real d)=0;
  virtual void Shrink(Real d){Expand(-d);}
  virtual void Scale(Real d)=0;

  // traverse mask hierarchy with a visitor
  virtual void Apply(MaskVisitor& v)=0;

};

}}//ns

#endif /*MASK_BASE_HH_*/
