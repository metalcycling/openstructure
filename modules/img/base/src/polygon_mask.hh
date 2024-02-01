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

#ifndef POLYGON_MASK_HH_
#define POLYGON_MASK_HH_

#include <ost/img/vecmat.hh>
#include "mask_base.hh"

namespace ost {

namespace img
{
///
/// \sa \ref create_split_image.py "Create Split Image Example"
class DLLEXPORT_OST_IMG_BASE PolygonMask : public MaskBase, public Polygon2
{
public:
  PolygonMask();
  PolygonMask(const Polygon2& p);
  virtual MaskPtr Clone();
  virtual bool IsInside(const Vec2& v);

  virtual void Shift(const Vec2& v);
  virtual void Expand(Real d);
  virtual void Scale(Real d);

  virtual void Apply(MaskVisitor& v);

  void AddNode(const Vec2& n);
  void SetNode(unsigned int i,const Vec2& v);
  void Erase(circular_iterator first,circular_iterator last);

  
protected:
  std::vector<Real> get_intercept_x_values_(Real y);
  bool intercepts_valid_;
  Real intercept_y_;
  std::vector<Real> intercepts_;
  Rectangle2 bounding_box_;
  

};

}}//ns

#endif /*POLYGON_MASK_HH_*/
