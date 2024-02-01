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

#include <algorithm>
#include "polygon_mask.hh"
#include "mask_visitor.hh"

namespace ost { namespace img {

PolygonMask::PolygonMask():
  MaskBase(),
  Polygon2(),
  intercepts_valid_(false),
  intercept_y_(),
  intercepts_(),
  bounding_box_()
{}

PolygonMask::PolygonMask(const Polygon2& p):
  MaskBase(),
  Polygon2(p),
  intercepts_valid_(false),
  intercept_y_(),
  intercepts_(),
  bounding_box_(p.GetBoundingBox())
{}

MaskPtr PolygonMask::Clone()
{
  return MaskPtr(new PolygonMask(*this));
}

bool PolygonMask::IsInside(const Vec2& v)
{
  if (GetNodeCount()<3 || !IsInRectangle(bounding_box_,v))
    return false;
  if(!intercepts_valid_ || v[1]!=intercept_y_){
    intercepts_=get_intercept_x_values_(v[1]);
    intercepts_valid_=true;
    intercept_y_=v[1];
  }
  unsigned int count=0;
  for(std::vector<Real>::const_iterator it=intercepts_.begin();it!=intercepts_.end();++it){
    if(v[0]>(*it)){
      ++count;
    }
  }
  return count%2==1;
}

void PolygonMask::Shift(const Vec2& v)
{
  intercepts_valid_=false;
  for(unsigned int i=0; i<GetNodeCount();++i){
    SetNode(i,Vec2(GetNode(i)+v));
  }
  bounding_box_=GetBoundingBox();
}

void PolygonMask::Expand(Real d)
{
  intercepts_valid_=false;
  Polygon2::Expand(d);
  bounding_box_=GetBoundingBox();
}
void PolygonMask::Scale(Real d)
{
  intercepts_valid_=false;
  for(unsigned int i=0; i<GetNodeCount();++i){
    SetNode(i,Vec2(GetNode(i)*d));
  }
  bounding_box_=GetBoundingBox();
}
void PolygonMask::AddNode(const Vec2& n)
{
  intercepts_valid_=false;
  Polygon2::AddNode(n);
  bounding_box_=GetBoundingBox();
}
void PolygonMask::SetNode(unsigned int i,const Vec2& v)
{
  intercepts_valid_=false;
  Polygon2::SetNode(i,v);
  bounding_box_=GetBoundingBox();
}
void PolygonMask::Erase(circular_iterator first,circular_iterator last)
{
  intercepts_valid_=false;
  Polygon2::Erase(first,last);
  bounding_box_=GetBoundingBox();
}

void PolygonMask::Apply(MaskVisitor& v)
{
  v.VisitPolygonMask(*this);
}

std::vector<Real> PolygonMask::get_intercept_x_values_(Real y)
{
  std::vector<Real> result;
  Vec2 old=GetNode(GetNodeCount()-1);
  bool oldcomp= y<old[1];
  for (const_iterator it=begin();it!=end(); ++it) {
    bool comp= y < (*it)[1];
    if( comp ^ oldcomp ) {
      result.push_back((old[0] - (*it)[0]) * (y - (*it)[1]) / (old[1] - (*it)[1]) + (*it)[0]); 
    }
    oldcomp=comp;
    old=(*it);
  }
  sort(result.begin(),result.end());
  return result;
}


}}//ns

