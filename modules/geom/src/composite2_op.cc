//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
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
#include <cmath>
#include "exc.hh"
#include "vec2.hh"
#include "mat2.hh"
#include "vecmat2_op.hh"
#include "composite2_op.hh"

namespace geom {


Real Angle(const Line2& l1, const Line2& l2)
{
  return Angle(l1.GetDirection(),l2.GetDirection());
}


Vec2 IntersectionPoint(const Line2& l1, const Line2& l2)
{
  if(AreParallel(l1,l2)) throw GeomException("Lines are parallel!");
  Vec2 p1=l1.At(1);
  Vec2 p2=l1.GetOrigin();
  Vec2 p3=l2.At(1);
  Vec2 p4=l2.GetOrigin();
  Vec2 dir1=l1.GetDirection();
  Vec2 dir2=l2.GetDirection();
  Real d1=Det(Mat2(p1[0],p1[1],p2[0],p2[1]));
  Real d2=Det(Mat2(p3[0],p3[1],p4[0],p4[1]));
  Real d3=Det(Mat2(dir1[0],dir1[1],dir2[0],dir2[1]));
  return Vec2(Det(Mat2(d1,dir1[0],d2,dir2[0]))/d3,Det(Mat2(d1,dir1[1],d2,dir2[1]))/d3);
  
}

Real Distance(const Line2& l, const Vec2& v)
{
  Vec2 d1=v-l.GetOrigin();
  Vec2 d2=l.GetDirection();
  return std::fabs(Det(Mat2(d2[0],d1[0],d2[1],d1[1]))/Length(d1));
}

bool IsOnLine(const Line2& l, const Vec2& v, Real ephilon)
{
  Vec2 d1=v-l.GetOrigin();
  Vec2 d2=l.GetDirection();
  return Length2(d1-Dot(d1, d2)*d2)<(ephilon*ephilon);
}


bool operator==(const Line2& l1, const Line2& l2)
{
  return Equal(l1,l2);
}
bool Equal(const Line2& l1, const Line2& l2,Real ephilon)
{
  return  Equal(l1.GetDirection(),l2.GetDirection(),ephilon) && Equal(l1.GetOrigin(),l2.GetOrigin(),ephilon);
}
bool EqualPosition(const Line2& l1, const Line2& l2,Real ephilon)
{
  return  IsOnLine(l1,l2.GetOrigin()) && AreParallel(l1,l2);
}


bool AreParallel(const Line2& l1, const Line2& l2,Real ephilon)
{
  return std::fabs(1-Dot(l1.GetDirection(),l1.GetDirection()))<ephilon;  
}

bool IsInPolygon(const Polygon2& p, const Vec2& v)
{
  if(p.GetNodeCount()<3){
    return false;
  }
  bool c = false;
  Vec2 old=p.GetNode(p.GetNodeCount()-1);
  bool oldcomp= v.y<old.y;
  for (Polygon2::const_iterator it=p.begin();it!=p.end(); ++it) {
    bool comp= v.y < (*it).y;
    if( comp ^ oldcomp ) {
      if( (v.x < (old.x - (*it).x) * (v.y - (*it).y) / (old.y - (*it).y) + (*it).x) ){
        c = !c;
      }
    }
  oldcomp=comp;
        old=(*it);
  }
  return c;
}
DLLEXPORT Polygon2 SimplifyPolygon(const Polygon2& p)
{
  Polygon2 result;
  Polygon2::const_circular_iterator it1=p.cend(),it2=p.cend(),it3=p.cbegin();
  --it1;
  for(unsigned int i=0;i<p.GetNodeCount();++i){
    if(!IsOnLine(Line2(*it1,*it3),*it2)){
      result.AddNode(*it2);
    }
    ++it1;
    ++it2;
    ++it3;
  }
  return result;
}



DLLEXPORT bool IsInRectangle(const Rectangle2& r, const Vec2& v)
{
  Vec2 start=r.GetStart();
  Vec2 end=r.GetEnd();
  return v.x>=start.x && v.x<=end.x &&  v.y>=start.y && v.y<=end.y;
}
DLLEXPORT bool IsInCircle(const Circle2& c, const Vec2& v)
{
  return Length(v-c.GetCenter())<=c.GetRadius();
}

} // ns

