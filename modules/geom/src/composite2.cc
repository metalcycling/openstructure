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
/*
  Author: Andreas Schenk
*/

#include "composite2.hh"
#include "composite2_op.hh"
#include "exc.hh"

#include "vecmat2_op.hh"
namespace geom {

/************************************************************
*  Line2                                                    *
*************************************************************/
Line2::Line2():ori_(Vec2(0,0)),dir_(Vec2(1,0)) {};

Line2::Line2(const Vec2& from, const Vec2& to):ori_(from)
{
        dir_=(to-from)/Length(to-from);
}

Vec2 Line2::At(Real r) const   {return ori_+r*dir_;};
Vec2 Line2::GetOrigin() const    {return ori_;};
Vec2 Line2::GetDirection() const {return dir_;};

/************************************************************
*  Ellipse2                                                  *
*************************************************************/
Ellipse2::Ellipse2(Vec2 ori, Real a, Real b, Real gamma):
  origin_(ori),
  a_(a),
  b_(b),
  gamma_(gamma)
{
}
Ellipse2::Ellipse2():
  origin_(Vec2(0.0,0.0)),
  a_(1.0),
  b_(1.0),
  gamma_(0.0)
{
}

Vec2 Ellipse2::AtAngle(Real angle) const
{
    Real t;
    angle=fmod(angle+M_PI,2*M_PI)-M_PI;
    if(a_==0.0 || b_==0.0)
        t=0.0;
    else
        t=atan2(a_*(sin(gamma_)-tan(angle)*cos(gamma_)),-b_*(sin(gamma_)*tan(angle)+cos(gamma_)));
    if(angle<-M_PI*0.5){t-=M_PI;}
    if(angle>M_PI*0.5){t+=M_PI;}
    t+=M_PI;
    return calc_(t);
}
Vec2 Ellipse2::At(Real t) const
{
    Real t2;
    if(a_==0.0 || b_==0.0)
        t2=0.0;
    else
        t2=t-atan(a_/b_*tan(gamma_));
    return calc_(t2);
}
Vec2 Ellipse2::calc_(Real t) const
{
    return Vec2(a_*cos(t)*cos(gamma_)-b_*sin(t)*sin(gamma_),a_*cos(t)*sin(gamma_)+b_*sin(t)*cos(gamma_))+origin_;
}
Rectangle2 Ellipse2::GetBoundingBox() const
{
  if(gamma_==0.0 || gamma_==M_PI){
    return Rectangle2(Vec2(origin_.x-a_,origin_.y-b_),Vec2(origin_.x+a_/2,origin_.y+b_));
  }
  if(gamma_==M_PI*0.5 || gamma_==3*M_PI*0.5){
    return Rectangle2(Vec2(origin_.x-b_,origin_.y-a_),Vec2(origin_.x+b_/2,origin_.y+a_));
  }
  Real t=atan(-b_/a_*tan(gamma_));
  Vec2 vx=calc_(t);
  t=atan(b_/a_/tan(gamma_));
  Vec2 vy=calc_(t);
  return Rectangle2(Vec2(origin_.x-vx.x,origin_.y-vy.y),Vec2(origin_.x+vx.x,origin_.y+vy.y));
}

Rectangle2::Rectangle2():
topleft_(Vec2(-1.0,-1.0)),
bottomright_(Vec2(1.0,1.0))
{
}
Rectangle2::Rectangle2(Vec2 topleft, Vec2 bottomright):
topleft_(topleft),
bottomright_(bottomright)
{
}
Real Rectangle2::GetWidth() const
{
    return bottomright_.x-topleft_.x;
}
Real Rectangle2::GetHeight() const
{
    return bottomright_.y-topleft_.y;
}
Real Rectangle2::GetArea() const
{
    return GetWidth()*GetHeight();
}
Vec2 Rectangle2::operator[] (unsigned index) const
{
    switch(index){
    case 0:
       return topleft_;
    case 1:
       return bottomright_;
    }
    throw(OutOfRangeException());    
}
Vec2& Rectangle2::operator[] (unsigned index)
{
    switch(index){
    case 0:
       return topleft_;
    case 1:
       return bottomright_;
    default:
       throw(OutOfRangeException());
    }
}

Vec2 Rectangle2::GetStart() const
{
  return topleft_;
}
Vec2 Rectangle2::GetEnd() const
{
  return bottomright_;
}
void Rectangle2::SetStart(const Vec2& v)
{
  topleft_=v;
}
void Rectangle2::SetEnd(const Vec2& v)
{
  bottomright_=v;
}

Vec2 Hyperbola2::At(Real t,bool righthalf) const
{
    Vec2 v(a_*cosh(t),b_*sinh(t));
    if(righthalf)
        return Rotate(v,gamma_)+origin_;
    else
        return origin_-Rotate(v,gamma_);

}

/************************************************************
*  Polygon2                                                  *
*************************************************************/
Polygon2::Polygon2():std::vector<Vec2>(){}
Polygon2::Polygon2(const Polygon2& p):std::vector<Vec2>(p){}
Polygon2::Polygon2(const std::vector<Vec2>& v):std::vector<Vec2>(v){}
void Polygon2::SetNode(unsigned int i,const Vec2& v)
{
    if(i<size()){
        operator[](i)=v;
    }
}
Polygon2 Polygon2::operator+(const Vec2& v) const
{
    Polygon2 result;
    for(Polygon2::const_iterator it=begin();it!=end();++it){
        result.AddNode((*it)+v);
    }
    return result;
}
Polygon2 Polygon2::operator*(Real d) const
{
    Polygon2 result;
    for(Polygon2::const_iterator it=begin();it!=end();++it){
        result.AddNode((*it)*d);
    }
    return result;
}
Real Polygon2::GetArea() const
{
  unsigned int nodecount=GetNodeCount();
  Real area=0;
  for(unsigned int i=0;i<nodecount;++i){
    unsigned int ip1=(i+1)%nodecount;
    Vec2 n1=GetNode(i);
    Vec2 n2=GetNode(ip1);
    area+=n1.x*n2.y-n2.x*n1.y;
  }
  return area/2.0;
}
Vec2 Polygon2::GetCentroid() const
{
  unsigned int nodecount=GetNodeCount();
  Real cx=0;
  Real cy=0;
  for(unsigned int i=0;i<nodecount;++i){
    unsigned int ip1=(i+1)%nodecount;
    Vec2 n1=GetNode(i);
    Vec2 n2=GetNode(ip1);
    cx+=(n1.x+n2.x)*(n1.x*n2.y-n2.x*n1.y);
    cy+=(n1.y+n2.y)*(n1.x*n2.y-n2.x*n1.y);
  }
  Real area=GetArea();
  return Vec2(cx/(6.0*area),cy/(6.0*area));
}
void Polygon2::Expand(Real val)
{
  Polygon2 result;
  unsigned int nodecount=GetNodeCount();
  for(unsigned int i=0;i<nodecount;++i){
    unsigned int ip1=(i+1)%nodecount;
    unsigned int im1=(i-1+nodecount)%nodecount;
    Vec2 node=GetNode(i);
    Vec2 nodep1=GetNode(ip1);
    Vec2 nodem1=GetNode(im1);
    Vec2 n=
      ((nodep1-node)/Length(nodep1-node)+(nodem1-node)/Length(nodem1-node))*0.5;
    n=n/Length(n)*val;
    if(!IsInPolygon(*this,node+n)){
      result.AddNode(node+n);
    }else{
      result.AddNode(node-n);
    }
  }
  for(unsigned int i=0;i<nodecount;++i){
    operator[](i)=result.GetNode(i);
  }
}
Rectangle2 Polygon2::GetBoundingBox() const
{
  if(GetNodeCount()==0)
  {
    return Rectangle2(Vec2(0.0,0.0),Vec2(0.0,0.0));
  }
  Real x_min=operator[](0).x;
  Real x_max=x_min;
  Real y_min=operator[](0).y;
  Real y_max=y_min;
  unsigned int nodecount=GetNodeCount();
  for(unsigned int i=0;i<nodecount;++i){
    Vec2 n=operator[](i);
    x_min=std::min<Real>(n.x,x_min);
    x_max=std::max<Real>(n.x,x_max);
    y_min=std::min<Real>(n.y,y_min);
    y_max=std::max<Real>(n.y,y_max);
  }
  return Rectangle2(Vec2(x_min,y_min),Vec2(x_max,y_max));
}

Polygon2::circular_iterator Polygon2::FindSegment(std::vector<Vec2>::const_iterator start_it, std::vector<Vec2>::const_iterator end_it)
{
  for(std::vector<Vec2>::iterator first=begin();first!=end();++first){
    if((*first)==(*start_it)){
      circular_iterator polygon_it(begin(),end(),first);
      bool flag=false;
      for(std::vector<Vec2>::const_iterator segment_it=start_it;segment_it!=end_it;++segment_it,++polygon_it){
        if((*segment_it)!=(*polygon_it)){
          flag=true;
          break;
        }
      }
      if(!flag){
        return circular_iterator(begin(),end(),first);
      }
    }
  }
  return cend();
}
void Polygon2::Erase(circular_iterator first,circular_iterator last)
{
  std::vector<Vec2>::iterator vfirst=first,vlast=last;
  if(first!=cbegin()){
    for(circular_iterator it=first;it!=last;++it){
      if(it==cbegin()){ //special case with wrap around
        std::vector<Vec2>::iterator vbegin=begin(),vend=end();
        std::vector<Vec2>::erase(vfirst,vend);
        std::vector<Vec2>::erase(vbegin,vlast);
      }
    }
  }
  std::vector<Vec2>::erase(vfirst,vlast);
}

Circle2::Circle2():
  center_(Vec2(0.0,0.0)),
  radius_(1.0)
{
}
Circle2::Circle2(const Circle2& c):
  center_(c.center_),
  radius_(c.radius_)
{
}
Circle2::Circle2(const Vec2& center, Real radius):
  center_(center),
  radius_(radius)
{
}
void Circle2::SetCenter(const Vec2& center)
{
  center_=center;
}
void Circle2::SetRadius(Real r)
{
  radius_=r;
}
Vec2 Circle2::GetCenter() const
{
  return center_;
}
Real Circle2::GetRadius() const
{
  return radius_;
}
Real Circle2::GetCircumference() const
{
  return 2*M_PI*radius_;
}
Real Circle2::GetArea() const
{
  return M_PI*radius_*radius_;
}


} // ns
