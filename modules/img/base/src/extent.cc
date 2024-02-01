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
  Extent implementation

  Author: Ansgar Philippsen
*/

#include <cstdlib>
#include <iostream>

#include <ost/img/module_config.hh>
#include "extent.hh"
#include "vecmat.hh"

namespace ost { namespace img {

Extent::Extent():
  start_(),
  end_()
{}

Extent::Extent(const Extent &r):
  start_(r.start_),
  end_(r.end_)
{}

Extent::Extent(const Point& p1, const Point& p2):
  start_(std::min(p1[0],p2[0]),std::min(p1[1],p2[1]),std::min(p1[2],p2[2])),
  end_(std::max(p1[0],p2[0]),std::max(p1[1],p2[1]),std::max(p1[2],p2[2]))
{}

Extent::Extent(const Size& size):
  start_(),
  end_(size - Point(1,1,1))
{
}

Extent::Extent(const Size& size, const Point& cen):
  start_(-size.GetHalf() + cen),
  end_(size + start_ - Point(1,1,1))
{
}

Extent::Extent(const Point& start, const Size& size):
  start_(start),
  end_(size + start_ - Point(1,1,1))
{
}

const Point& Extent::GetStart() const {return start_;}

void Extent::SetStart(const Point& o) 
{
  set(o,end_);
}

void Extent::SetEnd(const Point& e) 
{
  set(start_,e);
}

const Point& Extent::GetEnd() const {return end_;}

bool Extent::Contains(const Point& p) const
{
  return (p[0]>=start_[0] && p[0]<=end_[0] &&
          p[1]>=start_[1] && p[1]<=end_[1] &&
          p[2]>=start_[2] && p[2]<=end_[2]);
}

bool Extent::Contains(const Extent& e) const
{
  return Contains(e.GetStart()) && Contains(e.GetEnd());
}

Point Extent::GetCenter() const 
{
  return Size(start_,end_).GetHalf()+start_;
}

Size Extent::GetSize() const {return Size(start_,end_);}

int Extent::GetWidth() const {return end_[0]-start_[0]+1;}

int Extent::GetHeight() const {return end_[1]-start_[1]+1;}

int Extent::GetDepth() const {return end_[2]-start_[2]+1;}

int Extent::GetVolume() const {return Size(start_,end_).GetVolume();}

int Extent::GetDim() const {return Size(start_,end_).GetDim();}


Point Extent::WrapAround(const Point& p)
{
  Point r(p-start_);
  Size size=GetSize();
  for(int i=0;i<3;i++) {
    r[i] = start_[i] + (r[i]<0 ? size[i]+std::div(r[i],size[i]).rem : std::div(r[i],size[i]).rem);
  }
  return r;
}


Extent Extent::Mirror(int planes)
{
  Point new_start(planes & Plane::YZ ? -end_[0] : start_[0],
                  planes & Plane::XZ ? -end_[1] : start_[1],
                  planes & Plane::XY ? -end_[2] : start_[2]);
  return Extent(new_start,GetSize());
}

unsigned int Extent::Point2Offset(const Point& p)
{
  if(this->Contains(p)) {
    Point d(p-start_);
    return d[0]+GetWidth()*(d[1]+d[2]*GetHeight());
  } else {
    return 0;
  }
}

void Extent::Shift(const Point& p)
{
  set(start_+p,end_+p);
}

// private methods

void Extent::set(const Point& p1, const Point& p2) 
{
  // order them automatically
  for(int i=0;i<3;i++) {
    if(p1[i]<=p2[i]) {
      start_[i]=p1[i];
      end_[i]=p2[i];
    } else {
      start_[i]=p2[i];
      end_[i]=p1[i];
    }
  }
}

bool Extent::equal(const Extent& b) const 
{
  return (start_==b.start_ && end_==b.end_);
}

// 'global' functions


Extent Overlap(const Extent& e1, const Extent& e2)
{
  Point a,b;

  for(int k=0;k<3;++k) {
    a[k] = std::max(e1.GetStart()[k],e2.GetStart()[k]);
    b[k] = std::min(e1.GetEnd()[k],e2.GetEnd()[k]);
    if(b[k]<a[k]) throw InvalidExtentException("non-overlapping extents in call to Overlap");
  }
  return Extent(a,b);
}

bool HasOverlap(const Extent& e1, const Extent& e2)
{
  Point a,b;

  for(int k=0;k<3;++k) {
    a[k] = std::max(e1.GetStart()[k],e2.GetStart()[k]);
    b[k] = std::min(e1.GetEnd()[k],e2.GetEnd()[k]);
    if(b[k]<a[k]) return false;
  }
  return true;
}



std::ostream& operator<<(std::ostream& os, const Extent &b) 
{
  os << "[" << b.GetStart() << "," << b.GetEnd() << "]";
  return os;
}


}}// namespace
