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
  size

  Author: Ansgar Philippsen
*/

#include <ost/message.hh>
#include <iostream>
#include <sstream>

#include "size.hh"
#include "point.hh"

namespace ost { namespace img {

Size::Size():
  w_(1),
  h_(1),
  d_(1)
{}

Size::Size(const Size &s):
  w_(s.w_),
  h_(s.h_),
  d_(s.d_)
{}

Size::Size(unsigned int w):
  w_(w<1 ? 1 : w),
  h_(1),
  d_(1)
{}

Size::Size(unsigned int w, unsigned int h):
  w_(w<1 ? 1 : w),
  h_(h<1 ? 1 : h),
  d_(1)
{}

Size::Size(unsigned int w, unsigned int h, unsigned int d):
  w_(w<1 ? 1 : w),
  h_(h<1 ? 1 : h),
  d_(d<1 ? 1 : d)
{}

Size::Size(const Point& p1, const Point& p2):
  w_(std::abs(p2[0]-p1[0])+1<1 ? 1 : std::abs(p2[0]-p1[0])+1),
  h_(std::abs(p2[1]-p1[1])+1<1 ? 1 : std::abs(p2[1]-p1[1])+1),
  d_(std::abs(p2[2]-p1[2])+1<1 ? 1 : std::abs(p2[2]-p1[2])+1)
{}

Point Size::GetHalf() const 
{
  return Point((w_&0x1) ? (w_-1)/2 : w_/2-1,
	       (h_&0x1) ? (h_-1)/2 : h_/2-1,
	       (d_&0x1) ? (d_-1)/2 : d_/2-1);
}

unsigned int Size::GetWidth() const {return w_;}

unsigned int Size::GetHeight() const {return h_;}

unsigned int Size::GetDepth() const {return d_;}

unsigned int Size::GetSlab() const {return w_*h_;}

unsigned int Size::GetVol() const {return GetVolume();}

unsigned int Size::GetVolume() const {return w_*h_*d_;}


unsigned int Size::GetDim() const 
{
  return d_>1 ? 3 : (h_>1 ? 2 : 1);
}


Size& Size::operator=(const Size& s)
{
  w_=s.w_;
  h_=s.h_;
  d_=s.d_;
  return *this;
}

unsigned int& Size::operator[](unsigned int index) 
{
  if(index==0) return w_;
  else if(index==1) return h_;
  else if(index==2) return d_;
  else throw geom::OutOfRangeException("Size index out of range");
}

unsigned int Size::operator[](unsigned int index) const 
{
  if(index==0) return w_;
  else if(index==1) return h_;
  else if(index==2) return d_;
  else throw geom::OutOfRangeException("Size index out of range");
}

bool Size::operator==(const Size& s) const 
{
  return equal(s);
}

bool Size::operator!=(const Size& s) const 
{
  return !equal(s);
}

Size& Size::operator+=(const Size& s)
{
  w_+=s.w_;
  h_+=s.h_;
  d_+=s.d_;
  return *this;
}

// clamped to 1
Size& Size::operator-=(const Size& s)
{
  w_ = (w_>s.w_) ? w_-s.w_ : 1;
  h_ = (h_>s.h_) ? h_-s.h_ : 1;
  d_ = (d_>s.d_) ? d_-s.d_ : 1;
  return *this;
}

bool Size::equal(const Size& s) const
{
  return (w_==s.w_ &&
	  h_==s.h_ &&
	  d_==s.d_);
}

// global functions

Size operator+(const Size& s1, const Size& s2)
{
  Size tmp(s1);
  return tmp+=s2;
}

Size operator-(const Size& s1, const Size& s2)
{
  Size tmp(s1);
  return tmp-=s2;
}

std::ostream& operator<<(std::ostream& os, const img::Size &s) 
{
  os << "(" << s.GetWidth() << "," << s.GetHeight() << "," << s.GetDepth() << ")";
  return os;
}


}} // namespace
