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
  value holder index

  Author: Ansgar Philippsen
*/
#include <iostream>

#include "index.hh"

namespace ost { namespace img { namespace image_state {

// Index


bool Index::equal(const Index& i) const {
  return u==i.u && v==i.v && w==i.w;
}


// IndexIterator
IndexIterator::IndexIterator(int w, int h, int d):
  w_(w),
  h_(h),
  d_(d),
  indx_(0,0,0),
  at_end_(false)
{}

IndexIterator::IndexIterator(const Size& s):
  w_(s.GetWidth()),
  h_(s.GetHeight()),
  d_(s.GetDepth()),
  indx_(0,0,0),
  at_end_(false)
{}

bool IndexIterator::AtStart() const 
{
  return (indx_.u==0 && indx_.v==0 && indx_.w==0);
}

bool IndexIterator::AtEnd() const 
{
  return at_end_;
}

IndexIterator& IndexIterator::ToStart()
{
  indx_.u=0;
  indx_.v=0;
  indx_.w=0;
  at_end_=false;
  return *this;
}

IndexIterator& IndexIterator::ToEnd()
{
  indx_.u=w_;
  indx_.v=h_;
  indx_.w=d_;
  at_end_=true;
  return *this;
}

// prefix
IndexIterator& IndexIterator::operator++() 
{
  increment();
  return *this;
}

// postfix
IndexIterator IndexIterator::operator++(int) 
{
  IndexIterator tmp(*this);
  increment();
  return tmp;
}

// prefix
IndexIterator& IndexIterator::operator--() 
{
  decrement();
  return *this;
}

// postfix
IndexIterator IndexIterator::operator--(int) 
{
  IndexIterator tmp(*this);
  decrement();
  return tmp;
}

IndexIterator::operator Index() const 
{
  return indx_;
}

void IndexIterator::increment() 
{
  if(++indx_.u>=w_) {
    indx_.u=0;
    if(++indx_.v>=h_) {
      indx_.v=0;
      if(++indx_.w>=d_) {
	indx_.u=w_;
	indx_.v=h_;
	indx_.w=d_;
	at_end_=true;
      }
    }
  }
}

void IndexIterator::decrement() 
{
  if(at_end_) {
    at_end_=false;
    indx_.u=w_;
    indx_.v=h_-1;
    indx_.w=d_-1;
  }
  if(indx_.u>0) {
    --indx_.u;
  } else {
    indx_.u=w_-1;
    if(indx_.v>0) {
      --indx_.v;
    } else {
      indx_.v=h_-1;
      if(indx_.w>0) {
	--indx_.w;
      } else {
	indx_.u=0;
	indx_.v=0;
	indx_.w=0;
      }
    }
  }
}

// global operators

bool operator==(const Index& i1, const Index& i2) {
  return i1.equal(i2);
}

bool operator!=(const Index& i1, const Index& i2) {
  return !i1.equal(i2);
}

std::ostream& operator<<(std::ostream& out, const Index& i)
{
  out << "(" << i.u << "," << i.v << "," << i.w << ")"; 
  return out;
}

}}} // namespaces
