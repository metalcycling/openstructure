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
  extent iterator

  Author: Ansgar Philippsen
*/

#include "extent_iterator.hh"
#include "extent.hh"

namespace ost { namespace img {

ExtentIterator::ExtentIterator(const Extent& ext, unsigned int domain):
  start_(ext.GetStart()),
  end_(ext.GetEnd()+Point(1,1,1))
{
  if(domain==HALF_FREQUENCY){
    startpos_=Point(0,0,0);
  }else{
    startpos_=start_;
  }
  current_=startpos_;
    
}

ExtentIterator::ExtentIterator(const ExtentIterator& i):
  start_(i.start_),
  end_(i.end_),
  current_(i.current_),
  startpos_(i.startpos_)
{}

int ExtentIterator::operator[](int i) const
{
  return current_[i];
}

ExtentIterator& ExtentIterator::operator=(const ExtentIterator& i)
{
  start_=i.start_;
  startpos_=i.startpos_;
  end_=i.end_;
  current_=i.current_;
  return *this;
}

ExtentIterator& ExtentIterator::operator++() 
{
  Inc();
  return *this;
}

ExtentIterator ExtentIterator::operator++(int) 
{
  ExtentIterator tmp(*this);
  Inc();
  return tmp;
}

ExtentIterator & ExtentIterator::operator--() 
{
  Dec();
  return *this;
}

ExtentIterator ExtentIterator::operator--(int) 
{
  ExtentIterator tmp(*this);
  Dec();
  return tmp;
}


void ExtentIterator::Inc() 
{
  if(++current_[0]>=end_[0]) {
    current_[0]=start_[0];
    if(++current_[1]>=end_[1]) {
      current_[1]=start_[1];
      if(++current_[2]>=end_[2]) {
	current_=end_;
      }
    } 
  }
}

void ExtentIterator::Dec() 
{
  if(current_==startpos_)return;
  if(AtEnd()) {
    --current_[0];
    --current_[1];
    --current_[2];
  } else {
    if(--current_[0]<start_[0]) {
      current_[0]=end_[0]-1;
      if(--current_[1]<start_[1]) {
	current_[1]=end_[1]-1;
	if(--current_[2]<start_[2]) {
	  current_=start_;
	}
      } 
    }
  }
}

ExtentIterator& ExtentIterator::ToStart() 
{
  current_=startpos_;
  return *this;
}

ExtentIterator& ExtentIterator::ToEnd()
{
  current_=end_;
  return *this;
}

Extent ExtentIterator::GetExtent() const
{
  return Extent(start_,end_);
}

const Point& ExtentIterator::GetStart() const 
{
  return startpos_;
}

const Point &ExtentIterator::GetEnd() const
{
  return end_;
}

bool ExtentIterator::AtEnd() const
{
  return current_==end_;
}

bool ExtentIterator::AtStart() const
{
  return current_==startpos_;
}

ExtentIterator::operator Point() const
{
  return current_;
}


}} // namespace
