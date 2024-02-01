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
  Authors: Ansgar Philippsen, Giani Signorell
*/

#include "line_iterator.hh"

#include <ost/img/extent_iterator.hh>
#include <ost/message.hh>

namespace ost {
namespace img {
namespace alg {


namespace 
{
    bool setx( unsigned int flag) {return (flag & Axis::X)!=0;}
    bool sety( unsigned int flag) {return (flag & Axis::Y)!=0;}
    bool setz( unsigned int flag) {return (flag & Axis::Z)!=0;}
}

LineIterator::LineIterator(const Extent& ext, unsigned int flag):
  flag_(flag),
  ext_(ext),
  start_(ext_.GetStart()),
  end_(ext_.GetEnd()),
  indx_(start_)
{

   if (!setx(flag_) && sety(flag_))  {
    current_dir_=Axis::Y;
  } else if (!setx(flag_) && !sety(flag_)) {
    current_dir_=Axis::Z; 
  } else {
    current_dir_=Axis::X;
  }

  if (current_dir_==Axis::Y && ext_.GetHeight()==1 && ext_.GetDepth()==1) {
    indx_=end_;
  }

  if (current_dir_==Axis::Z && ext_.GetDepth()==1) {
    indx_=end_;
  }

}

ExtentIterator LineIterator::GetLine()
{
  if( current_dir_ == Axis::X ){
    return ExtentIterator(Extent(Point(start_[0],indx_[1],indx_[2]),
				 Size(ext_.GetWidth(),0,0)));
  } else if( current_dir_ == Axis::Y ) {
    return ExtentIterator(Extent(Point(indx_[0],start_[1],indx_[2]),
				 Size(0,ext_.GetHeight(),0)));
  } else {
    return ExtentIterator(Extent(Point(indx_[0],indx_[1],start_[2]),
				 Size(0,0,ext_.GetDepth())));
  }
}


int LineIterator::GetSize() const
{
  if( current_dir_ == Axis::X ) {
    return ext_.GetWidth();
  } else if ( current_dir_ == Axis::Y) {
    return ext_.GetHeight();
  } else{
    return ext_.GetDepth();
  }
}

void LineIterator::Inc()
{
  if(current_dir_==Axis::X) {
    indx_[0]=start_[0];
    if(++indx_[1]>end_[1]) {
      indx_[1]=start_[1];
      if(++indx_[2]>end_[2]) {
	if(ext_.GetHeight()==1) {
	  if(ext_.GetDepth()==1) {
	    indx_=end_;
	  } 
        } else if (sety(flag_) && ext_.GetHeight()!=1) {
	  indx_=start_;
	  current_dir_=Axis::Y;
	} else if (setz(flag_) && ext_.GetDepth()!=1) {
	  indx_=start_;
	  current_dir_=Axis::Z;
        } else {
          indx_=end_;
        }  
      }
    }
  } else if(current_dir_==Axis::Y) {
    indx_[1]=start_[1];
    if(++indx_[0]>end_[0]) {
      indx_[0]=start_[0];
      if(++indx_[2]>end_[2]) {
	if(ext_.GetDepth()==1) {
	  indx_=end_;
	} else if (setz(flag_) && ext_.GetDepth()!=1) {
	  indx_=start_;
	  current_dir_=Axis::Z;
	} else {
          indx_=end_;
        }          
      }
    }
  } else {
    indx_[2]=start_[2];    
    if(++indx_[0]>end_[0]) {
      indx_[0]=start_[0];
      if(++indx_[1]>end_[1]) {
	indx_=end_;
      } 
    }
  }
}


bool LineIterator::AtEnd()
{
  return indx_==end_;
}

LineIterator& LineIterator::operator++() 
{
  Inc();
  return *this;
}

LineIterator LineIterator::operator++(int) 
{
  LineIterator tmp(*this);
  this->operator++();
  return tmp;
}

ExtentIterator LineIterator::GetNextLine() 
{
  ExtentIterator ti=GetLine();
  Inc();
  return ti;
}

}}} // ns
