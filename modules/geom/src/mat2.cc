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
#include "mat.hh"
#include "exc.hh"

#include "vecmat2_op.hh"

namespace geom {

Mat2::Mat2()
{
  this->set(1.0,0.0, 0.0,1.0);
}

Mat2::Mat2(Real i00, Real i01, Real i10, Real i11)
{
  this->set(i00,i01,i10,i11);
}

Mat2::Mat2(const Mat2& m)
{
  this->set(m(0,0),m(0,1),m(1,0),m(1,1));
}

Mat2::Mat2(const Real arr[4])
{
  this->set(arr[0],arr[1],arr[2],arr[3]);
}

Mat2& Mat2::operator=(const Mat2& m)
{
  if(&m!=this) {
    this->set(m(0,0),m(0,1),m(1,0),m(1,1));
  }
  return *this;
}

Mat2 Mat2::Identity() 
{
  static Mat2 i(1.0,0.0,0.0,1.0);
  return i;
}

bool Mat2::operator==(const Mat2& rhs) const
{
  return data_[0][0] == rhs.data_[0][0] &&
    data_[1][0] == rhs.data_[1][0] &&
    data_[0][1] == rhs.data_[0][1] &&
    data_[1][1] == rhs.data_[1][1];
}

Mat2& Mat2::operator+=(const Mat2& rhs)
{
  data_[0][0]+=rhs(0,0);
  data_[0][1]+=rhs(0,1);
  data_[1][0]+=rhs(1,0);
  data_[1][1]+=rhs(1,1);
  return *this;
}

Mat2& Mat2::operator-=(const Mat2& rhs)
{
  data_[0][0]-=rhs(0,0);
  data_[0][1]-=rhs(0,1);
  data_[1][0]-=rhs(1,0);
  data_[1][1]-=rhs(1,1);
  return *this;
}

Mat2& Mat2::operator*=(const Real d)
{
  data_[0][0]*=d;
  data_[0][1]*=d;
  data_[1][0]*=d;
  data_[1][1]*=d;
  return *this;
}

Mat2& Mat2::operator/=(const Real d)
{
  data_[0][0]/=d;
  data_[0][1]/=d;
  data_[1][0]/=d;
  data_[1][1]/=d;
  return *this;
}

Mat2& Mat2::operator*=(const Mat2& m)
{
  (*this)=operator*(*this,m);
  return *this;
}

////////////////////////////////////////////
// private member functions
void Mat2::set(Real i00, Real i01, Real i10, Real i11)
{
  data_[0][0]=i00; data_[0][1]=i01;
  data_[1][0]=i10; data_[1][1]=i11;
}


std::ostream& operator<<(std::ostream& os, const Mat2& m)
{
  os << "{{" << m(0,0) << "," << m(0,1) << "},{" << m(1,0) << "," << m(1,1) << "}}";
  return os;
}

} // ns
