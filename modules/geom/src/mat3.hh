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
#ifndef GEOM_MAT3_HH
#define GEOM_MAT3_HH

#include <cstddef> // for size_t
#include <ostream>
#include <cassert>
#include <stdexcept>

#include <boost/operators.hpp>

#include <ost/geom/module_config.hh>
#include <ost/geom/mat2.hh>
namespace geom {

class Vec3;

class DLLEXPORT_OST_GEOM Mat3:
    private boost::equality_comparable<Mat3>,
    private boost::additive1<Mat3>,
    private boost::multiplicative2<Mat3, Real>
{
public:  
  //! Default initialization, identity matrix  
  Mat3()
  {
    this->set(1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0);
  }
  //! In with 9 values in row-major order
  /*!
    row-major order means that the matrix

    | a b c |
    | d e f |
    | g h i |

    is initialized with (a,b,c, d,e,f, g,h,i)
  */
  Mat3(Real i00, Real i01, Real i02,
       Real i10, Real i11, Real i12,
       Real i20, Real i21, Real i22)
  {
    this->set(i00,i01,i02,i10,i11,i12,i20,i21,i22);
  }

  Mat3(const Mat3& m)
  {
    this->set(m(0,0),m(0,1),m(0,2),m(1,0),m(1,1),m(1,2),m(2,0),m(2,1),m(2,2));
  }

  Mat3(const Mat2& m)
  {
    this->set(m(0, 0), m(0, 1), Real(0.0),
              m(1, 0), m(1, 1), Real(0.0),
              Real(0.0), Real(0.0), Real(1.0));
  }
  
  explicit Mat3(const Real arr[9])
  {
    this->set(arr[0],arr[1],arr[2],arr[3],arr[4],arr[5],arr[6],arr[7],arr[8]);
  }

  explicit Mat3(Real x, Real y, Real z)
  {
    this->set(x, 0.0, 0.0, 0.0, y, 0.0, 0.0, 0.0, z);
  }
  //! element access
  Real& At(std::size_t r, std::size_t c)
  {
    if (r>2 || c>2) {
      throw std::out_of_range("indices must be smaller than 3");
    }
    return data_[r][c];
  }
  //! element access
  const Real& At(std::size_t r, std::size_t c) const
  {
    if (r>2 || c>2) {
      throw std::out_of_range("indices must be smaller than 3");
    }
    return data_[r][c];
  }
  //! element access
  Real& operator()(std::size_t r, std::size_t c)
  {
    assert(r<3 && c<3);
    return data_[r][c];
  }
  //! const element access
  const Real& operator()(std::size_t r, std::size_t c) const
  {
    assert(r<3 && c<3);
    return data_[r][c];
  }

  Mat3& operator=(const Mat3& m)
  {
    if(&m!=this) {
      this->set(m(0,0),m(0,1),m(0,2),m(1,0),m(1,1),m(1,2),m(2,0),m(2,1),m(2,2));
    }
    return *this;
  }

  static Mat3 Identity() 
  {
    static Mat3 i(1.0,0.0,0.0,
                  0.0,1.0,0.0,
                  0.0,0.0,1.0);
    return i;
  }

  bool operator==(const Mat3& rhs) const
  {
    return data_[0][0] ==  rhs.data_[0][0] &&
      data_[1][0] == rhs.data_[1][0] &&
      data_[2][0] == rhs.data_[2][0] &&
      data_[0][1] == rhs.data_[0][1] &&
      data_[1][1] ==  rhs.data_[1][1] &&
      data_[2][1] ==  rhs.data_[2][1] &&
      data_[0][2] ==  rhs.data_[0][2] &&
      data_[1][2] ==  rhs.data_[1][2] &&
      data_[2][2] ==  rhs.data_[2][2];
  }

  Mat3& operator+=(const Mat3& rhs)
  {
    data_[0][0]+=rhs(0,0);
    data_[0][1]+=rhs(0,1);
    data_[0][2]+=rhs(0,2);
    data_[1][0]+=rhs(1,0);
    data_[1][1]+=rhs(1,1);
    data_[1][2]+=rhs(1,2);
    data_[2][0]+=rhs(2,0);
    data_[2][1]+=rhs(2,1);
    data_[2][2]+=rhs(2,2);
    return *this;
  }
  Mat3& operator-=(const Mat3& rhs)
  {
    data_[0][0]-=rhs(0,0);
    data_[0][1]-=rhs(0,1);
    data_[0][2]-=rhs(0,2);
    data_[1][0]-=rhs(1,0);
    data_[1][1]-=rhs(1,1);
    data_[1][2]-=rhs(1,2);
    data_[2][0]-=rhs(2,0);
    data_[2][1]-=rhs(2,1);
    data_[2][2]-=rhs(2,2);
    return *this;
  }
  Mat3& operator*=(const Real d)
  {
    data_[0][0]*=d;
    data_[0][1]*=d;
    data_[0][2]*=d;
    data_[1][0]*=d;
    data_[1][1]*=d;
    data_[1][2]*=d;
    data_[2][0]*=d;
    data_[2][1]*=d;
    data_[2][2]*=d;
    return *this;
  }
  Mat3& operator/=(const Real d)
  {
    data_[0][0]/=d;
    data_[0][1]/=d;
    data_[0][2]/=d;
    data_[1][0]/=d;
    data_[1][1]/=d;
    data_[1][2]/=d;
    data_[2][0]/=d;
    data_[2][1]/=d;
    data_[2][2]/=d;
    return *this;
  }

  Mat3& operator*=(const Mat3& m)
  {
    (*this)=Mat3((*this)(0,0)*m(0,0)+(*this)(0,1)*m(1,0)+(*this)(0,2)*m(2,0),
                 (*this)(0,0)*m(0,1)+(*this)(0,1)*m(1,1)+(*this)(0,2)*m(2,1),
                 (*this)(0,0)*m(0,2)+(*this)(0,1)*m(1,2)+(*this)(0,2)*m(2,2),
                 (*this)(1,0)*m(0,0)+(*this)(1,1)*m(1,0)+(*this)(1,2)*m(2,0),
                 (*this)(1,0)*m(0,1)+(*this)(1,1)*m(1,1)+(*this)(1,2)*m(2,1),
                 (*this)(1,0)*m(0,2)+(*this)(1,1)*m(1,2)+(*this)(1,2)*m(2,2),                                   
                 (*this)(2,0)*m(0,0)+(*this)(2,1)*m(1,0)+(*this)(2,2)*m(2,0),
                 (*this)(2,0)*m(0,1)+(*this)(2,1)*m(1,1)+(*this)(2,2)*m(2,1),
                 (*this)(2,0)*m(0,2)+(*this)(2,1)*m(1,2)+(*this)(2,2)*m(2,2));
    return *this;
  }

  Real* Data() {return &data_[0][0];}
  const Real* Data() const {return &data_[0][0];}
  
  geom::Vec3 GetCol(int index) const;
  geom::Vec3 GetRow(int index) const;
private:
  Real data_[3][3];

  void set(Real i00, Real i01, Real i02,
           Real i10, Real i11, Real i12,
           Real i20, Real i21, Real i22)
  {
    data_[0][0]=i00; data_[0][1]=i01; data_[0][2]=i02;
    data_[1][0]=i10; data_[1][1]=i11; data_[1][2]=i12;
    data_[2][0]=i20; data_[2][1]=i21; data_[2][2]=i22;
  }
};
  
DLLEXPORT_OST_GEOM std::ostream& operator<<(std::ostream& o, const Mat3& m);

} // ns geom


#endif
