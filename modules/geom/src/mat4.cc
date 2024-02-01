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
#include "vec.hh"
#include "exc.hh"

#include "vecmat4_op.hh"

namespace geom {

Mat4::Mat4()
{
  this->set(1.0,0.0,0.0,0.0,
      0.0,1.0,0.0,0.0,
      0.0,0.0,1.0,0.0,
      0.0,0.0,0.0,1.0);
}

Mat4::Mat4(Real i00, Real i01, Real i02, Real i03, Real i10, Real i11, Real i12, Real i13,
           Real i20, Real i21, Real i22, Real i23, Real i30, Real i31, Real i32, Real i33)
{
  this->set(i00,i01,i02,i03,
      i10,i11,i12,i13,
      i20,i21,i22,i23,
      i30,i31,i32,i33);
}

Mat4::Mat4(const Mat4& m)
{
  this->set(m(0,0),m(0,1),m(0,2),m(0,3),
            m(1,0),m(1,1),m(1,2),m(1,3),
            m(2,0),m(2,1),m(2,2),m(2,3),
            m(3,0),m(3,1),m(3,2),m(3,3));
}

Mat4::Mat4(const Mat2& m)
{
  this->set(m(0,0),m(0,1),Real(0.0),Real(0.0),
            m(1,0),m(1,1),Real(0.0), Real(0.0),
            Real(0.0),Real(0.0),Real(1.0),Real(0.0),
            Real(0.0),Real(0.0),Real(0.0),Real(1.0));
}

Mat4::Mat4(const Mat3& m)
{
  this->set(m(0,0),m(0,1),m(0,2),Real(0.0),
            m(1,0),m(1,1),m(1,2),Real(0.0),
            m(2,0),m(2,1),m(2,2),Real(0.0),
            Real(0.0),Real(0.0),Real(0.0),Real(1.0));
}

Mat4::Mat4(const float arr[16])
{
  this->set(arr[0],arr[1],arr[2],arr[3],
      arr[4],arr[5],arr[6],arr[7],
      arr[8],arr[9],arr[10],arr[11],
      arr[12],arr[13],arr[14],arr[15]);
}

Mat4::Mat4(const double arr[16])
{
  this->set(arr[0],arr[1],arr[2],arr[3],
      arr[4],arr[5],arr[6],arr[7],
      arr[8],arr[9],arr[10],arr[11],
      arr[12],arr[13],arr[14],arr[15]);
}

Mat4& Mat4::operator=(const Mat4& m)
{
  if(&m!=this) {
    this->set(m(0,0),m(0,1),m(0,2),m(0,3),
        m(1,0),m(1,1),m(1,2),m(1,3),
        m(2,0),m(2,1),m(2,2),m(2,3),
        m(3,0),m(3,1),m(3,2),m(3,3));
  }
  return *this;
}

Mat4 Mat4::Identity() 
{
  static Mat4 i(1.0,0.0,0.0,0.0,
                0.0,1.0,0.0,0.0,
                0.0,0.0,1.0,0.0,
                0.0,0.0,0.0,1.0);
  return i;
}

bool Mat4::operator==(const Mat4& rhs) const
{
  return data_[0][0] ==  rhs.data_[0][0] &&
    data_[1][0] == rhs.data_[1][0] &&
    data_[2][0] == rhs.data_[2][0] &&
    data_[3][0] == rhs.data_[3][0] &&
    data_[0][1] == rhs.data_[0][1] &&
    data_[1][1] ==  rhs.data_[1][1] &&
    data_[2][1] ==  rhs.data_[2][1] &&
    data_[3][1] ==  rhs.data_[3][1] &&
    data_[0][2] ==  rhs.data_[0][2] &&
    data_[1][2] ==  rhs.data_[1][2] &&
    data_[2][2] ==  rhs.data_[2][2] &&
    data_[3][2] ==  rhs.data_[3][2] &&
    data_[0][3] ==  rhs.data_[0][3] &&
    data_[1][3] ==  rhs.data_[1][3] &&
    data_[2][3] ==  rhs.data_[2][3] &&
    data_[3][3] ==  rhs.data_[3][3];
}



Mat4& Mat4::operator+=(const Mat4& rhs)
{
  data_[0][0]+=rhs(0,0);
  data_[0][1]+=rhs(0,1);
  data_[0][2]+=rhs(0,2);
  data_[0][3]+=rhs(0,3);
  data_[1][0]+=rhs(1,0);
  data_[1][1]+=rhs(1,1);
  data_[1][2]+=rhs(1,2);
  data_[1][3]+=rhs(1,3);
  data_[2][0]+=rhs(2,0);
  data_[2][1]+=rhs(2,1);
  data_[2][2]+=rhs(2,2);
  data_[2][3]+=rhs(2,3);
  data_[3][0]+=rhs(3,0);
  data_[3][1]+=rhs(3,1);
  data_[3][2]+=rhs(3,2);
  data_[3][3]+=rhs(3,3);
  return *this;
}
Mat4& Mat4::operator-=(const Mat4& rhs)
{
  data_[0][0]-=rhs(0,0);
  data_[0][1]-=rhs(0,1);
  data_[0][2]-=rhs(0,2);
  data_[0][3]-=rhs(0,3);
  data_[1][0]-=rhs(1,0);
  data_[1][1]-=rhs(1,1);
  data_[1][2]-=rhs(1,2);
  data_[1][3]-=rhs(1,3);
  data_[2][0]-=rhs(2,0);
  data_[2][1]-=rhs(2,1);
  data_[2][2]-=rhs(2,2);
  data_[2][3]-=rhs(2,3);
  data_[3][0]-=rhs(3,0);
  data_[3][1]-=rhs(3,1);
  data_[3][2]-=rhs(3,2);
  data_[3][3]-=rhs(3,3);
  return *this;
}
Mat4& Mat4::operator*=(const Real d)
{
  data_[0][0]*=d;
  data_[0][1]*=d;
  data_[0][2]*=d;
  data_[0][3]*=d;
  data_[1][0]*=d;
  data_[1][1]*=d;
  data_[1][2]*=d;
  data_[1][3]*=d;
  data_[2][0]*=d;
  data_[2][1]*=d;
  data_[2][2]*=d;
  data_[2][3]*=d;
  data_[3][0]*=d;
  data_[3][1]*=d;
  data_[3][2]*=d;
  data_[3][3]*=d;
  return *this;
}
Mat4& Mat4::operator/=(const Real d)
{
  data_[0][0]/=d;
  data_[0][1]/=d;
  data_[0][2]/=d;
  data_[0][3]/=d;
  data_[1][0]/=d;
  data_[1][1]/=d;
  data_[1][2]/=d;
  data_[1][3]/=d;
  data_[2][0]/=d;
  data_[2][1]/=d;
  data_[2][2]/=d;
  data_[2][3]/=d;
  data_[3][0]/=d;
  data_[3][1]/=d;
  data_[3][2]/=d;
  data_[3][3]/=d;
  return *this;
}

Mat4& Mat4::operator*=(const Mat4& m)
{
  (*this)=operator*(*this,m);
  return *this;
}

Mat3 Mat4::ExtractRotation() const
{
  Mat3 nrvo;
  for(int i=0;i<3;++i) {
    for(int k=0;k<3;++k) {
      nrvo(i,k)=data_[i][k];
    }
  }
  return nrvo;
}

void Mat4::PasteRotation(const Mat3& m)
{
  for(int i=0;i<3;++i) {
    for(int k=0;k<3;++k) {
      data_[i][k]=m(i,k);
    }
  }
}

Vec3 Mat4::ExtractTranslation() const
{
  Vec3 nrvo;
  nrvo[0]=data_[0][3];
  nrvo[1]=data_[1][3];
  nrvo[2]=data_[2][3];
  return nrvo;
}

void Mat4::PasteTranslation(const Vec3& v)
{
  data_[0][3]=v[0];
  data_[1][3]=v[1];
  data_[2][3]=v[2];
}

////////////////////////////////////////////
// private member functions
void Mat4::set(double i00, double i01, double i02, double i03, double i10, double i11,
double i12, double i13, double i20, double i21, double i22, double i23, double i30, double
i31, double i32, double i33)
{
  data_[0][0]=i00; data_[0][1]=i01; data_[0][2]=i02; data_[0][3]=i03;
  data_[1][0]=i10; data_[1][1]=i11; data_[1][2]=i12; data_[1][3]=i13;
  data_[2][0]=i20; data_[2][1]=i21; data_[2][2]=i22; data_[2][3]=i23;
  data_[3][0]=i30; data_[3][1]=i31; data_[3][2]=i32; data_[3][3]=i33;
}

void Mat4::set(float i00, float i01, float i02, float i03, float i10, float i11, float
i12, float i13, float i20, float i21, float i22, float i23, float i30, float i31, float
i32, float i33)
{
  data_[0][0]=i00; data_[0][1]=i01; data_[0][2]=i02; data_[0][3]=i03;
  data_[1][0]=i10; data_[1][1]=i11; data_[1][2]=i12; data_[1][3]=i13;
  data_[2][0]=i20; data_[2][1]=i21; data_[2][2]=i22; data_[2][3]=i23;
  data_[3][0]=i30; data_[3][1]=i31; data_[3][2]=i32; data_[3][3]=i33;
}


std::ostream& operator<<(std::ostream& os, const Mat4& m)
{
  os << "{{" << m(0,0) << "," << m(0,1) << "," << m(0,2) << "," << m(0,3) << "},";
  os << "{" << m(1,0) << "," << m(1,1) << "," << m(1,2) << "," << m(1,3) << "},";
  os << "{" << m(2,0) << "," << m(2,1) << "," << m(2,2) << "," << m(2,3) << "},";
  os << "{" << m(3,0) << "," << m(3,1) << "," << m(3,2) << "," << m(3,3) << "}}";
  return os;
}

} // ns
