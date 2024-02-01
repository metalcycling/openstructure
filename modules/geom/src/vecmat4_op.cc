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

#include <ost/base.hh>
#include "vecmat4_op.hh"
#include "vecmat3_op.hh"

#include "exc.hh"


namespace geom {

Real Comp(const Mat4& m, unsigned int i_in, unsigned int j_in)
{
  return Minor(m,i_in,j_in)*( (i_in+j_in) & 0x1 ? -1.0 : 1.0);
}

Real Minor(const Mat4& m,unsigned int i_in, unsigned int j_in)
{
  if(i_in>3 || j_in>3){
    throw OutOfRangeException();
  }
  Mat3 r;
  unsigned int itmp=0;
  for (unsigned int i=0;i<4;i++){
    if(i==i_in) continue;
    unsigned int jtmp=0;
    for (unsigned int j=0;j<4;j++){
      if(j==j_in) continue;
      r(itmp,jtmp)=m(i,j);
      jtmp++;
    }
    itmp++;
  }
  
  return Det(r);
}

Real Det(const Mat4& m)
{
  Real d=0;
  for(int i=0;i<4;i++)
  {
    d+=m(0,i)*Comp(m,0,i);
  }
  return d;
}

Mat4 Transpose(const Mat4& m)
{
  Mat4 r;
  for (int i=0;i<4;i++){
    for (int j=0;j<4;j++){
      r(j,i)=m(i,j);
    }
  }
  return r;
}

Mat4 Invert(const Mat4& in)
{
  Mat4 m(in);
  int indxc[]={0,0,0,0};
  int indxr[]={0,0,0,0};
  int ipiv[]={0,0,0,0};
  int icol=0, irow=0;
  Real dum=0.0, pivinv=0.0;
  Real b[][1]={{0.0},{0.0},{0.0},{0.0}};
  
  for(int i=0;i<4;++i) {
    Real big=0.0;
    for(int j=0;j<4;++j) {
      if(ipiv[j]!=1) {
  for(int k=0;k<4;++k) {
    if(ipiv[k]==0) {
      if(std::fabs(m(j,k)) >= big) {
        big=std::fabs(m(j,k));
        irow=j;
        icol=k;
      }
    } else if(ipiv[k] >1) {
      throw GeomException("error in matrix inversion");
    }
  }
      }
    }
    ++(ipiv[icol]);
    if(irow!=icol) {
      for(int l=0;l<4;++l) std::swap(m(irow,l),m(icol,l));
      for(int l=0;l<1;++l) std::swap(b[irow][l], b[icol][l]);
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if(m(icol,icol)==0.0) {
      throw GeomException("error in matrix inversion");
    }
    pivinv=1.0/m(icol,icol);
    m(icol,icol)=1.0;
    for(int l=0;l<4;++l) m(icol,l)*=pivinv;
    for(int l=0;l<1;++l) b[icol][l]*=pivinv;

    for(int ll=0;ll<4;++ll) {
      if(ll!=icol) {
  dum=m(ll,icol);
  m(ll,icol)=0.0;
  for(int l=0;l<4;++l) m(ll,l) -=m(icol,l)*dum;
  for(int l=0;l<1;++l) b[ll][l] -=b[icol][l]*dum;
      }
    }
  }
  for(int l=3;l>=0;--l) {
    if(indxr[l]!=indxc[l]) {
      for(int k=0;k<4;++k) {
  std::swap(m(k,indxr[l]),m(k,indxc[l]));
      }
    }
  }

  return m;
}

Mat4 operator*(const Mat4& m1, const Mat4& m2)
{
  Mat4 nrvo(m1(0,0)*m2(0,0)+m1(0,1)*m2(1,0)+m1(0,2)*m2(2,0)+m1(0,3)*m2(3,0),
      m1(0,0)*m2(0,1)+m1(0,1)*m2(1,1)+m1(0,2)*m2(2,1)+m1(0,3)*m2(3,1),
      m1(0,0)*m2(0,2)+m1(0,1)*m2(1,2)+m1(0,2)*m2(2,2)+m1(0,3)*m2(3,2),
      m1(0,0)*m2(0,3)+m1(0,1)*m2(1,3)+m1(0,2)*m2(2,3)+m1(0,3)*m2(3,3),

      m1(1,0)*m2(0,0)+m1(1,1)*m2(1,0)+m1(1,2)*m2(2,0)+m1(1,3)*m2(3,0),
      m1(1,0)*m2(0,1)+m1(1,1)*m2(1,1)+m1(1,2)*m2(2,1)+m1(1,3)*m2(3,1),
      m1(1,0)*m2(0,2)+m1(1,1)*m2(1,2)+m1(1,2)*m2(2,2)+m1(1,3)*m2(3,2),
      m1(1,0)*m2(0,3)+m1(1,1)*m2(1,3)+m1(1,2)*m2(2,3)+m1(1,3)*m2(3,3),

      m1(2,0)*m2(0,0)+m1(2,1)*m2(1,0)+m1(2,2)*m2(2,0)+m1(2,3)*m2(3,0),
      m1(2,0)*m2(0,1)+m1(2,1)*m2(1,1)+m1(2,2)*m2(2,1)+m1(2,3)*m2(3,1),
      m1(2,0)*m2(0,2)+m1(2,1)*m2(1,2)+m1(2,2)*m2(2,2)+m1(2,3)*m2(3,2),
      m1(2,0)*m2(0,3)+m1(2,1)*m2(1,3)+m1(2,2)*m2(2,3)+m1(2,3)*m2(3,3),

      m1(3,0)*m2(0,0)+m1(3,1)*m2(1,0)+m1(3,2)*m2(2,0)+m1(3,3)*m2(3,0),
      m1(3,0)*m2(0,1)+m1(3,1)*m2(1,1)+m1(3,2)*m2(2,1)+m1(3,3)*m2(3,1),
      m1(3,0)*m2(0,2)+m1(3,1)*m2(1,2)+m1(3,2)*m2(2,2)+m1(3,3)*m2(3,2),
      m1(3,0)*m2(0,3)+m1(3,1)*m2(1,3)+m1(3,2)*m2(2,3)+m1(3,3)*m2(3,3));
  return nrvo;
}

Real Angle(const Vec4& v1, const Vec4& v2)
{
  Real t=std::min(std::max(Dot(v1,v2)/(Length(v1)*Length(v2)),Real(-1.0)),Real(1.0));
  return std::acos(t);
}

} // ns
