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
#include <iostream>

#include "vecmat3_op.hh"
#include "vecmat2_op.hh"

#include "vec3.hh"
#include "mat3.hh"
#include "mat2.hh"
#include "exc.hh"

namespace geom {

Real Comp(const Mat3& m,unsigned int i_in, unsigned int j_in)
{
  return Minor(m,i_in,j_in)* static_cast<Real>( ((i_in+j_in)&0x1) ? -1 : 1);
}

Real Minor(const Mat3& m, unsigned int i_in, unsigned int j_in)
{
  if(i_in>2 || j_in>2){ 
    throw OutOfRangeException();
  }
  Mat2 r;
  unsigned int itmp=0;
  for (unsigned int i=0;i<3;i++) {
    if(i==i_in) continue;
    unsigned int jtmp=0;
    for (unsigned int j=0;j<3;j++){
      if(j==j_in) continue;
      r(itmp,jtmp)=m(i,j);
      jtmp++;
    }
    itmp++;
  }
  
  return Det(r);
}

Mat3 Transpose(const Mat3& m)
{
  Mat3 r;
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      r(j,i)=m(i,j);
    }
  }
  return r;
}

Mat3 Invert(const Mat3& m)
{
  Mat3 r;
  Real d=Det(m);
  if(d<1.0e-30) d=1.0e30;
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      r(j,i)=Comp(m,i,j)/d;
    }
  }
  return r;
}

Real Det(const Mat3& m)
{
#if 0
  return m(0,0)*m(1,1)*m(2,2)-m(0,0)*m(1,2)*m(2,1)-m(0,1)*m(1,0)*m(2,2)+m(0,2)*m(1,0)*m(2,1)-m(0,2)*m(1,1)*m(2,0);
#else
  Real d=0;
  for(int i=0;i<3;i++)
  {
    d+=m(0,i)*Comp(m,0,i);
  }
  return d;
#endif
}

Real Angle(const Vec3& v1, const Vec3& v2)
{
  Real dot_product = Dot(Normalize(v1), Normalize(v2));
  dot_product=std::max(static_cast<Real>(-1.0),dot_product);
  dot_product=std::min(static_cast<Real>(1.0),dot_product);
  return std::acos(dot_product);
}

Real SignedAngle(const Vec3& v1, const Vec3& v2, const Vec3& ref_normal )
{
  Vec3 a=Normalize(v1);
  Vec3 b=Normalize(v2);
  Vec3 c=Normalize(ref_normal);
  return std::atan2(Dot(c,Cross(a,b)), Dot(a,b));
}

Mat3 EulerTransformation(Real theta, Real phi, Real xi)
{
  Real costheta=cos(theta);
  Real sintheta=sin(theta);
  Real cosphi=cos(phi);
  Real sinphi=sin(phi);
  Real cosxi=cos(xi);
  Real sinxi=sin(xi);

  // rotation of xi around z-axis
  Mat3 rmat1(cosxi, sinxi, 0.0,
       -sinxi, cosxi, 0.0,
       0.0, 0.0, 1.0);

  // out of plane rotation of phi around y-axis
  Mat3 rmat2(cosphi, 0.0, sinphi,
       0.0, 1.0, 0.0,
       -sinphi, 0.0, cosphi);

  // xy plane orientation of theta - again rotation around z-axis
  Mat3 rmat3(costheta, -sintheta, 0.0,
       sintheta, costheta, 0.0,
       0.0, 0.0, 1.0);

  /*
    note column vector multiplication ala v2=A B C v
    where C is the first transformation, followed by B and A
  */
  Mat3 nrvo=rmat3*rmat2*rmat1;
  return nrvo;
}

Vec3 OrthogonalVector(const Vec3& vec) {
  if (std::abs(vec[0]) < std::abs(vec[1])) {
    if (std::abs(vec[0]) < std::abs(vec[2]))
      return Normalize(Cross(vec, Vec3(1, 0, 0)+vec));
    else
      return Normalize(Cross(vec, Vec3(0, 0, 1)+vec));
  } else {
    if (std::abs(vec[1]) < std::abs(vec[2]))
      return Normalize(Cross(vec, Vec3(0, 1, 0)+vec));
    else
      return Normalize(Cross(vec, Vec3(0, 0, 1)+vec));
  }
}

Mat3 AxisRotation(const Vec3& axis, Real angle)
{
  Real sa=sin(angle);
  Real ca=cos(angle);
  Real v=Length(axis);
  if(v==0.0) return Mat3();

  Vec3 vv=1.0/v * axis;
  
  Real x=vv[0];
  Real y=vv[1];
  Real z=vv[2];
  Real xx=x*x;
  Real xy=x*y;
  Real xz=x*z;
  Real yy=y*y;
  Real yz=y*z;
  Real zz=z*z;

  return Mat3(xx+ca-xx*ca, xy-ca*xy-sa*z, xz-ca*xz+sa*y,
          xy-ca*xy+sa*z, yy+ca-ca*yy, yz-ca*yz-sa*x,
          xz-ca*xz-sa*y, yz-ca*yz+sa*x,zz+ca-ca*zz);
}
  
Real MinDistance(const Vec3List& l1, const Vec3List& l2)
{ 
  // returns the minimal distance between two sets of points (Vec3List)
  if (l1.size()==0 || l2.size()==0){throw GeomException("cannot calculate minimal distance: empty Vec3List");}
  Real min=Length2(*l1.begin()-*l2.begin());
  Real d;
  for (Vec3List::const_iterator p1=l1.begin(),e1=l1.end(); p1!=e1; p1++) {
    for (Vec3List::const_iterator p2=l2.begin(),e2=l2.end(); p2!=e2; p2++) {
      d=Length2(*p1-*p2);
      if (d<min) min=d;
    }
  }
  return std::sqrt(min);
}

Real MinDistanceWithPBC(const Vec3List& l1, const Vec3List& l2, Vec3& ucell_size)
{ 
  // returns the minimal distance between two sets of points (Vec3List)
  // given the periodic boundary condition along x,y,z given in the ucell_size
  if (l1.size()==0 || l2.size()==0){throw GeomException("cannot calculate minimal distance: empty Vec3List");}
  Real min=Length2(*l1.begin()-*l2.begin());
  Real d;
  Vec3 v;
  for (int i=0; i<3; i++) {
    ucell_size[i]=std::fabs(ucell_size[i]);
  }
  for (Vec3List::const_iterator p1=l1.begin(),e1=l1.end(); p1!=e1; p1++) {
    for (Vec3List::const_iterator p2=l2.begin(),e2=l2.end(); p2!=e2; p2++) {
      d=Distance2WithPBC(*p1, *p2, ucell_size);
      if (d<min) min=d;
    }
  }
  return std::sqrt(min);
}  

std::vector<unsigned int> MinDistanceIndices(const Vec3List& l1, const Vec3List& l2)
{ 
  // returns the indices index1, index2 corresponding to the points in
  // the Vec3List l1 and l2 having the minimal distance.
  if (l1.size()==0 || l2.size()==0){throw GeomException("cannot calculate minimal distance: empty Vec3List");}
  Real min=Length2(*l1.begin()-*l2.begin());
  Real d;
  std::vector<unsigned int> il;
  il.push_back(0);
  il.push_back(0);
  for (unsigned int i=0;i!=l1.size();i++){
    for (unsigned int j=0;j!=l2.size();j++){
      d=Length2(l1[i]-l2[j]);
      if (d<min){
        min=d;
        il[0]=i;
        il[1]=j;
      }
    }
  }
  return il;
}

Vec3List CalculateUnitCellVectors(const Vec3& ucell_size, const Vec3& ucell_angles){
  // Calculates the unit cell vectors from their sizes and angles.
  // The angles are given as Vec3(gamma,beta,alpha)
  geom::Vec3List ucell_vec;
  geom::Vec3 ucell_vec1,ucell_vec2,ucell_vec3;
  Real cosa=cos(ucell_angles[2]),cosb=cos(ucell_angles[1]);
  Real cosg=cos(ucell_angles[0]),sing=sin(ucell_angles[0]);
  ucell_vec1=geom::Vec3(ucell_size[0],0,0);
  ucell_vec2=ucell_size[1]*geom::Vec3(cosg,sing,0);
  ucell_vec3=ucell_size[2]*geom::Vec3(cosb,(cosa-cosb*cosg)/sing,\
                                       pow(1-(cosa*cosa+cosb*cosb-2.*cosa*cosb*cosg)/(sing*sing),0.5));
  ucell_vec.push_back(ucell_vec1);
  ucell_vec.push_back(ucell_vec2);
  ucell_vec.push_back(ucell_vec3);
  return ucell_vec;  
}
  
Vec3 WrapVec3(const Vec3& v1,const Vec3& box_center,const Vec3& ucell_size){
  Vec3 v;
  Real r;
  for (int i=0; i<3; i++) {
    r=(v1[i]-box_center[i])/ucell_size[i];
    r=(r > 0.0) ? std::floor(r + 0.5) : std::ceil(r - 0.5);
    v[i]=v1[i]-ucell_size[i]*r;
  }
  return v;
}

Vec3List WrapVec3List(const Vec3List& vl, const Vec3& box_center,const Vec3& ucell_size){
  Vec3List vl_out;
  vl_out.reserve(vl_out.size());
  for (Vec3List::const_iterator v1=vl.begin(),e=vl.end();v1!=e ; v1++) {
    vl_out.push_back(WrapVec3(*v1,box_center,ucell_size));
  }
  return vl_out;
}
  
Vec3 WrapVec3(const Vec3& v1,const Vec3& box_center,const Vec3& ucell_size,const Vec3& ucell_angles){
  Vec3List vl=CalculateUnitCellVectors(ucell_size,ucell_angles);
  Vec3 v=v1-box_center,v_wrapped=v1,r;
  for (int i=0; i<3; i++) {
    r[2-i]=v[2-i]/vl[2-i][2-i];
    v=v-r[2-i]*vl[2-i];
    r[2-i]=(r[2-i] > 0.0) ? std::floor(r[2-i] + 0.5) : std::ceil(r[2-i] - 0.5);
    v_wrapped=v_wrapped-vl[2-i]*r[2-i];
  }
  return v_wrapped;
}

Vec3List WrapVec3List(const Vec3List& vl, const Vec3& box_center,const Vec3& ucell_size,const Vec3& ucell_angles){
  Vec3List vl_out;
  vl_out.reserve(vl_out.size());
  for (Vec3List::const_iterator v1=vl.begin(),e=vl.end();v1!=e ; v1++) {
    vl_out.push_back(WrapVec3(*v1,box_center,ucell_size,ucell_angles));
  }
  return vl_out;
}  
  
} // ns
