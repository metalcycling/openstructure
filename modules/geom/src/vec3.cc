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
#include <algorithm>

#include <Eigen/SVD>

#include "vec3.hh"
// TODO: these are for the (misplaced) Vec3List algorithm functions
#include "vecmat3_op.hh"
#include "composite3.hh"
#include "composite3_op.hh"
#include "mat4.hh"

namespace geom {


#if OST_DOUBLE_PRECISION
typedef Eigen::Matrix3d EMat3;
#else
typedef Eigen::Matrix3f EMat3;
#endif


Mat3 Vec3List::GetInertia() const
{
  Mat3 cov(0,0,0,0,0,0,0,0,0);
  Vec3 center=this->GetCenter();
  for (Vec3List::const_iterator i=this->begin(),e=this->end(); i!=e; ++i) {
    Vec3 p=*i-center;
    cov(0,0)+=p.y*p.y+p.z*p.z;
    cov(1,1)+=p.x*p.x+p.z*p.z;
    cov(2,2)+=p.x*p.x+p.y*p.y;
    cov(0,1)-=p.x*p.y;
    cov(1,2)-=p.y*p.z;
    cov(0,2)-=p.x*p.z;
  }  
  cov(1,0)=cov(0,1);    
  cov(2,1)=cov(1,2);    
  cov(2,0)=cov(0,2);  
  return cov;
}

Mat3 Vec3List::GetPrincipalAxes() const
{
  Mat3 inertia=this->GetInertia();  
  EMat3 inertia_mat(inertia.Data());

  Eigen::JacobiSVD<EMat3> svd(inertia_mat,Eigen::ComputeFullU);
  EMat3 rot=svd.matrixU();
  Mat3 axes(rot.data());
  return axes;
}

Vec3 Vec3List::GetCenter() const
{
  Vec3 center;
  if (this->empty()) {
    return center;
  }
  for (Vec3List::const_iterator i=this->begin(),e=this->end(); i!=e; ++i) {
    center+=*i;
  }
  return center/=this->size();
}

Line3 Vec3List::GetODRLine() const
{
  Vec3 center=this->GetCenter();
  Vec3 direction=this->GetPrincipalAxes().GetRow(2);
  return Line3(center,center+direction);
}

Plane Vec3List::GetODRPlane() const
{
  Vec3 origin=this->GetCenter();
  Vec3 normal=this->GetPrincipalAxes().GetRow(0);
  return Plane(origin,normal);
}

void Vec3List::ApplyTransform(const Mat4& m) 
{
  Real x,y,z;
  for(uint i = 0; i < this->size(); ++i) {
    x = (*this)[i][0];
    y = (*this)[i][1];
    z = (*this)[i][2];
    (*this)[i][0] = m(0,0)*x+m(0,1)*y+m(0,2)*z+m(0,3);
    (*this)[i][1] = m(1,0)*x+m(1,1)*y+m(1,2)*z+m(1,3);
    (*this)[i][2] = m(2,0)*x+m(2,1)*y+m(2,2)*z+m(2,3);
  }
}

Real Vec3List::GetSummedSquaredDistances(const Vec3List& other) const
{
  if(this->size() != other.size()) {
    String m = "Inconsistent sizes in Vec3List::GetSummedSquaredDistances";
    throw GeomException(m);
  }
  Real summed_squared_dist = 0.0;
  for(uint i = 0; i < this->size(); ++i) {
    summed_squared_dist += geom::Length2((*this)[i] - other[i]);
  }
  return summed_squared_dist;
}

Real Vec3List::GetRMSD(const Vec3List& other) const
{
  if(this->empty()) {
    return 0.0;
  } else {
    Real summed_squared_dist = this->GetSummedSquaredDistances(other);
    return std::sqrt(summed_squared_dist/this->size());
  }
}

Real Vec3List::GetGDTHA(const Vec3List& other, bool norm) const
{
  if(this->size() != other.size()) {
    String m = "Inconsistent sizes in Vec3List::GetNWithin";
    throw GeomException(m);
  }
  int n = 0;
  for(uint i = 0; i < this->size(); ++i) {
    Real squared_d = geom::Length2((*this)[i] - other[i]);
    // GDTHA default thresholds: [4,2,1,0.5]
    if(squared_d < static_cast<Real>(16)) {
      ++n;
    }
    if(squared_d < static_cast<Real>(4)) {
      ++n;
    }
    if(squared_d < static_cast<Real>(1)) {
      ++n;
    }
    if(squared_d < static_cast<Real>(0.25)) {
      ++n;
    }
  }
  return norm && !this->empty() ? static_cast<Real>(n)/(this->size()*4) : n;
}

Real Vec3List::GetGDTTS(const Vec3List& other, bool norm) const
{
  if(this->size() != other.size()) {
    String m = "Inconsistent sizes in Vec3List::GetNWithin";
    throw GeomException(m);
  }
  int n = 0;
  for(uint i = 0; i < this->size(); ++i) {
    Real squared_d = geom::Length2((*this)[i] - other[i]);
    // GDTTS default thresholds: [8,4,2,1]
    if(squared_d < static_cast<Real>(64)) {
      ++n;
    }
    if(squared_d < static_cast<Real>(16)) {
      ++n;
    }
    if(squared_d < static_cast<Real>(4)) {
      ++n;
    }
    if(squared_d < static_cast<Real>(1)) {
      ++n;
    }
  }
  return norm && !this->empty() ? static_cast<Real>(n)/(this->size()*4) : n;
}

Real Vec3List::GetGDT(const Vec3List& other, Real thresh, bool norm) const
{
  if(this->size() != other.size()) {
    String m = "Inconsistent sizes in Vec3List::GetNWithin";
    throw GeomException(m);
  }
  int n = 0;
  Real squared_thresh = thresh * thresh;
  for(uint i = 0; i < this->size(); ++i) {
    Real squared_d = geom::Length2((*this)[i] - other[i]);
    if(squared_d < squared_thresh) {
      ++n;
    }
  }
  return norm && !this->empty() ? static_cast<Real>(n)/(this->size()) : n;
}

Real Vec3List::GetMinDist(const Vec3List& other) const {
  Real min = std::numeric_limits<Real>::max();
  for(size_t i = 0; i < this->size(); ++i) {
    for(size_t j = 0; j < other.size(); ++j) {
      min = std::min(min, geom::Length2((*this)[i] - other[j]));
    }
  }
  return std::sqrt(min);
}

bool Vec3List::IsWithin(const Vec3List& other, Real dist) const {
  Real squared_dist = dist*dist;
  for(size_t i = 0; i < this->size(); ++i) {
    for(size_t j = 0; j < other.size(); ++j) {
      if(geom::Length2((*this)[i] - other[j]) <= squared_dist) {
        return true;
      }
    }
  }
  return false;
}

std::pair<Line3, Real> Vec3List::FitCylinder(const Vec3& initial_direction) const
  { 
    Vec3 center=this->GetCenter();
    Line3 axis=Line3(center,center+initial_direction), axis_old;
    Real radius,res_sum_old,res_sum,delta_0=0.01,prec=0.0000001,err,norm,delta;
    unsigned long n_step=1000, n_res=this->size();
    Vec3 v,gradient_dir,gradient_center;
    
    radius=0.0;
    delta=delta_0;
    for (Vec3List::const_iterator i=this->begin(),e=this->end(); i!=e; ++i) {
      radius+=geom::Distance(axis,(*i));
    }
    radius/=Real(n_res);
    res_sum=0.0;
    for (Vec3List::const_iterator i=this->begin(),e=this->end(); i!=e; ++i) {
      Real r=Distance(axis,(*i))-radius;
      res_sum+=r*r;
    }
    unsigned long k=0;
    err=2.0*prec;
    while (err>prec && k<n_step) {
      res_sum_old=res_sum;
      axis_old=axis;
      //radius=0.0;
      if (k>50) {
        delta=delta_0*50.0*50.0/(k*k);
      }
      //for (Vec3List::const_iterator i=this->begin(),e=this->end(); i!=e; ++i) {
      //  radius+=Distance(axis,(*i));
      //}
      radius/=Real(n_res);
      for (int j=0; j!=3; ++j){
        res_sum=0.0;
        v=Vec3(0.0,0.0,0.0);
        v[j]=delta;
        axis=Line3(axis_old.GetOrigin(),axis_old.GetOrigin()+axis_old.GetDirection()+v);
        radius=0.0;
        for (Vec3List::const_iterator i=this->begin(),e=this->end(); i!=e; ++i) {
          radius+=Distance(axis,(*i));
        }
        radius/=Real(n_res);
        for (Vec3List::const_iterator i=this->begin(),e=this->end(); i!=e; ++i) {
          Real r=Distance(axis,(*i))-radius;
          res_sum+=r*r;
        }
        gradient_dir[j]=(res_sum-res_sum_old)/delta;
      }
      norm=Dot(gradient_dir,gradient_dir);
      if (norm>1.) {
        gradient_dir=Normalize(gradient_dir);
      }
      for (int j=0; j!=3; ++j){
        res_sum=0.0;
        v=Vec3(0.0,0.0,0.0);
        v[j]=delta;
        axis=Line3(axis_old.GetOrigin()+v,axis_old.GetOrigin()+axis_old.GetDirection()+v);
        radius=0.0;
        for (Vec3List::const_iterator i=this->begin(),e=this->end(); i!=e; ++i) {
          radius+=Distance(axis,(*i));
        }
        radius/=Real(n_res);
        for (Vec3List::const_iterator i=this->begin(),e=this->end(); i!=e; ++i) {
          Real r=Distance(axis,(*i))-radius;
          res_sum+=r*r;
        }
        gradient_center[j]=(res_sum-res_sum_old)/delta;
      }
      norm=Dot(gradient_center,gradient_center);
      if (norm>1.) {
        gradient_center=Normalize(gradient_center);
      }      
      axis=Line3(axis_old.GetOrigin()-50*delta*gradient_center,axis_old.GetOrigin()-50*delta*gradient_center+axis_old.GetDirection()-delta*gradient_dir);
      radius=0.0;
      for (Vec3List::const_iterator i=this->begin(),e=this->end(); i!=e; ++i) {
        radius+=Distance(axis,(*i));
      }
      radius/=Real(n_res);
      res_sum=0.0;
      for (Vec3List::const_iterator i=this->begin(),e=this->end(); i!=e; ++i) {
        Real r=Distance(axis,(*i))-radius;
        res_sum+=r*r;
      }
      err=fabs((res_sum-res_sum_old)/float(n_res));
      k++;
    }
    if (err>prec) {
      std::cout<<"axis fitting did not converge"<<std::endl;
    }
    return std::make_pair(axis,radius);
  }
}
