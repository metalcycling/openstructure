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

#include <ost/config.hh>

#include "transform.hh"
#include "vecmat3_op.hh"
#include "vecmat4_op.hh"

namespace geom {

Transform::Transform():
  rot_(Mat3(1,0,0, 0,1,0, 0,0,1)),
  trans_(Vec3(0,0,0)),
  cen_(Vec3(0,0,0)),
  tm_(),
  ttm_()
{
  update_tm();
}
    
void Transform::SetMatrix(const Mat4& m)
{
  tm_=m;
  ttm_ = Transpose(tm_);
  try {
    itm_ = Invert(tm_);
  } catch (GeomException& e) {
    std::cerr << "caught GeomException in Transform::SetMatrix: " << e.what() << std::endl;
    std::cerr << m << std::endl;
    itm_=geom::Mat4();
  }
  update_components();
}

void Transform::SetTrans(const Vec3& t) 
{
  trans_=t;
  update_tm();
}
    
Vec3 Transform::GetTrans() const 
{
  return trans_;
}
    
void Transform::SetCenter(const Vec3& c) 
{
  cen_=c;
  update_tm();
}
    
Vec3 Transform::GetCenter() const 
{
  return cen_;
}
    
Mat3 Transform::GetRot() const 
{
  return rot_;
}
    
void Transform::SetRot(const Mat3& r) 
{
  rot_=r;
  update_tm();
}

namespace {

static const float P_180 = M_PI/180.0;

}

geom::Mat3 Transform::GetXAxisRotation(float delta)
{
  return AxisRotation(Vec3(rot_(0,0),rot_(0,1),rot_(0,2)),delta*P_180);
}

geom::Mat3 Transform::GetYAxisRotation(float delta)
{
  return AxisRotation(Vec3(rot_(1,0),rot_(1,1),rot_(1,2)),delta*P_180);
}

geom::Mat3 Transform::GetZAxisRotation(float delta)
{
  return AxisRotation(Vec3(rot_(2,0),rot_(2,1),rot_(2,2)),-delta*P_180);
}

void Transform::ApplyXAxisRotation(float delta)
{
  rot_ = rot_*GetXAxisRotation(delta);
  update_tm();
}

void Transform::ApplyYAxisRotation(float delta)
{
  rot_ = rot_*GetYAxisRotation(delta);
  update_tm();
}
    
void Transform::ApplyZAxisRotation(float delta)
{
  rot_ = rot_*GetZAxisRotation(delta);
  update_tm();
}

void Transform::ApplyAxisRotation(float delta, const Vec3& axis)
{
  rot_=rot_*AxisRotation(axis, delta*P_180);
  update_tm();
}

void Transform::ApplyXAxisTranslation(float delta)
{
  trans_[0]+=delta;
  update_tm();
}

void Transform::ApplyYAxisTranslation(float delta)
{
  trans_[1]+=delta;
  update_tm();
}

void Transform::ApplyZAxisTranslation(float delta)
{
  trans_[2]+=delta;
  update_tm();
}
    
Vec3 Transform::Apply(const Vec3& v) const
{
  Vec3 nrvo(tm_*Vec4(v));
  return nrvo;
}

Vec4 Transform::Apply(const Vec4& v) const
{
  Vec4 nrvo=tm_*v;
  return nrvo;
}

Vec3 Transform::ApplyInverse(const Vec3& v) const
{
  Vec3 nrvo(itm_*Vec4(v));
  return nrvo;
}

Vec4 Transform::ApplyInverse(const Vec4& v) const
{
  Vec4 nrvo=itm_*v;
  return nrvo;
}

geom::AlignedCuboid Transform::Apply(const geom::AlignedCuboid& c) const
{
  geom::Vec3 cmin=c.GetMin();
  geom::Vec3 cmax=c.GetMax();
  Vec3 t1 = Apply(Vec3(cmin[0],cmin[1],cmin[2]));
  Vec3 t2 = Apply(Vec3(cmin[0],cmax[1],cmin[2]));
  Vec3 t3 = Apply(Vec3(cmax[0],cmax[1],cmin[2]));
  Vec3 t4 = Apply(Vec3(cmax[0],cmin[1],cmin[2]));
  Vec3 t5 = Apply(Vec3(cmin[0],cmin[1],cmax[2]));
  Vec3 t6 = Apply(Vec3(cmin[0],cmax[1],cmax[2]));
  Vec3 t7 = Apply(Vec3(cmax[0],cmax[1],cmax[2]));
  Vec3 t8 = Apply(Vec3(cmax[0],cmin[1],cmax[2]));
  geom::Vec3 minc = Min(t1,Min(t2,Min(t3,Min(t4,Min(t5,Min(t6,Min(t7,t8)))))));
  geom::Vec3 maxc = Max(t1,Max(t2,Max(t3,Max(t4,Max(t5,Max(t6,Max(t7,t8)))))));
  return geom::AlignedCuboid(minc,maxc);
}

Transform Transform::Apply(const Transform& tf) const
{
  Transform nrvo(*this);
  nrvo.SetMatrix(tf.GetMatrix()*nrvo.GetMatrix());
  return nrvo;
}

/*
  The order of the transformations given herein is conceptually
  "backward" as they are applied to a vertex, because the left-right
  notation is opposite to the order that the vertex experiences
  the transformations (see Ch 3 in the Red Book for more detail)

    v' = T R S C v

     v' := transformed vector
     T  := translation matrix
     R  := rotation matrix
     S  := scaling matrix
     C  := centering matrix
     v  := original vector

*/

void Transform::update_tm()
{
  tm_ =
    Mat4(1.0,0.0,0.0,trans_[0],
         0.0,1.0,0.0,trans_[1],
         0.0,0.0,1.0,trans_[2],
         0.0,0.0,0.0,1.0) *
    Mat4(rot_) *
    Mat4(1.0,0.0,0.0,-cen_[0],
         0.0,1.0,0.0,-cen_[1],
         0.0,0.0,1.0,-cen_[2],
         0.0,0.0,0.0,1.0);
  ttm_ = Transpose(tm_);
  // TODO: calculate from rot, cen and trans
  try {
    itm_ = Invert(tm_);
  } catch (GeomException& e) {
    std::cerr << "caught GeomException in Transform::update_tm: " << e.what() << std::endl;
    itm_=geom::Mat4();
  }
}

void Transform::update_components()
{
  // there is no way to extract the centering component
  // so we just get a rotation and translation
  rot_ = tm_.ExtractRotation();
  trans_ = tm_.ExtractTranslation();
  cen_ = Vec3(0,0,0);
}

} // ns
