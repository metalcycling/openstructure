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

#include "vecmat2_op.hh"
#include "vecmat3_op.hh"
#include "vec3.hh"

#include "vec2.hh"
#include "mat2.hh"

namespace geom {
// determinant
Real Det(const Mat2& m)
{
  return m(0,0)*m(1,1)-m(0,1)*m(1,0);
}

Mat2 Transpose(const Mat2& m)
{
  Mat2 r;
  for (int i=0;i<2;i++){
    for (int j=0;j<2;j++){
      r(j,i)=m(i,j);
    }
  }
  return r;
}
Mat2 Invert(const Mat2& m)
{
  Real d=Det(m);
  return Mat2(m(1,1),-m(0,1),-m(1,0),m(0,0))/d;    
}

Real Angle(const Vec2& v1, const Vec2& v2)
{
  return std::acos(std::max(static_cast<Real>(-1.0),std::min(static_cast<Real>(1.0),Dot(v1,v2)/Length(v1)/Length(v2))));
}

Real SignedAngle(const Vec2& v1, const Vec2& v2)
{
  if(Length(v1)==0.0 || Length(v2)==0.0){
    return 0.0;
  }
  Vec3 vc=Cross(Vec3(v1),Vec3(v2));
  Real sign=1.0;
  if(vc.z!=0.0){
    sign=vc.z/std::fabs(vc.z);
  }
  return acos(std::max(static_cast<Real>(-1.0),std::min(static_cast<Real>(1.0),Dot(v1,v2)/Length(v1)/Length(v2))))*sign;
}

Vec2 Rotate(const Vec2& v,Real ang)
{
  return Mat2(cos(ang),-sin(ang),sin(ang),cos(ang))*v;
}

} // ns
