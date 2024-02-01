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
  Author: Ansgar Philippsen
*/

#include "transformation_base.hh"

namespace ost { namespace img { namespace alg {

Transformation::Transformation():
  tmat_(),imat_()
{}

Transformation::Transformation(const Mat4& tmat):
  tmat_(tmat),imat_(Invert(tmat)) 
{}

Transformation::Transformation(const Mat4& tmat, const Mat4& imat):
  tmat_(tmat),imat_(imat) 
{}

Point Transformation::Apply(const Point& p)
{
  return Point(Vec4(this->Matrix()*p.ToVec4()));
}

Extent Transformation::Apply(const Extent& e)
{
  /*
    apply matrix to all 8 points of extent and
    return minimal and maximal values in new extent
  */
  Vec4 corner[8], tcorner[8];

  // assume index is middle of pixel

  Vec3 p1(e.GetStart().ToVec3());
  Vec3 p2(e.GetEnd().ToVec3());

  Vec3 offset(0,0,0);

  switch(e.GetDim()) {
  case 3:
    offset[2]=0.5;
    // fallthrough on purpose
  case 2:
    offset[1]=0.5;
    // fallthrough on purpose
  case 1:
    offset[0]=0.5;
    break;
  }

  p1-=offset;
  p2+=offset;

  corner[0] = Vec4(p1[0],p1[1],p1[2],1.0);
  corner[1] = Vec4(p2[0],p1[1],p1[2],1.0);
  corner[2] = Vec4(p1[0],p2[1],p1[2],1.0);
  corner[3] = Vec4(p1[0],p1[1],p2[2],1.0);
  corner[4] = Vec4(p2[0],p2[1],p1[2],1.0);
  corner[5] = Vec4(p2[0],p1[1],p2[2],1.0);
  corner[6] = Vec4(p1[0],p2[1],p2[2],1.0);
  corner[7] = Vec4(p2[0],p2[1],p2[2],1.0);
 
  for(int i=0;i<8;i++) {
    tcorner[i] = Vec4(this->Matrix()*corner[i]);
  }

  Vec3 cmin((tcorner[0])[0],(tcorner[0])[1],(tcorner[0])[2]);
  Vec3 cmax = cmin;
  for(int i=1;i<8;i++) {
    for(int k=0;k<3;++k) {
      cmin[k] = std::min(tcorner[i][k],cmin[k]);
      cmax[k] = std::max(tcorner[i][k],cmax[k]);
    }
  }

  return Extent(Point(Vec3(cmin+offset)),Point(Vec3(cmax-offset)));
}

}}}

