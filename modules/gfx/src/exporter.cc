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

#include <ost/gfx/scene.hh>
#include "exporter.hh"

namespace ost { namespace gfx {

  void Exporter::SetupTransform(const Scene* scene)
  {
    vertex_tf_=geom::Mat4();
    normal_tf_=geom::Mat3();
    geom::Mat4 Tmat,Rmat,Smat,Cmat;
    if(to_origin_) {
      geom::Vec3 cen=scene->GetTransform().GetCenter();
      vertex_tf_=geom::Mat4(scene->GetTransform().GetRot()) *
                 geom::Mat4(scale_, 0.0, 0.0, 0.0,
                            0.0, scale_, 0.0, 0.0,
                            0.0, 0.0, scale_, 0.0,
                            0.0, 0.0, 0.0, 1.0)*
                 geom::Mat4(1.0,0.0,0.0,-cen[0],
                           0.0,1.0,0.0,-cen[1],
                           0.0,0.0,1.0,-cen[2],
                           0.0,0.0,0.0,1.0);
      normal_tf_=scene->GetTransform().GetRot();
    } else {
      vertex_tf_=geom::Mat4(scale_, 0.0, 0.0, 0.0,
                            0.0, scale_, 0.0, 0.0,
                            0.0, 0.0, scale_, 0.0,
                            0.0, 0.0, 0.0, 1.0);
    }
  }

  void Exporter::TransformPosition(float* v) const
  {
    geom::Vec3 result(vertex_tf_*geom::Vec4(v[0],v[1],v[2],1.0));
    v[0]=result[0]; v[1]=result[1]; v[2]=result[2];
  }

  void Exporter::TransformNormal(float* n) const
  {
    geom::Vec3 result=normal_tf_*geom::Vec3(n[0],n[1],n[2]);
    n[0]=result[0]; n[1]=result[1]; n[2]=result[2];
  }
}} // ns
