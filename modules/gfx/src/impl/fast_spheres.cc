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

/*
  Author: Ansgar Philippsen
*/

#include <ost/config.hh>
#if OST_SHADER_SUPPORT_ENABLED
#include <ost/gfx/shader.hh>
#endif
#include <ost/gfx/scene.hh>

#include "fast_spheres.hh"

namespace ost { namespace gfx { namespace impl {

  FastSphereRenderer::FastSphereRenderer(size_t reserve):
    data_()
  {
    data_.reserve(reserve);
  }

  FastSphereRenderer::FastSphereRenderer(const DataList& d):
    data_(d)
  {}

  void FastSphereRenderer::Clear()
  {
    data_.clear();
  }

  void FastSphereRenderer::Add(float pos[3], float col[4], float rad)
  {
    static VData data;
    data.pos[0]=pos[0]; data.pos[1]=pos[1]; data.pos[2]=pos[2];
    data.col[0]=col[0]; data.col[1]=col[1]; data.col[2]=col[2]; data.col[3]=col[3];
    data.rad=rad;
    data_.push_back(data);
  }

  void FastSphereRenderer::Add(const geom::Vec3& pos, const Color& col, float rad)
  {
    static VData data;
    data.pos[0]=pos[0]; data.pos[1]=pos[1]; data.pos[2]=pos[2];
    data.col[0]=col.Red(); data.col[1]=col.Green(); data.col[2]=col.Blue(); data.col[3]=col.Alpha();
    data.rad=rad;
    data_.push_back(data);
  }

  void FastSphereRenderer::RenderGL()
  {
#if OST_SHADER_SUPPORT_ENABLED
    geom::Mat3 irot=geom::Transpose(Scene::Instance().GetTransform().GetRot());
    geom::Vec3 cx=irot*geom::Vec3(1.0,0.0,0.0);
    geom::Vec3 cy=irot*geom::Vec3(0.0,1.0,0.0);

    uint write_normals = Shader::Instance().GetCurrentName()=="dumpnorm" ? 1 : 0;

    Shader::Instance().PushProgram();
    if(Shader::Instance().GetCurrentName()=="hemilight" ||
       Shader::Instance().GetCurrentName()=="hf") {
      Shader::Instance().Activate("fast_sphere_hf");
      GLuint cpr=Shader::Instance().GetCurrentProgram();
      GLint loc=glGetUniformLocation(cpr,"hemi_param");
      if(loc>=0) {
        geom::Vec4 hp=Scene::Instance().GetHemiParams();
        glUniform4f(loc,hp[0],hp[1],hp[2],hp[3]);
      }
    } else if(Shader::Instance().GetCurrentName()=="toon1") {
      Shader::Instance().Activate("fast_sphere_toon1");
    } else if(Shader::Instance().GetCurrentName()=="toon2") {
      Shader::Instance().Activate("fast_sphere_toon2");
    } else {
      Shader::Instance().Activate("fast_sphere_phong");
    }
    Shader::Instance().UpdateState();
    glUniform1i(glGetUniformLocation(Shader::Instance().GetCurrentProgram(),"write_normals"),write_normals);

    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);

    GLdouble gl_mmat[16];
    glGetDoublev(GL_MODELVIEW_MATRIX,gl_mmat);
    GLdouble gl_pmat[16];
    glGetDoublev(GL_PROJECTION_MATRIX,gl_pmat);

    glBegin(GL_QUADS);
    for(DataList::const_iterator it=data_.begin();it!=data_.end();++it) {
      // eye z
      float zz=gl_mmat[2]*it->pos[0]+gl_mmat[6]*it->pos[1]+gl_mmat[10]*it->pos[2]+gl_mmat[14];
      // difference in z of sphere center and sphere pole
      float scale=0.5*gl_pmat[14]*it->rad/(zz*(-it->rad+zz));
      
      // draw correctly scaled billboards
      glColor3fv(it->col);
      float px=it->rad*cx[0];
      float py=it->rad*cx[1];
      float pz=it->rad*cx[2];
      float qx=it->rad*cy[0];
      float qy=it->rad*cy[1];
      float qz=it->rad*cy[2];
      glTexCoord4f(0.0,0.0,-it->rad,scale);
      glVertex3f(it->pos[0]-px-qx,it->pos[1]-py-qy,it->pos[2]-pz-qz);
      glTexCoord4f(0.0,1.0,-it->rad,scale);
      glVertex3f(it->pos[0]-px+qx,it->pos[1]-py+qy,it->pos[2]-pz+qz);
      glTexCoord4f(1.0,1.0,-it->rad,scale);
      glVertex3f(it->pos[0]+px+qx,it->pos[1]+py+qy,it->pos[2]+pz+qz);
      glTexCoord4f(1.0,0.0,-it->rad,scale);
      glVertex3f(it->pos[0]+px-qx,it->pos[1]+py-qy,it->pos[2]+pz-qz);
    }
    glEnd();

    glPopAttrib();

    Shader::Instance().PopProgram();
    Shader::Instance().UpdateState();
#endif
  }


}}} // ns
