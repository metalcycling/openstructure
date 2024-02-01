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
  volatile test object for current hacking purposes

  Author: Ansgar Philippsen
*/

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

#include <ost/platform.hh>
#include "texture.hh"
#include "glext_include.hh"
#include "gfx_test_object.hh"
#include "scene.hh"
#include "bitmap_io.hh"

namespace ost { namespace gfx {

namespace bf = boost::filesystem;

namespace {

static const GLfloat coords[][3] = {
  {-5.0,-3.0,-1.0},
  {-5.0,3.0,-1.0},
  {5.0,3.0,-1.0},
  {5.0,-3.0,-1.0},
  {-5.0,-3.0,1.0},
  {-5.0,3.0,1.0},
  {5.0,3.0,1.0},
  {5.0,-3.0,1.0},
};

}

GfxTestObj::GfxTestObj():
  GfxObj("GfxTestObj")
{
  AddLabel("The Quick Brown Fox Jumps Over The Lazy Dog", geom::Vec3(0,0,0),Color(1,1,0));

  glGenTextures(1,&tex_id);
  String ost_root=GetSharedDataPath();
  bf::path ost_root_dir(ost_root);
  bf::path tex_file(ost_root_dir / "textures/test_texture.png");

  Texture tex(ImportBitmap(tex_file.string()));
  if(!tex.IsValid()) {
    LOG_ERROR("error loading " << tex_file.string());
  } else {
    LOG_DEBUG("importing tex with id " << tex_id);
    glBindTexture(GL_TEXTURE_2D, tex_id);
    glPixelStorei(GL_UNPACK_ALIGNMENT,1);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D,0,GL_RGB,tex.width(),tex.height(),tex.border(),tex.format(),tex.type(),tex.data());
  }
}

geom::Vec3 GfxTestObj::GetCenter() const
{
  return geom::Vec3(0,0,0);
}

geom::AlignedCuboid GfxTestObj::GetBoundingBox() const
{
  return geom::AlignedCuboid(geom::Vec3(-5,-5,-5),geom::Vec3(5,5,5));
}

void GfxTestObj::CustomRenderGL(RenderPass pass)
{
  if(pass>0) return;
#if 1
  GLint quad[][4] = {
    {0,1,2,3},
    {4,7,6,5},
    {0,4,5,1},
    {2,6,7,3},
    {0,3,7,4},
    {1,5,6,2}
  };

  GLfloat norm[][3] = {
    {0,0,-1},
    {0,0,1},
    {-1,0,0},
    {1,0,0},
    {0,-1,0},
    {0,1,0}
  };

  GLfloat col[][3] = {
    {1,0,0},
    {0,1,0},
    {0,0,1},
    {1,1,0},
    {1,0,1},
    {0,1,1}
  };

  GLfloat tex[][2] = {
    {0.0,0.0},{1.0,0.0},{1.0,1.0},{0.0,1.0}
  };

  GLfloat shift[][3] = {
    {0.0,1,-3},
    {0.0,-1.0,3}
  };

  glPushAttrib(GL_ENABLE_BIT);
  glEnable(GL_LIGHTING);
  glDisable(GL_CULL_FACE);
  glDisable(GL_NORMALIZE);
  glEnable(GL_TEXTURE_2D);
  glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
  glBindTexture(GL_TEXTURE_2D,tex_id);
  
  glBegin(GL_QUADS);
  for(int h=0;h<2;++h) {
    for(int i=0;i<6;++i) {
      glNormal3fv(norm[i]);
      glColor3fv(col[i]);
      for(int k=0;k<4;++k) {
        glTexCoord2fv(tex[k]);
        glVertex3f(coords[quad[i][k]][0]+shift[h][0],
                   coords[quad[i][k]][1]+shift[h][1],
                   coords[quad[i][k]][2]+shift[h][2]);
      }
    }
  }
  glEnd();
  glPopAttrib();
#endif
}

}}
