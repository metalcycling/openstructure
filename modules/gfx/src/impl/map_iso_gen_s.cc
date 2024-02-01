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
#include <cstring>

#include <boost/shared_array.hpp>

#include <ost/profile.hh>

#include <ost/img/image_state.hh>

#include <ost/gfx/glext_include.hh>
#include <ost/gfx/map_iso_prop.hh>
#include "map_iso_gen_s.hh"
#if OST_SHADER_SUPPORT_ENABLED
# define MAP_ISO_USE_SHADER 0
# include <ost/gfx/shader.hh>
# include <ost/gfx/scene.hh>
#endif

namespace ost { namespace gfx { namespace impl {

// singleton implementation
IsosurfaceGeneratorS& IsosurfaceGeneratorS::Instance()
{
  static IsosurfaceGeneratorS inst;
  return inst;
}

IsosurfaceGeneratorS::IsosurfaceGeneratorS()
{
#if MAP_ISO_USE_SHADER
  slab_tex_[0]=0;
#endif
}

namespace {

struct Box
{
  float fr3,fr0,fr4;
#if MAP_ISO_USE_SHADER
  float code;
#else
  unsigned char code;
#endif
};

typedef boost::shared_array<Box> BSlab;
typedef boost::shared_array<float> XYSlab;
typedef boost::shared_array<unsigned char> CSlab;
typedef boost::shared_array<uint> VSlab;

void set_slab(img::RealSpatialImageState* mh, XYSlab& slab, uint ucount, 
              uint vcount, uint w)
{
  for(uint v=0;v<vcount-1;++v) {
    uint off=v*ucount;
    for(uint u=0;u<ucount-1;++u,++off) {
      slab[off]=mh->Value(img::Point(u,v,w));
    }
    slab[off]=0.0;
  }

  for(uint u=0;u<ucount;++u) {
    slab[ucount*(vcount-1)+u]=0.0;
  }
}

BSlab process_edges(const XYSlab& slab1, const XYSlab& slab2, 
                    uint ucount, uint vcount,float level, 
                    uint slab_tex[2], uint slab_id[2], uint w)
{
  BSlab bslab(new Box[ucount*vcount]);
#if MAP_ISO_USE_SHADER
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D,slab_tex[0]);
  glTexImage2D(GL_TEXTURE_2D,0,GL_LUMINANCE,ucount,vcount,0,
               GL_LUMINANCE,GL_FLOAT,slab1.get());
  glUniform1i(slab_id[0],0);

  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D,slab_tex[1]);
  glTexImage2D(GL_TEXTURE_2D,0,GL_LUMINANCE,ucount,vcount,0,
               GL_LUMINANCE,GL_FLOAT,slab2.get());
  glUniform1i(slab_id[1],1);

  glClear(GL_COLOR_BUFFER_BIT);
  glBegin(GL_QUADS);
  glVertex2i(0,0);
  glVertex2i(1,0);
  glVertex2i(1,1);
  glVertex2i(0,1);
  glEnd();

  glReadPixels(0,0,ucount,vcount,GL_RGBA,GL_FLOAT,bslab.get());

#else
  memset(bslab.get(),0,sizeof(Box)*ucount*vcount);
  for(uint v=0;v<vcount-1;++v) {
    uint off=v*ucount;
    for(uint u=0;u<ucount-1;++u,++off) {
      float c0=slab1[off];
      float c1=slab1[off+ucount];
      float c3=slab1[off+1];
      float c4=slab2[off];
     
      if(level>c0 && level<=c3) {
        bslab[off].fr3=(level-c0)/(c3-c0);
      } else if(level>c3 && level<=c0) {
        bslab[off].fr3=(level-c0)/(c3-c0);
      } 

      if(level>c0 && level<=c1) {
        bslab[off].fr0=(level-c0)/(c1-c0);
      } else if(level>c1 && level<=c0) {
        bslab[off].fr0=(level-c0)/(c1-c0);
      } 
      if(level>c0 && level<=c4) {
        bslab[off].fr4=(level-c0)/(c4-c0);
      } else if(level>c4 && level<=c0) {
        bslab[off].fr4=(level-c0)/(c4-c0);
      } 

      bslab[off].code=0;
      if(c0<level) bslab[off].code+=1;
      if(c1<level) bslab[off].code+=2;
      if(slab1[off+ucount+1]<level) bslab[off].code+=4;
      if(c3<level) bslab[off].code+=8;
      if(c4<level) bslab[off].code+=16;
      if(slab2[off+ucount]<level) bslab[off].code+=32;
      if(slab2[off+ucount+1]<level) bslab[off].code+=64;
      if(slab2[off+1]<level) bslab[off].code+=128;
    }
  }
#endif

  return bslab;
}

void process_vids(const img::MapHandle& mh, IndexedVertexArray& va,
                  const BSlab& bslab, VSlab& vidu, VSlab& vidv, VSlab& vidw,
                  uint ucount, uint vcount, uint w)
{
  // u edges
  for(uint v=0;v<vcount;++v) {
    uint off=v*ucount;
    for(uint u=0;u<ucount-1;++u,++off) {
      if(bslab[off].fr3>0.0) {
        vidu[off]=va.Add(mh.FractionalIndexToCoord(geom::Vec3(static_cast<float>(u)+bslab[off].fr3,
                                               static_cast<float>(v),
                                               static_cast<float>(w))),
                         geom::Vec3(),Color(1,1,0));
      }
    }
  }

  // v edges
  for(uint v=0;v<vcount-1;++v) {
    uint off=v*ucount;
    for(uint u=0;u<ucount;++u,++off) {
      if(bslab[off].fr0>0.0) {
        vidv[off]=va.Add(mh.FractionalIndexToCoord(geom::Vec3(static_cast<float>(u),
                                               static_cast<float>(v)+bslab[off].fr0,
                                               static_cast<float>(w))),
                         geom::Vec3(),Color(1,1,0));
      }
    }
  }

  // w edges
  for(uint v=0;v<vcount-1;++v) {
    uint off=v*ucount;
    for(uint u=0;u<ucount-1;++u,++off) {
      if(bslab[off].fr4>0.0) {
        vidw[off]=va.Add(mh.FractionalIndexToCoord(geom::Vec3(static_cast<float>(u),
                                               static_cast<float>(v),
                                               static_cast<float>(w)+bslab[off].fr4)),
                         geom::Vec3(),Color(1,1,0));
      }
    }
  }
}

void assign_lf(IndexedVertexArray& va, const BSlab& bslab,
               const VSlab& vidu1, const VSlab& vidv1, const VSlab& vidw1,
               const VSlab& vidu2, const VSlab& vidv2, const VSlab& vidw2,
               uint ucount, uint vcount)
{
  uint vid[12];
  for(uint v=0;v<vcount-1;++v) {
    for(uint u=0;u<ucount-1;++u) {
      vid[0]=vidv1[(v+0)*ucount+(u+0)];
      vid[1]=vidu1[(v+1)*ucount+(u+0)];
      vid[2]=vidv1[(v+0)*ucount+(u+1)];
      vid[3]=vidu1[(v+0)*ucount+(u+0)];
      
      vid[4]=vidw1[(v+0)*ucount+(u+0)];
      vid[5]=vidw1[(v+1)*ucount+(u+0)];
      vid[6]=vidw1[(v+1)*ucount+(u+1)];
      vid[7]=vidw1[(v+0)*ucount+(u+1)];
      
      vid[8]=vidv2[(v+0)*ucount+(u+0)];
      vid[9]=vidu2[(v+1)*ucount+(u+0)];
      vid[10]=vidv2[(v+0)*ucount+(u+1)];
      vid[11]=vidu2[(v+0)*ucount+(u+0)];
      
#if MAP_ISO_USE_SHADER
      IsosurfaceGenerator::Instance().AddLF(static_cast<unsigned char>(bslab[v*ucount+u].code*255.0),va,vid);
#else
      IsosurfaceGenerator::Instance().AddLF(bslab[v*ucount+u].code,va,vid);
#endif
    }
  }
}

} // anon ns

void IsosurfaceGeneratorS::Generate(const img::MapHandle& mh,
                                   IndexedVertexArray& va, float level)
{
  uint ucount=mh.GetSize()[0]+1;
  uint vcount=mh.GetSize()[1]+1;
  uint wcount=mh.GetSize()[2]+1;
  
  XYSlab slab1(new float[ucount*vcount]);
  memset(slab1.get(),0,sizeof(float)*ucount*vcount);
  XYSlab slab2(new float[ucount*vcount]);
  memset(slab2.get(),0,sizeof(float)*ucount*vcount);

  VSlab vidu1(new uint[ucount*vcount]);
  memset(vidu1.get(),0,sizeof(uint)*ucount*vcount);
  VSlab vidv1(new uint[ucount*vcount]);
  memset(vidv1.get(),0,sizeof(uint)*ucount*vcount);
  VSlab vidw1(new uint[ucount*vcount]);
  memset(vidw1.get(),0,sizeof(uint)*ucount*vcount);

  VSlab vidu2(new uint[ucount*vcount]);
  memset(vidu2.get(),0,sizeof(uint)*ucount*vcount);
  VSlab vidv2(new uint[ucount*vcount]);
  memset(vidv2.get(),0,sizeof(uint)*ucount*vcount);
  VSlab vidw2(new uint[ucount*vcount]);
  memset(vidw2.get(),0,sizeof(uint)*ucount*vcount);

  va.Clear();
  va.Add(geom::Vec3(0,0,0),geom::Vec3(),Color(1,0,0));

  try {
#if MAP_ISO_USE_SHADER
    if(slab_tex_[0]==0) {
      glGenTextures(2,slab_tex_);
    }

    {
    Profile pp("switch GL context");


    OffscreenBuffer::Instance().Resize(ucount,vcount);
    Scene::Instance().InitGL();
    OffscreenBuffer::Instance().Begin();

    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glViewport(0,0,ucount,vcount);

    Shader::Instance().PushProgram();
    Shader::Instance().Activate("iso");

    glDepthMask(0);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
    glEnable(GL_TEXTURE_2D);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glShadeModel(GL_SMOOTH);
    glDrawBuffer(GL_FRONT);
    glReadBuffer(GL_FRONT);

    for(int i=0;i<2;++i) {
      glActiveTexture(GL_TEXTURE0+i);
      glBindTexture(GL_TEXTURE_2D,slab_tex_[i]);
      glPixelStorei(GL_UNPACK_ALIGNMENT,1);
      glPixelStorei(GL_PACK_ALIGNMENT,1);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);  
    }
    slab_id_[0]=glGetUniformLocation(Shader::Instance().GetCurrentProgram(),"slab1");
    slab_id_[1]=glGetUniformLocation(Shader::Instance().GetCurrentProgram(),"slab2");
    glUniform1f(glGetUniformLocation(Shader::Instance().GetCurrentProgram(),"level"),
                level);
    glUniform1f(glGetUniformLocation(Shader::Instance().GetCurrentProgram(),"tx"),
                1.0f/static_cast<float>(ucount));
    glUniform1f(glGetUniformLocation(Shader::Instance().GetCurrentProgram(),"ty"),
                1.0f/static_cast<float>(vcount));

    } // profile
    
#endif
    img::ImageStateBasePtr is=mh.ImageStatePtr();
    img::RealSpatialImageState* rsis=dynamic_cast<img::RealSpatialImageState*>(is.get());
    assert(rsis);
    set_slab(rsis,slab2,ucount,vcount,0);
    set_slab(rsis,slab1,ucount,vcount,1);
    BSlab bslab1 = process_edges(slab2,slab1,ucount,vcount,level,slab_tex_,slab_id_,0);
    process_vids(mh,va,bslab1,vidu1,vidv1,vidw1,ucount,vcount,0);
    set_slab(rsis,slab2,ucount,vcount,2);
    BSlab bslab2 = process_edges(slab1,slab2,ucount,vcount,level,slab_tex_,slab_id_,1);
    process_vids(mh,va,bslab2,vidu2,vidv2,vidw2,ucount,vcount,1);
    assign_lf(va,bslab1,vidu1,vidv1,vidw1,vidu2,vidv2,vidw2,ucount,vcount);
    
    for(uint w=2;w<wcount-1;++w) {
      slab1.swap(slab2);
      vidu1.swap(vidu2);
      vidv1.swap(vidv2);
      vidw1.swap(vidw2);
      bslab2.swap(bslab1);
      set_slab(rsis,slab2,ucount,vcount,w+1);
      bslab2 = process_edges(slab1,slab2,ucount,vcount,level,slab_tex_,slab_id_,w);
      process_vids(mh,va,bslab2,vidu2,vidv2,vidw2,ucount,vcount,w);
      assign_lf(va,bslab1,vidu1,vidv1,vidw1,vidu2,vidv2,vidw2,ucount,vcount);
    }
    
#if MAP_ISO_USE_SHADER
    Shader::Instance().PopProgram();
    glPopAttrib();
    OffscreenBuffer::Instance().End();
#endif
  } catch(...) {
#if MAP_ISO_USE_SHADER
    Shader::Instance().PopProgram();
    glPopAttrib();
    OffscreenBuffer::Instance().End();
#endif
    throw;
  }
  

}


}}} // ns
