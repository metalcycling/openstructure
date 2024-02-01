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

#include <boost/format.hpp>

#include "vertex_array.hh"

#include "glext_include.hh"
#include "color.hh"
#include "material.hh"
#include "scene.hh"
#include "vertex_array_helper.hh"
#include "povray.hh"
#include "exporter.hh"

#if OST_SHADER_SUPPORT_ENABLED
#include "shader.hh"
#endif


#define VERTEX_ARRAY_CHECK_GL_ERROR(m) \
  if((glerr=glGetError())!=0) {  \
    LOG_ERROR("Error during va buffer prep: " << m << " : " << gluErrorString(glerr)); \
    return false; \
  }

namespace ost { 

using namespace geom;

namespace gfx {

// the header and this file contain several magic number 7 uses related to the max buffer
static const int VA_VERTEX_BUFFER=0;
static const int VA_LINEINDEX_BUFFER=1;
static const int VA_TRIINDEX_BUFFER=2;
static const int VA_QUADINDEX_BUFFER=3;
static const int VA_AMBIENT_BUFFER=4;
static const int VA_NORMAL_BUFFER=5;
static const int VA_COLOR_BUFFER=6;



IndexedVertexArray::IndexedVertexArray()
{
  initialized_=false;
  Reset(); // replaces ctor initialization list
  //glGenTextures(1,&tex_id_);
}

IndexedVertexArray::~IndexedVertexArray()
{
}

IndexedVertexArray::IndexedVertexArray(const IndexedVertexArray& va)
{
  copy(va);
  //glGenTextures(1,&tex_id_);
}

IndexedVertexArray& IndexedVertexArray::operator=(const IndexedVertexArray& va)
{
  copy(va);
  // keep already allocated tex id
  return *this;
}

unsigned int IndexedVertexArray::GetFormat()
{
  // hardcoded for now, may be refactored and moved into an impl
  return GL_T2F_C4F_N3F_V3F;
}


void IndexedVertexArray::Cleanup() 
{
  if(initialized_) {
    //glDeleteTextures(1,&tex_id_);
    glDeleteLists(outline_mat_dlist_,1);
#if OST_SHADER_SUPPORT_ENABLED
    glDeleteBuffers(7,buffer_id_);
#endif
    initialized_=false;
  }
}

void IndexedVertexArray::SetMode(int m) {mode_=m;}
int IndexedVertexArray::GetMode() const {return mode_;}

void IndexedVertexArray::SetPolyMode(int m) {poly_mode_=m;}

void IndexedVertexArray::SetLighting(bool f) {lighting_=f;}
void IndexedVertexArray::SetTwoSided(bool f) {two_sided_=f;}
void IndexedVertexArray::SetCullFace(bool f) {cull_face_=f;}
void IndexedVertexArray::SetColorMaterial(bool f) {color_mat_=f;}
void IndexedVertexArray::SetLineWidth(float lw) {line_width_=lw;}
void IndexedVertexArray::SetAALines(bool f) {aalines_flag_=f;}
void IndexedVertexArray::SetPointSize(float ps) {point_size_=ps;}
void IndexedVertexArray::SetLineHalo(float lh) {line_halo_=lh;}
void IndexedVertexArray::SetOutlineMode(int m) {outline_mode_=m;}
void IndexedVertexArray::SetOutlineWidth(float f) {outline_width_=f;}
void IndexedVertexArray::SetOutlineMaterial(const Material& m) 
{
  outline_mat_=m;
  outline_mat_update_=true;
}
void IndexedVertexArray::SetOutlineExpandFactor(float f) { outline_exp_factor_=f;}
void IndexedVertexArray::SetOutlineExpandColor(const Color& c) {outline_exp_color_=c;}


unsigned int IndexedVertexArray::GetVertexCount() const
{
  return entry_list_.size();
}

void IndexedVertexArray::DumpVertices() const
{
  for(uint i=0;i<entry_list_.size();++i) {
    LOG_INFO("id=" << i << " v=" << entry_list_[i].v << " n=" << entry_list_[i].n << " c=" << entry_list_[i].c);
  }
}

LineID IndexedVertexArray::AddLine(VertexID id0, VertexID id1) 
{
  assert(id0<entry_list_.size() && id1<entry_list_.size());
  dirty_=true;
  line_index_list_.push_back(id0);
  line_index_list_.push_back(id1);
  return line_index_list_.size()-2;
}

TriID IndexedVertexArray::AddTri(VertexID id0, VertexID id1, VertexID id2) 
{
  assert(id0<entry_list_.size() && id1<entry_list_.size() && id2<entry_list_.size());
  dirty_=true;
  tri_index_list_.push_back(id0);
  tri_index_list_.push_back(id1);
  tri_index_list_.push_back(id2);
  return tri_index_list_.size()-3;
}

TriID IndexedVertexArray::AddTriN(VertexID id0, VertexID id1, VertexID id2)
{
  TriID tid = AddTri(id0,id1,id2);
  Vec3 d1= Vec3(entry_list_[id1].v)-Vec3(entry_list_[id0].v);
  Vec3 d2= Vec3(entry_list_[id2].v)-Vec3(entry_list_[id0].v);
  Vec3 x = Cross(d1,d2);
  NormalizerTriEntry nte = {
    tid, id0,id1,id2,
    Normalize(x),
#if 0
    1.0/(1.0+Length(x))
#else
    1.0
#endif
  };
  ntentry_list_.push_back(nte);
  return tid;
}

QuadID IndexedVertexArray::AddQuad(VertexID id0, VertexID id1, VertexID id2, VertexID id3) 
{
  assert(id0<entry_list_.size() && id1<entry_list_.size() && id2<entry_list_.size() && id3<entry_list_.size());
  dirty_=true;
  quad_index_list_.push_back(id0);
  quad_index_list_.push_back(id1);
  quad_index_list_.push_back(id2);
  quad_index_list_.push_back(id3);
  return quad_index_list_.size()-4;
}

void IndexedVertexArray::AddSphere(const SpherePrim& prim, unsigned int detail) 
{
  dirty_=true;

  unsigned int level= std::min(VA_SPHERE_MAX_DETAIL,detail);

  const detail::PrebuildSphereEntry &se = detail::GetPrebuildSphere(level);

  // use prebuild list of vertices which define given unit sphere
  std::vector<VertexID> vid_table;
  vid_table.reserve(se.vlist.size());
  for(std::vector<Vec3>::const_iterator it=se.vlist.begin();it!=se.vlist.end();++it) {
    VertexID id = Add(prim.radius*(*it)+prim.position,(*it),prim.color);
    vid_table.push_back(id);
  }
  for(unsigned int cc=0;cc<se.ilist.size();cc+=3) {
    AddTri(vid_table[se.ilist[cc+0]],vid_table[se.ilist[cc+1]],vid_table[se.ilist[cc+2]]);
  }
}

// add an icosahedral based sphere with the given params to the va
void IndexedVertexArray::AddIcoSphere(const SpherePrim& prim, unsigned int detail) 
{
  dirty_=true;
  
  unsigned int level= std::min(VA_ICO_SPHERE_MAX_DETAIL,detail);
  
  std::vector<Vec3> vlist = detail::GetPrebuildIcoSphere(level);
  
  for(std::vector<Vec3>::const_iterator it=vlist.begin();it!=vlist.end();it+=3) {
    VertexID id1 = Add(prim.radius*(*(it+0))+prim.position,(*(it+0)),prim.color);
    VertexID id2 = Add(prim.radius*(*(it+1))+prim.position,(*(it+1)),prim.color);
    VertexID id3 = Add(prim.radius*(*(it+2))+prim.position,(*(it+2)),prim.color);
    AddTri(id1,id2,id3);
  }
}

void IndexedVertexArray::AddCylinder(const CylinderPrim& prim, unsigned int detail, bool cap)
{
  dirty_=true;
  
  unsigned int level = std::min(VA_CYL_MAX_DETAIL,detail);
  
  const std::vector<Vec3>& vlist = detail::GetPrebuildCyl(level);
  
  Vec3 off(0.0,0.0,prim.length);

  Vec3 cn0 = cap ? prim.rotmat* geom::Vec3(0.0,0.0,-1.0) : Vec3();
  Vec3 cn1 = -cn0;
  VertexID cid0 = cap ? Add(prim.start, cn0 , prim.color1) : 0;
  VertexID cid7 = cap ? Add(prim.rotmat * off + prim.start, cn1, prim.color2) : 0;
  
  // prepare first vertices to add
  std::vector<Vec3>::const_iterator it=vlist.begin();
  Vec3 v0 = (*it);
  bool slant=(prim.radius1!=prim.radius2);
  // adjust for slant
  float beta = slant ? atan2(prim.radius1-prim.radius2,prim.length) : 0.0;
  float cosb = slant ? cos(beta) : 1.0;
  float sinb = slant ? sin(beta) : 0.0;
  Vec3 n0 = slant ? prim.rotmat * (cosb*v0+geom::Vec3(0.0,0.0,sinb)) : prim.rotmat*v0;

  v0*=prim.radius1;
  Vec3 v1 = (*it)*prim.radius2+off;
  VertexID id1 = Add(prim.rotmat * v0 + prim.start, n0, prim.color1);
  VertexID id2 = Add(prim.rotmat * v1 + prim.start, n0, prim.color2);
  VertexID cid1 = cap ? Add(prim.rotmat * v0 + prim.start, cn0, prim.color1) : 0;
  VertexID cid2 = cap ? Add(prim.rotmat * v1 + prim.start, cn1, prim.color2) : 0;
  
  // now for the loop around the circle
  VertexID id3=id1;
  VertexID id4=id2;
  VertexID cid3=cid1;
  VertexID cid4=cid2;
  ++it;
  for(;it!=vlist.end();++it) {
    v0 = (*it);
    Vec3 n0 = slant ? prim.rotmat * (cosb*v0+geom::Vec3(0.0,0.0,sinb)) : prim.rotmat*v0;
    v0 *= prim.radius1;
    Vec3 v1 = (*it)*prim.radius2+off;
    VertexID id5 = Add(prim.rotmat * v0 + prim.start, n0, prim.color1);
    VertexID id6 = Add(prim.rotmat * v1 + prim.start, n0, prim.color2);
    AddTri(id3,id5,id4);
    AddTri(id5,id6,id4);
    if(cap) {
      VertexID cid5 = Add(prim.rotmat * v0 + prim.start, cn0, prim.color1);
      VertexID cid6 = Add(prim.rotmat * v1 + prim.start, cn1, prim.color2);
      AddTri(cid0,cid5,cid3);
      AddTri(cid7,cid4,cid6);
      cid3=cid5;
      cid4=cid6;
    }
    id3=id5;
    id4=id6;
  }
  // and finally close the circle
  AddTri(id3,id1,id4);
  AddTri(id1,id2,id4);
  if(cap) {
    AddTri(cid0,cid1,cid3);
    AddTri(cid7,cid4,cid2);
  }
}

void IndexedVertexArray::SetOpacity(float o)
{
  o=std::max(0.0f,std::min(1.0f,o));
  if (std::abs(opacity_-o)<1e-6) {
    opacity_=o;
    return;
  }
  // this should really just set a value in the shader... 
  // but we need to support the fixed function pipeline as well
  for(EntryList::iterator it=entry_list_.begin();it!=entry_list_.end();++it) {
    it->c[3]=o;
  }
  opacity_=o;
  LOG_ERROR("new opacity (slow)" << opacity_);
  FlagRefresh();
}

void IndexedVertexArray::RenderGL() 
{
  static bool use_buff=false;
  
  if(!initialized_) {
    LOG_DEBUG("initializing vertex array lists");
#if OST_SHADER_SUPPORT_ENABLED
    glGenBuffers(7,buffer_id_);
#endif
    outline_mat_dlist_=glGenLists(1);
    initialized_=true;
  }

  if(dirty_) {
    dirty_=false;
#if OST_SHADER_SUPPORT_ENABLED

    LOG_DEBUG("checking buffer object availability");
    if(mode_&0x2 && aalines_flag_) {
      use_buff=false;
    } else {
      use_buff=prep_buff();
    }
    if(!use_buff) {
      LOG_TRACE("buffer not available");
      glBindBuffer(GL_ARRAY_BUFFER,0);
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
    } else {
      LOG_DEBUG("using buffer objects for vertex array");
    }
#endif
  }
  
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT);

  set_clip_offset(clip_offset_);

  if(use_tex_) {
    glEnable(GL_TEXTURE_2D);
  } else {
    glDisable(GL_TEXTURE_2D);
  }

  if(outline_mode_>0) {
    LOG_TRACE("outline rendering");
    if(outline_mat_update_) {
      glNewList(outline_mat_dlist_,GL_COMPILE);
      outline_mat_.RenderGL();
      glEndList();
      outline_mat_update_=false;
    }
    glDisable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
    glCullFace(GL_FRONT);
    glEnable(GL_CULL_FACE);
    glShadeModel(GL_FLAT);
    if(outline_mode_==1) {
      glCallList(outline_mat_dlist_);
      glEnable(GL_POLYGON_OFFSET_LINE);
      glEnable(GL_POLYGON_OFFSET_POINT);
      glPolygonOffset(10.0,1.0);
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
      glEnable(GL_LINE_SMOOTH);
      glDisable(GL_POINT_SMOOTH); // kills selfx fragment shader if enabled
      glDisable(GL_POLYGON_SMOOTH);
    } else if(outline_mode_==2) {
      glCallList(outline_mat_dlist_);
      glEnable(GL_POLYGON_OFFSET_LINE);
      glPolygonOffset(10.0,1.0);
      glDisable(GL_BLEND);
      glDisable(GL_LINE_SMOOTH);
      glDisable(GL_POINT_SMOOTH);
      glDisable(GL_POLYGON_SMOOTH);
      glCullFace(GL_FRONT);
    } else {
      glDisable(GL_POLYGON_OFFSET_LINE);
      glDisable(GL_POLYGON_OFFSET_POINT);
      glDisable(GL_LINE_SMOOTH);
      if(opacity_<1.0) {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
      } else {
        glDisable(GL_BLEND);
      }
    }
  } else {
    // not in outline mode
    if(lighting_) {
      glEnable(GL_LIGHTING); 
    } else {
      glDisable(GL_LIGHTING);
    }
    if(two_sided_) {
      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    } else {
      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
    }
    if(cull_face_ && !solid_ && clip_offset_<=0.0) {
      glEnable(GL_CULL_FACE);
    } else { 
      glDisable(GL_CULL_FACE); 
    }
  }
  if(color_mat_) {
    glEnable(GL_COLOR_MATERIAL); 
  } else {
    glDisable(GL_COLOR_MATERIAL);
  }
  
  if(mode_&0x1) {
    if(outline_mode_>0) {
      glPointSize(outline_width_);
    } else {
      glPointSize(point_size_);
    }
#if OST_SHADER_SUPPORT_ENABLED      
    Shader::Instance().UpdateState();
#endif      
    draw_p(use_buff);
  }
  if(mode_&0x6) { // 0x2 and 0x4
    if(outline_mode_==1) {
      // draw points beneath lines to reduce aliasing
      glLineWidth(outline_width_);
      glPointSize(outline_width_-1);
      glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
      draw_ltq(use_buff);
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      draw_ltq(use_buff);
    } else if(outline_mode_==2) {
      glLineWidth(outline_width_);
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      draw_ltq(use_buff);
    } else if(outline_mode_==3) {
#if OST_SHADER_SUPPORT_ENABLED
      Shader::Instance().PushProgram();
      Shader::Instance().Activate("outline");
      glUniform1f(glGetUniformLocation(Shader::Instance().GetCurrentProgram(),"scalef"),outline_exp_factor_);
      glUniform4f(glGetUniformLocation(Shader::Instance().GetCurrentProgram(),"color"),
                  outline_exp_color_[0],outline_exp_color_[1],
                  outline_exp_color_[2],opacity_);
      set_clip_offset(clip_offset_);
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      draw_ltq(use_buff);

      Shader::Instance().PopProgram();
#endif
    } else {
      glLineWidth(line_width_);
      if(poly_mode_==0) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
      } else if(poly_mode_==1) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      } else {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      }
#if OST_SHADER_SUPPORT_ENABLED
      Shader::Instance().UpdateState();
      if(aalines_flag_) {
        draw_aalines();
      } else {
#endif
        if((mode_&0x2 || (mode_&0x4 && poly_mode_==1)) && line_halo_>0.0) {
          draw_line_halo(use_buff);
        }
        draw_ltq(use_buff);

#if OST_SHADER_SUPPORT_ENABLED
      }
#endif
    }
  }

  if(draw_normals_) {
    //glColor3f(1,0,0);
    glBegin(GL_LINES);
    for(EntryList::const_iterator it=entry_list_.begin();it!=entry_list_.end();++it) {
      glNormal3fv(it->n);
      glVertex3fv(it->v);
      glNormal3fv(it->n);
      glVertex3f(it->v[0]+it->n[0]*0.75,
                 it->v[1]+it->n[1]*0.75,
                 it->v[2]+it->n[2]*0.75);
    }
    glEnd();
  }

  glPopClientAttrib();
  glPopAttrib();
}

namespace {
  unsigned int col_to_index(float* c)
  {
    // don't look too closely - I am lacking sufficient caffeine to do this more elegantly
    int ret= std::max(0,std::min<int>(511,static_cast<int>(round(c[0]*7.0f))*64+static_cast<int>(round(c[1]*7.0f))*8+static_cast<unsigned int>(round(c[2]*7.0f))));
    return static_cast<unsigned int>(ret);
  }
}

void IndexedVertexArray::RenderPov(PovState& pov, const std::string& name)
{
  if(entry_list_.empty()) return;

  pov.inc() << "mesh2 {\n";

  pov.inc() << " vertex_vectors { " << entry_list_.size() << ",\n";
  for(unsigned int i=0;i<entry_list_.size();++i) {
    pov.inc() << "  " << pov.write_coord(entry_list_[i].v);
    if(i+1<entry_list_.size()) pov.inc() << ",";
    pov.inc() << "\n";
  }
  pov.inc() << " }\n";

  pov.inc() << " normal_vectors { " << entry_list_.size() << ",\n";
  for(unsigned int i=0;i<entry_list_.size();++i) {
    pov.inc() << "  " << pov.write_norm(entry_list_[i].n);
    if(i+1<entry_list_.size()) pov.inc() << ",";
    pov.inc() << "\n";
  }
  pov.inc() << " }\n";

  pov.inc() << " texture_list { " << 8*8*8 << ",\n";
  for(unsigned int r=0;r<8;++r) {
    float fr=static_cast<float>(r)/7.0;
    for(unsigned int g=0;g<8;++g) {
      float fg=static_cast<float>(g)/7.0;
      for(unsigned int b=0;b<8;++b) {
        float fb=static_cast<float>(b)/7.0;
        pov.inc() << boost::format("  texture {_%s_tex pigment {color rgbft <%.4f,%.4f,%.4f,_%s_fi,_%s_tp>}}") % name % fr % fg % fb % name %name ;
        pov.inc() << ((r<7 && g<7 && b<7) ? ",\n" : "\n");
      }
    }
  }
  pov.inc() << " }\n";

  pov.inc() << " face_indices { " << tri_index_list_.size()/3 + quad_index_list_.size()/2;
  for(unsigned int c=0;c<tri_index_list_.size();c+=3) {
    pov.inc() << ",\n";
    unsigned int i0 = tri_index_list_[c];
    unsigned int i1 = tri_index_list_[c+1];
    unsigned int i2 = tri_index_list_[c+2];
    unsigned int x0 = col_to_index(entry_list_[i0].c);
    unsigned int x1 = col_to_index(entry_list_[i1].c);
    unsigned int x2 = col_to_index(entry_list_[i2].c);
    pov.inc() << boost::format("  <%u,%u,%u>, %u, %u, %u") % i0 % i1 % i2 % x0 % x1 % x2;
  }
  for(unsigned int c=0;c<quad_index_list_.size();c+=4) {
    pov.inc() << ",\n";
    unsigned int i0 = tri_index_list_[c];
    unsigned int i1 = tri_index_list_[c+1];
    unsigned int i2 = tri_index_list_[c+2];
    unsigned int i3 = tri_index_list_[c+3];
    unsigned int x0 = col_to_index(entry_list_[i0].c);
    unsigned int x1 = col_to_index(entry_list_[i1].c);
    unsigned int x2 = col_to_index(entry_list_[i2].c);
    unsigned int x3 = col_to_index(entry_list_[i3].c);
    pov.inc() << boost::format("  <%u,%u,%u>, %u, %u, %u") % i0 % i1 % i2 % x0 % x1 % x2;
    pov.inc() << boost::format("  <%u,%u,%u>, %u, %u, %u") % i1 % i2 % i3 % x1 % x2 % x3;
  }
  pov.inc() << "\n }\n";

  pov.inc() << "}\n";
}

void IndexedVertexArray::Export(Exporter* ex) const
{
  ex->WriteVertexData(entry_list_[0].v,entry_list_[0].n, entry_list_[0].c, entry_list_[0].t, sizeof(Entry), entry_list_.size());
  ex->WriteLineData(&line_index_list_[0],line_index_list_.size()/2);
  ex->WriteTriData(&tri_index_list_[0],tri_index_list_.size()/3);
  ex->WriteQuadData(&quad_index_list_[0],quad_index_list_.size()/4);
}

void IndexedVertexArray::Clear()
{
  dirty_=true;
  entry_list_.clear();
  quad_index_list_.clear();
  tri_index_list_.clear();
  line_index_list_.clear();
  ntentry_list_.clear();
} 

void IndexedVertexArray::Reset() 
{
  Clear();
  mode_=0x4;
  poly_mode_=2;
  lighting_=true;
  two_sided_=false;
  cull_face_=true;;
  color_mat_=true;
  line_width_=2.0;
  aalines_flag_=false;
  point_size_=2.0;
  line_halo_=0.0;
  opacity_=1.0;
  outline_mode_=0;
  outline_width_=4.0;
  outline_mat_=Material(0.3,1.0,0.0,0.0,0.0);
  outline_mat_update_=true;
  outline_exp_factor_=0.1;
  outline_exp_color_=Color(0,0,0);
  solid_=false;
  #undef RGB
  solid_color_=RGB(1,1,1);
  clip_offset_=0.0;
  draw_normals_=false;
  use_tex_=false;
}

void IndexedVertexArray::FlagRefresh()
{
  dirty_=true;
}

void IndexedVertexArray::CalcNormals(float smoothf) 
{
  if(ntentry_list_.empty()) return;

  std::vector<std::pair<bool,Vec3> > norm_list(entry_list_.size());

  // is this necessary or is the zeroing handled by std::vector ctor
  for(uint c=0;c<norm_list.size();++c) {
    norm_list[c].first=false;
  }

  for(NTEntryList::const_iterator it=ntentry_list_.begin();it!=ntentry_list_.end();++it) {
    norm_list[it->id0].first=true;
    norm_list[it->id1].first=true;
    norm_list[it->id2].first=true;
    // is this necessary or is the zeroing handled by std::vector ctor
    norm_list[it->id0].second=Vec3(0,0,0);
    norm_list[it->id1].second=Vec3(0,0,0);
    norm_list[it->id2].second=Vec3(0,0,0);
  }

  for(NTEntryList::const_iterator it=ntentry_list_.begin();it!=ntentry_list_.end();++it) {
    norm_list[it->id0].second+=it->weight*it->norm;
    norm_list[it->id1].second+=it->weight*it->norm;
    norm_list[it->id2].second+=it->weight*it->norm;
  }

  if(smoothf>0.0) {
    std::vector<std::pair<bool,Vec3> > norm_list2=norm_list;
    for(NTEntryList::const_iterator it=ntentry_list_.begin();it!=ntentry_list_.end();++it) {
      norm_list[it->id0].second+=smoothf*norm_list2[it->id1].second;
      norm_list[it->id0].second+=smoothf*norm_list2[it->id2].second;
      norm_list[it->id1].second+=smoothf*norm_list2[it->id0].second;
      norm_list[it->id1].second+=smoothf*norm_list2[it->id2].second;
      norm_list[it->id2].second+=smoothf*norm_list2[it->id0].second;
      norm_list[it->id2].second+=smoothf*norm_list2[it->id1].second;
    }
  }

  for(uint c=0;c<norm_list.size();++c) {
    if(norm_list[c].first) {
      Vec3 n=Normalize(norm_list[c].second);
      entry_list_[c].n[0]=n[0];
      entry_list_[c].n[1]=n[1];
      entry_list_[c].n[2]=n[2];
    }
  }

}

void IndexedVertexArray::CalcFullNormals() 
{
  LOG_DEBUG("calculating normals for vertex array");
  
  // book keeping setup
  static NormalizerTriEntry nte={0,0,0,0,Vec3(),1.0};
  
  NVEntryList nv_entry_list(entry_list_.size());
  unsigned int count=0;
  for(EntryList::const_iterator it = entry_list_.begin();
      it!=entry_list_.end(); ++it) {
    NormalizerVertexEntry nve;
    nve.pos=Vec3(it->v[0],it->v[1],it->v[2]);
    nve.tri_count=0;
    nve.weight=0.0;
    nv_entry_list[count++] = nve;
  }
  
  NTEntryList nt_entry_list(tri_index_list_.size()/3);
  count=0;
  for(IndexList::const_iterator it = tri_index_list_.begin();
      it!=tri_index_list_.end();) {
    nte.id0=*(it++);
    nte.id1=*(it++);
    nte.id2=*(it++);
    unsigned int nt_entry_id = count;
    nt_entry_list[count++]=nte;
    nv_entry_list[nte.id0].tri_list[nv_entry_list[nte.id0].tri_count++]=nt_entry_id;
    nv_entry_list[nte.id1].tri_list[nv_entry_list[nte.id1].tri_count++]=nt_entry_id;
    nv_entry_list[nte.id2].tri_list[nv_entry_list[nte.id2].tri_count++]=nt_entry_id;
  }
  
  // normals for all faces
  for(NTEntryList::iterator it=nt_entry_list.begin();
      it!=nt_entry_list.end();++it) {
    Vec3& p0=nv_entry_list[it->id0].pos;
    Vec3& p1=nv_entry_list[it->id1].pos;
    Vec3& p2=nv_entry_list[it->id2].pos;
    Vec3 d=Cross(Normalize(p1-p0),Normalize(p2-p0));
    //it->norm=Normalize(d);
    it->norm=d;
    it->weight=Dot(d,Cross(p0,p1)+Cross(p1,p2)+Cross(p2,p0));
    nv_entry_list[it->id0].weight+=acos(Dot(Normalize(p1-p0),Normalize(p2-p0)));
    nv_entry_list[it->id1].weight+=acos(Dot(Normalize(p0-p1),Normalize(p2-p1)));
    nv_entry_list[it->id2].weight+=acos(Dot(Normalize(p0-p2),Normalize(p1-p2)));
  }
  
  // now the normals for each vertex
  for(uint vc=0;vc<nv_entry_list.size();++vc) {
    NormalizerVertexEntry& nve=nv_entry_list[vc];
    if(nve.tri_count>0) {
      Vec3 norm(0.0,0.0,0.0);
      for(uint tc=0;tc<nve.tri_count;++tc) {
        TriID id = nve.tri_list[tc];
        if(id<nt_entry_list.size()) {
          NormalizerTriEntry& nte=nt_entry_list[id];
          //norm+=nte.norm;
          norm+=nte.norm*nte.weight;
          //norm+=nte.norm*nve.weight;
        } else {
          LOG_DEBUG("faulty vertex lookup in VA Normalize");
        }
      }
      norm=Normalize(norm);
      entry_list_[vc].n[0]=norm[0];
      entry_list_[vc].n[1]=norm[1];
      entry_list_[vc].n[2]=norm[2];
    }
  }
}

void IndexedVertexArray::DrawNormals(bool f)
{
  draw_normals_=f;
}


void IndexedVertexArray::SmoothNormals(float smoothf)
{
  if(smoothf==0.0) return;
  std::vector<std::pair<Vec3, int> > nlist(entry_list_.size());
#if 1
  // smooth by adjoining face normals
  for(uint c=0;c<tri_index_list_.size();) {
    VertexID id0 = tri_index_list_[c++];
    VertexID id1 = tri_index_list_[c++];
    VertexID id2 = tri_index_list_[c++];

    Vec3 n= Cross(Vec3(entry_list_[id2].v)-Vec3(entry_list_[id1].v),
                  Vec3(entry_list_[id0].v)-Vec3(entry_list_[id1].v));
    nlist[id0].first+=n;
    nlist[id0].second+=1;
    nlist[id1].first+=n;
    nlist[id1].second+=1;
    nlist[id2].first+=n;
    nlist[id2].second+=1;
  }
  for(uint c=0;c<quad_index_list_.size();) {
    VertexID id0 = quad_index_list_[c++];
    VertexID id1 = quad_index_list_[c++];
    VertexID id2 = quad_index_list_[c++];
    VertexID id3 = quad_index_list_[c++];
    Vec3 n= Cross(Vec3(entry_list_[id2].v)-Vec3(entry_list_[id1].v),
                  Vec3(entry_list_[id0].v)-Vec3(entry_list_[id1].v));
    nlist[id0].first+=n;
    nlist[id0].second+=1;
    nlist[id1].first+=n;
    nlist[id1].second+=1;
    nlist[id2].first+=n;
    nlist[id2].second+=1;
    nlist[id3].first+=n;
    nlist[id3].second+=1;
  }
#else
  // smooth by adjoining vertice normals
  for(uint c=0;c<tri_index_list_.size();) {
    VertexID id0 = tri_index_list_[c++];
    VertexID id1 = tri_index_list_[c++];
    VertexID id2 = tri_index_list_[c++];
    nlist[id0].first += Vec3(entry_list_[id1].n);
    nlist[id0].first += Vec3(entry_list_[id2].n);
    nlist[id0].second+=2;
    nlist[id1].first += Vec3(entry_list_[id2].n);
    nlist[id1].first += Vec3(entry_list_[id0].n);
    nlist[id1].second+=2;
    nlist[id2].first += Vec3(entry_list_[id0].n);
    nlist[id2].first += Vec3(entry_list_[id1].n);
    nlist[id2].second+=2;
  }
  for(uint c=0;c<quad_index_list_.size();) {
    VertexID id0 = quad_index_list_[c++];
    VertexID id1 = quad_index_list_[c++];
    VertexID id2 = quad_index_list_[c++];
    VertexID id3 = quad_index_list_[c++];
    nlist[id0].first += Vec3(entry_list_[id1].n);
    nlist[id0].first += Vec3(entry_list_[id2].n);
    nlist[id0].first += Vec3(entry_list_[id3].n);
    nlist[id0].second+=3;
    nlist[id1].first += Vec3(entry_list_[id2].n);
    nlist[id1].first += Vec3(entry_list_[id3].n);
    nlist[id1].first += Vec3(entry_list_[id0].n);
    nlist[id1].second+=3;
    nlist[id2].first += Vec3(entry_list_[id3].n);
    nlist[id2].first += Vec3(entry_list_[id0].n);
    nlist[id2].first += Vec3(entry_list_[id1].n);
    nlist[id2].second+=3;
    nlist[id3].first += Vec3(entry_list_[id0].n);
    nlist[id3].first += Vec3(entry_list_[id1].n);
    nlist[id3].first += Vec3(entry_list_[id2].n);
    nlist[id3].second+=3;
  }
#endif
  for(uint c=0;c<entry_list_.size();++c) {
    if(nlist[c].second>0) {
      float f=1.0/static_cast<float>(nlist[c].second);
      Vec3 n=Normalize((1.0-smoothf)*Vec3(entry_list_[c].n)+f*smoothf*nlist[c].first);
      entry_list_[c].n[0]=n[0];
      entry_list_[c].n[1]=n[1];
      entry_list_[c].n[2]=n[2];
    }
  }
}

void IndexedVertexArray::SmoothVertices(float smoothf)
{
  if(smoothf==0.0) return;
  std::vector<std::pair<Vec3, int> > dlist(entry_list_.size());

  for(uint c=0;c<tri_index_list_.size();) {
    VertexID id0 = tri_index_list_[c++];
    VertexID id1 = tri_index_list_[c++];
    VertexID id2 = tri_index_list_[c++];
    Vec3 p0(entry_list_[id0].v);
    Vec3 p1(entry_list_[id1].v);
    Vec3 p2(entry_list_[id2].v);
    dlist[id0].first+=p1-p0;
    dlist[id0].first+=p2-p0;
    dlist[id0].second+=2;
    dlist[id1].first+=p0-p1;
    dlist[id1].first+=p2-p1;
    dlist[id1].second+=2;
    dlist[id2].first+=p0-p2;
    dlist[id2].first+=p1-p2;
    dlist[id2].second+=2;
  }
  for(uint c=0;c<quad_index_list_.size();) {
    VertexID id0 = quad_index_list_[c++];
    VertexID id1 = quad_index_list_[c++];
    VertexID id2 = quad_index_list_[c++];
    VertexID id3 = quad_index_list_[c++];
    Vec3 p0(entry_list_[id0].v);
    Vec3 p1(entry_list_[id1].v);
    Vec3 p2(entry_list_[id2].v);
    Vec3 p3(entry_list_[id3].v);
    dlist[id0].first+=p1-p0;
    dlist[id0].first+=p3-p0;
    dlist[id0].second+=2;
    dlist[id1].first+=p0-p1;
    dlist[id1].first+=p2-p1;
    dlist[id1].second+=2;
    dlist[id2].first+=p3-p2;
    dlist[id2].first+=p1-p2;
    dlist[id2].second+=2;
    dlist[id3].first+=p2-p3;
    dlist[id3].first+=p0-p3;
    dlist[id3].second+=2;
  }
  for(uint c=0;c<entry_list_.size();++c) {
    Vec3 d = smoothf*dlist[c].first/static_cast<float>(dlist[c].second);
    entry_list_[c].v[0]+=d[0];
    entry_list_[c].v[1]+=d[1];
    entry_list_[c].v[2]+=d[2];
  }
}

namespace {

uint npatch_tab_id(uint u, uint v, uint N)
{
  return (N+1)*(N+2)/2-(N+1-v)*(N+2-v)/2+u;
}

}

void IndexedVertexArray::NPatch()
{
  uint N=3;

  EntryList entry_list;
  IndexList tri_index_list;
  IndexList line_index_list;

  entry_list.push_back(Entry(Vec3(),Vec3(),Color(1,0,0),Vec2()));

  for(uint c=0;c<tri_index_list_.size();) {
    VertexID id0 = tri_index_list_[c++];
    VertexID id1 = tri_index_list_[c++];
    VertexID id2 = tri_index_list_[c++];
    Vec3 p300(entry_list_[id0].v);
    Vec3 p030(entry_list_[id1].v);
    Vec3 p003(entry_list_[id2].v);
    Vec3 n200(entry_list_[id0].n);
    Vec3 n020(entry_list_[id1].n);
    Vec3 n002(entry_list_[id2].n);
    
    Vec3 p210 = (2.0*p300+p030-Dot(n200,p030-p300)*n200)/3.0;
    Vec3 p120 = (p300+2.0*p030-Dot(n020,p300-p030)*n020)/3.0;
    Vec3 p201 = (2.0*p300+p003-Dot(n200,p003-p300)*n200)/3.0;
    Vec3 p102 = (p003+2.0*p300-Dot(n002,p300-p003)*n002)/3.0;
    Vec3 p021 = (2.0*p030+p003-Dot(n020,p003-p030)*n020)/3.0;
    Vec3 p012 = (p030+2.0*p003-Dot(n002,p030-p003)*n002)/3.0;
    Vec3 p111 = (p210+p120+p102+p201+p021+p012)/4.0-(p300+p030+p003)/6.0;

    Vec3 d330 = p030-p300;
    Vec3 d033 = p003-p030;
    Vec3 d303 = p300-p003;
    Vec3 n110 = Normalize(n200+n020-2.0*Dot(d330,n200+n020)/Dot(d330,d330)*d330);
    Vec3 n011 = Normalize(n020+n002-2.0*Dot(d033,n020+n002)/Dot(d033,d033)*d033);
    Vec3 n101 = Normalize(n002+n200-2.0*Dot(d303,n002+n200)/Dot(d303,d303)*d303);

    std::vector<TriID> tab((N+1)*(N+2)/2);
    for(uint vc=0;vc<=N;++vc) {
      float v = static_cast<float>(vc)/static_cast<float>(N);
      for(uint uc=0;uc<=N-vc;++uc) {
        float u = static_cast<float>(uc)/static_cast<float>(N);
        float w = 1.0-u-v;
        Vec3 p = u*u*u*p300+v*v*v*p030+w*w*w*p003+
          3.0*(u*u*v*p210+u*u*w*p201+u*v*v*p120+v*v*w*p021+v*w*w*p012+u*w*w*p102)+
          6.0*u*v*w*p111;
        Vec3 n = Normalize(u*u*n200+v*v*n020+w*w*n002+u*v*n110+u*w*n101+v*w*n011);
        // TODO: patch texture coord as well
        entry_list.push_back(Entry(p,n,Color(),Vec2()));
        tab[npatch_tab_id(uc,vc,N)] = entry_list.size()-1;
      }
    }

    line_index_list.push_back(tab[npatch_tab_id(0,0,N)]);
    line_index_list.push_back(tab[npatch_tab_id(0,N,N)]);

    line_index_list.push_back(tab[npatch_tab_id(0,N,N)]);
    line_index_list.push_back(tab[npatch_tab_id(N,0,N)]);

    line_index_list.push_back(tab[npatch_tab_id(N,0,N)]);
    line_index_list.push_back(tab[npatch_tab_id(0,0,N)]);

    for(uint vc=0;vc<N;++vc) {
      for(uint uc=0;uc<N-vc;++uc) {
        tri_index_list.push_back(tab[npatch_tab_id(uc,vc,N)]);
        tri_index_list.push_back(tab[npatch_tab_id(uc+1,vc,N)]);
        tri_index_list.push_back(tab[npatch_tab_id(uc,vc+1,N)]);
        if(uc+vc<N-1) {
          tri_index_list.push_back(tab[npatch_tab_id(uc+1,vc,N)]);
          tri_index_list.push_back(tab[npatch_tab_id(uc+1,vc+1,N)]);
          tri_index_list.push_back(tab[npatch_tab_id(uc,vc+1,N)]);
        }
      }
    }
  }

  entry_list_.swap(entry_list);
  line_index_list_.swap(line_index_list);
  tri_index_list_.swap(tri_index_list);
  quad_index_list_.clear();

  FlagRefresh();
}


////////////////////////////////////////
// private methods

    
void IndexedVertexArray::copy(const IndexedVertexArray& va)
{
  initialized_=false;
  entry_list_=va.entry_list_;
  quad_index_list_=va.quad_index_list_;
  tri_index_list_=va.tri_index_list_;
  line_index_list_=va.line_index_list_;
  ntentry_list_=va.ntentry_list_;
  dirty_=true;
  mode_=va.mode_;
  poly_mode_=va.poly_mode_;
  lighting_=va.lighting_;
  two_sided_=va.two_sided_;
  cull_face_=va.cull_face_;
  color_mat_=va.color_mat_;
  line_width_=va.line_width_;
  aalines_flag_=va.aalines_flag_;
  point_size_=va.point_size_;
  line_halo_=va.line_halo_;
  outline_mode_=va.outline_mode_;
  outline_width_=va.outline_width_;
  outline_mat_=va.outline_mat_;
  outline_mat_update_=true;
  outline_exp_factor_=va.outline_exp_factor_;
  outline_exp_color_=va.outline_exp_color_;
  solid_=va.solid_;
  solid_color_=va.solid_color_;
  clip_offset_=va.clip_offset_;
  draw_normals_=va.draw_normals_;
  use_tex_=va.use_tex_;
}
  
bool IndexedVertexArray::prep_buff()
{
#if OST_SHADER_SUPPORT_ENABLED
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
  glEnableClientState(GL_INDEX_ARRAY);

  int glerr=glGetError(); // clear error flag
  glBindBuffer(GL_ARRAY_BUFFER, buffer_id_[VA_VERTEX_BUFFER]);
  VERTEX_ARRAY_CHECK_GL_ERROR("bind vertex buf");
  glBufferData(GL_ARRAY_BUFFER, 
               sizeof(Entry) * entry_list_.size(),
               &entry_list_[0],
               GL_STATIC_DRAW);
  VERTEX_ARRAY_CHECK_GL_ERROR("set vertex buf");

  // this will be used later, when refactoring v,c,n to live in separate lists
  /*
  glBindBuffer(GL_ARRAY_BUFFER, buffer_id_[VA_COLOR_BUFFER]);
  VERTEX_ARRAY_CHECK_GL_ERROR("bind buf1");
  glBufferData(GL_ARRAY_BUFFER, 
               sizeof(Entry) * entry_list_.size(),
               &entry_list_[0].c[0],
               GL_STATIC_DRAW);
  VERTEX_ARRAY_CHECK_GL_ERROR("set buf1");

  glBindBuffer(GL_ARRAY_BUFFER, buffer_id_[VA_NORMAL_BUFFER]);
  VERTEX_ARRAY_CHECK_GL_ERROR("bind buf2");
  glBufferData(GL_ARRAY_BUFFER, 
               sizeof(Entry) * entry_list_.size(),
               &entry_list_[0].n[0],
               GL_STATIC_DRAW);
  VERTEX_ARRAY_CHECK_GL_ERROR("set buf2");
  */

  if(!line_index_list_.empty()) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer_id_[VA_LINEINDEX_BUFFER]);
    VERTEX_ARRAY_CHECK_GL_ERROR("bind lindex buf");
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 
                 sizeof(unsigned int) * line_index_list_.size(),
                 &line_index_list_[0],
                 GL_STATIC_DRAW);
    VERTEX_ARRAY_CHECK_GL_ERROR("set lindex buf");
  }

  if(!tri_index_list_.empty()) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer_id_[VA_TRIINDEX_BUFFER]);
    VERTEX_ARRAY_CHECK_GL_ERROR("bind tindex buf");
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 
                 sizeof(unsigned int) * tri_index_list_.size(),
                 &tri_index_list_[0],
                 GL_STATIC_DRAW);
    VERTEX_ARRAY_CHECK_GL_ERROR("set tindex buf");
  }

  if(!quad_index_list_.empty()) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer_id_[VA_QUADINDEX_BUFFER]);
    VERTEX_ARRAY_CHECK_GL_ERROR("bind qindex buf");
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 
                 sizeof(unsigned int) * quad_index_list_.size(),
                 &quad_index_list_[0],
                 GL_STATIC_DRAW);
    VERTEX_ARRAY_CHECK_GL_ERROR("set qindex buf");
  }

  return true;
#else
  return false;
#endif
}

void IndexedVertexArray::draw_ltq(bool use_buff)
{
  if(solid_) {
    glClearStencil(0x0);
    glClear(GL_STENCIL_BUFFER_BIT);
    glEnable(GL_STENCIL_TEST);
    glStencilFunc(GL_ALWAYS,0x0,0x1);
    glStencilOp(GL_KEEP,GL_INVERT,GL_INVERT);
  }

  if(use_buff) {
#if OST_SHADER_SUPPORT_ENABLED
#if 0
    /*
      for now, since v,n,c live in a packed format (formerly used
      with glInterleavedArrays), only a single buffer is
      used, with the gl*Pointer calls giving the byte offset
      in place of the absolute pointer
    */
    glBindBuffer(GL_ARRAY_BUFFER, buffer_id_[VA_VERTEX_BUFFER]);
    glVertexPointer(3, GL_FLOAT, sizeof(Entry), reinterpret_cast<void*>(sizeof(float)*9));
    glNormalPointer(GL_FLOAT, sizeof(Entry), reinterpret_cast<void*>(sizeof(float)*6));
    glColorPointer(4, GL_FLOAT, sizeof(Entry), reinterpret_cast<void*>(sizeof(float)*2));
    if(use_tex_) {
      glTexCoordPointer(2, GL_FLOAT, sizeof(Entry), 0);
    }
#else
    glBindBuffer(GL_ARRAY_BUFFER, buffer_id_[VA_VERTEX_BUFFER]);
    glInterleavedArrays(GetFormat(),sizeof(Entry),NULL);
#endif

    glBindBuffer(GL_ARRAY_BUFFER,0);

    if(!tri_index_list_.empty() && (mode_ & 0x4)) {
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer_id_[VA_TRIINDEX_BUFFER]);
      glDrawElements(GL_TRIANGLES,tri_index_list_.size(),GL_UNSIGNED_INT,NULL);
    }
    if(!quad_index_list_.empty() && (mode_ & 0x4)) {
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer_id_[VA_QUADINDEX_BUFFER]);
      glDrawElements(GL_QUADS,quad_index_list_.size(),GL_UNSIGNED_INT,NULL);
    }
    if(!line_index_list_.empty() && (mode_ & 0x2)) {
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer_id_[VA_LINEINDEX_BUFFER]);
      glDrawElements(GL_LINES,line_index_list_.size(),GL_UNSIGNED_INT,NULL);
    }
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
#endif
  } else {
    
    LOG_TRACE("setting up vertex array");
    glInterleavedArrays(GetFormat(),sizeof(Entry),&entry_list_[0]);
    
    if(!tri_index_list_.empty() && (mode_ & 0x4)) {
      LOG_TRACE("rendering vertex arras tris");
      glDrawElements(GL_TRIANGLES,tri_index_list_.size(),GL_UNSIGNED_INT,&tri_index_list_[0]);
    }
    if(!quad_index_list_.empty() && (mode_ & 0x4)) {
      LOG_TRACE("rendering vertex arras quads");
      glDrawElements(GL_QUADS,quad_index_list_.size(),GL_UNSIGNED_INT,&quad_index_list_[0]);
    }
    if(!line_index_list_.empty() && (mode_ & 0x2)) {
      LOG_TRACE("rendering vertex arras lines");
      glDrawElements(GL_LINES,line_index_list_.size(),GL_UNSIGNED_INT,&line_index_list_[0]);
    }
  }

  if(solid_) {
    glStencilFunc(GL_EQUAL,0x1,0x1);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    float fov=Scene::Instance().GetFOV()*M_PI/180.0;
    float znear=Scene::Instance().GetNear();
    float aspect=Scene::Instance().GetAspect();
    float rh=2.0*fabs(tan(fov)*znear);
    float rw=rh*aspect;
    float rz=-(znear+clip_offset_)-0.05;

    glDisable(GL_LIGHTING);

    set_clip_offset(0.0);

    glBegin(GL_TRIANGLE_STRIP);
    glColor3fv(solid_color_);
    glNormal3f(0,0,1);
    glVertex3f(-rw,-rh,rz);
    glVertex3f(rw,-rh,rz);
    glVertex3f(-rw,rh,rz);
    glVertex3f(rw,rh,rz);
    glEnd();

    set_clip_offset(clip_offset_);

    glDisable(GL_STENCIL_TEST);
    glPopMatrix();
  }
}

void IndexedVertexArray::draw_p(bool use_buff)
{
  if(use_buff) {
#if OST_SHADER_SUPPORT_ENABLED
    glBindBuffer(GL_ARRAY_BUFFER, buffer_id_[VA_VERTEX_BUFFER]);
    glVertexPointer(3, GL_FLOAT, sizeof(Entry), reinterpret_cast<void*>(sizeof(float)*9));
    glNormalPointer(GL_FLOAT, sizeof(Entry), reinterpret_cast<void*>(sizeof(float)*6));
    glColorPointer(4, GL_FLOAT, sizeof(Entry), reinterpret_cast<void*>(sizeof(float)*2));
    if(use_tex_) {
      glTexCoordPointer(2, GL_FLOAT, sizeof(Entry), 0);
    }
    glDrawArrays(GL_POINTS,0,entry_list_.size());
#endif
  } else {
    LOG_TRACE("calling vertex array");
    glInterleavedArrays(GetFormat(),sizeof(Entry),&entry_list_[0]);
    glDrawArrays(GL_POINTS,0,entry_list_.size());
  }
}

namespace {

  struct AALineEntry {
    float p0[3], p1[3];
    float edge0[3],edge1[3],edge2[3],edge3[3];
    float c0[4], c1[4];
    float z;
  };

  typedef std::vector<AALineEntry> AALineList;

  struct AALineEntryLess
  {
    bool operator()(const AALineEntry& e1, const AALineEntry& e2)
    {
      return e1.z<e2.z;
    }
  };

  /*
     Adapted after Chan and Durand, "Fast Prefiltered Lines",
     in GPU Gems 2
     */
  void render_aalines(AALineList& line_list, float w, float r)
  {
#if OST_SHADER_SUPPORT_ENABLED
    static GLuint table_tex_id=0;
    std::vector<unsigned char> tex_table(32);

    if(table_tex_id==0) {
      for(int i=0;i<32;++i) {
        float x=static_cast<float>(i)/31.0;
        tex_table[31-i]=static_cast<unsigned char>(255.0*exp(-4.0*x*x));
      }
      glGenTextures(1,&table_tex_id);
      glBindTexture(GL_TEXTURE_1D,table_tex_id);
      glPixelStorei(GL_UNPACK_ALIGNMENT,1);
      glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexImage1D(GL_TEXTURE_1D,0,GL_LUMINANCE,32,0,GL_LUMINANCE,GL_UNSIGNED_BYTE,&tex_table[0]);
    }
    std::sort(line_list.begin(),line_list.end(),AALineEntryLess());

    Shader::Instance().PushProgram();
    Shader::Instance().Activate("aaline");

    glLineWidth(ceil((2.0f*r+w)*1.5));

    glPushAttrib(GL_ENABLE_BIT | GL_DEPTH_BUFFER_BIT);
    glDisable(GL_NORMALIZE);
    glDisable(GL_CULL_FACE);
    glDisable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
    //glDisable(GL_DEPTH_TEST);
    glDisable(GL_LINE_SMOOTH);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    glDepthMask(0);
    Shader::Instance().UpdateState();

    GLuint cpr = Shader::Instance().GetCurrentProgram();
    GLint edge0_id = glGetUniformLocation(cpr,"edge0");
    GLint edge1_id = glGetUniformLocation(cpr,"edge1");
    GLint edge2_id = glGetUniformLocation(cpr,"edge2");
    GLint edge3_id = glGetUniformLocation(cpr,"edge3");
    GLint table_id = glGetUniformLocation(cpr,"table");

    glActiveTexture(GL_TEXTURE0);
    glDisable(GL_TEXTURE_2D);
    glEnable(GL_TEXTURE_1D);
    glBindTexture(GL_TEXTURE_1D,table_tex_id);
    glUniform1i(table_id,0);

    for(AALineList::iterator it=line_list.begin();
        it!=line_list.end();++it) {
      glUniform3fv(edge0_id,1,it->edge0);
      glUniform3fv(edge1_id,1,it->edge1);
      glUniform3fv(edge2_id,1,it->edge2);
      glUniform3fv(edge3_id,1,it->edge3);
      // glUniform cannot be inside a glBegin / glEnd
      glBegin(GL_LINES);
      glColor4fv(it->c0);
      glVertex3fv(it->p0);
      glColor4fv(it->c1);
      glVertex3fv(it->p1); 
      glEnd();
    }

    glPopAttrib();

    Shader::Instance().PopProgram();
#endif
  }

  Vec3 make_aaline_edge(const Vec2& c1, const Vec2& c0, float s)
  {
    Vec3 nrvo(c1[1]-c0[1],c0[0]-c1[0],c1[0]*c0[1]-c0[0]*c1[1]);
    nrvo*=1.0/(s*Length(c1-c0));
    return nrvo;
  }

}

void IndexedVertexArray::draw_aalines()
{
#if OST_SHADER_SUPPORT_ENABLED
  float w=0.5*line_width_;
  float r=w;
  float hwr = 0.5*w+r;

  AALineList line_list;
  line_list.reserve(line_index_list_.size() >> 1);
  for(unsigned int i=0;i<line_index_list_.size();i+=2) {
    Entry& ve0 = entry_list_[line_index_list_[i]];
    Entry& ve1 = entry_list_[line_index_list_[i+1]];
    Vec3 q0(Scene::Instance().Project(Vec3(ve0.v),false));
    Vec3 q1(Scene::Instance().Project(Vec3(ve1.v),false));
    Vec2 p0(q0[0],q0[1]);
    Vec2 p1(q1[0],q1[1]);
    Vec2 d=Normalize(p1-p0);
    Vec2 n(d[1],-d[0]);
    Vec2 c0=p0-(r)*d+(hwr)*n;
    Vec2 c1=p1+(r)*d+(hwr)*n;
    Vec2 c2=p1+(r)*d-(hwr)*n;
    Vec2 c3=p0-(r)*d-(hwr)*n;
    Vec3 e0=make_aaline_edge(c0,c1,r);
    Vec3 e1=make_aaline_edge(c1,c2,r);
    Vec3 e2=make_aaline_edge(c2,c3,r);
    Vec3 e3=make_aaline_edge(c3,c0,r);

    // TODO: it should be possible to avoid the unproject and
    // pass the projected coordinates directly to GL, setting
    // modelview and projection to unity; requires handling of
    // the viewport though!
    Vec3 np1=Scene::Instance().UnProject(q0-Vec3(hwr*d),false);
    Vec3 np2=Scene::Instance().UnProject(q1+Vec3(hwr*d),false);

    AALineEntry le={{np1[0],np1[1],np1[2]},
      {np2[0],np2[1],np2[2]},
      {e0[0],e0[1],e0[2]},
      {e1[0],e1[1],e1[2]},
      {e2[0],e2[1],e2[2]},
      {e3[0],e3[1],e3[2]},
      {ve0.c[0],ve0.c[1],ve0.c[2],ve0.c[3]},
      {ve1.c[0],ve1.c[1],ve1.c[2],ve1.c[3]},
      -0.5f*(q0[2]+q1[2])};
    line_list.push_back(le);
  }

  render_aalines(line_list,w,r);
#endif
}

void IndexedVertexArray::draw_line_halo(bool use_buff)
{
  glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_FALSE);
  glLineWidth(line_width_+line_halo_);
  glPushAttrib(GL_ENABLE_BIT);
  glDisable(GL_LIGHTING);
  glDisable(GL_COLOR_MATERIAL);
  glEnable(GL_LINE_SMOOTH);
  glDisable(GL_POINT_SMOOTH);
  glDisable(GL_POLYGON_SMOOTH);
  glEnable(GL_BLEND);
  glDisable(GL_FOG);
  float epsilon=0.05;
  glDepthRange(0.0+epsilon,1.0+epsilon);

  draw_ltq(use_buff);

  glDepthRange(0.0,1.0);
  glPopAttrib();
  glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
  glLineWidth(line_width_);
}

geom::AlignedCuboid IndexedVertexArray::GetBoundingBox() const
{
  if(entry_list_.empty()) {
    return geom::AlignedCuboid(geom::Vec3(0,0,0),geom::Vec3(0,0,0));
  } else {
    geom::Vec3 minc(std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max());
    geom::Vec3 maxc(-std::numeric_limits<float>::max(),
        -std::numeric_limits<float>::max(),
        -std::numeric_limits<float>::max());
    for(EntryList::const_iterator it=entry_list_.begin();it!=entry_list_.end();++it) {
      geom::Vec3 p(it->v[0],it->v[1],it->v[2]);
      minc=geom::Min(minc,p);
      maxc=geom::Max(maxc,p);
    }
    return geom::AlignedCuboid(minc-1.0,maxc+1.0);
  }
}

void IndexedVertexArray::set_clip_offset(float o)
{
#if OST_SHADER_SUPPORT_ENABLED
  float n=Scene::Instance().GetNear();
  float f=Scene::Instance().GetFar();
  float t=(f*(-n+o)+n*(n+o))/((f-n)*(n+o));
  t=(t+1.0)*0.5;
  glUniform1f(glGetUniformLocation(Shader::Instance().GetCurrentProgram(),"clip_offset"),t);
#endif
}

}} // ns
