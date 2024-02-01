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
#include <ost/log.hh>
#include <ost/dyn_cast.hh>

#include <ost/mol/entity_property_mapper.hh>

#include "gl_helper.hh"
#include "glext_include.hh"
#include "scene.hh"
#include "gfx_object.hh"

#include "povray.hh"
#include "impl/mapped_property.hh"
#include "exporter.hh"
#include <ost/img/alg/stat.hh>

#if OST_SHADER_SUPPORT_ENABLED
#include "shader.hh"
#endif

namespace ost { namespace gfx {

GfxObj::GfxObj(const String& name):
  GfxObjBase(name),
  va_(),
  render_mode_(RenderMode::SIMPLE),
  debug_flags_(0),
  transform_(),
  rebuild_(true),
  refresh_(false),
  line_width_(2.0),
  poly_mode_(2),
  aalines_flag_(false),
  line_halo_(0.0),
  mat_(0.0,1.0,0.5,64.0,0.0),
  mat_dlist_(0),
  mat_update_(true),
  opacity_(1.0),
  smoothf_(0.0),
  outline_flag_(false),
#if defined (__APPLE__)
  outline_mode_(2),
#else
  outline_mode_(1),
#endif
  solid_(false),
  solid_color_(RGB(0.7,0.7,0.7)),
  clip_flag_(false),
  clip_plane_(),
  clip_offset_(0.0),
  c_ops_(),
  labels_(),
  use_occlusion_(false)
{
}

GfxObj::GfxObj(const GfxObj&):
  GfxObjBase("") // to make the compiler happy
{}

GfxObj& GfxObj::operator=(const GfxObj&) {return *this;}

////////////////////////////////////////
// the GfxNode interface

GfxNodeP GfxObj::Copy() const
{
  return GfxNodeP(new GfxObj(*this));
}

void GfxObj::DeepSwap(GfxObj& go)
{
  GfxNode::DeepSwap(*this);
  std::swap(transform_,go.transform_);
  std::swap(rebuild_,go.rebuild_);
  std::swap(refresh_,go.refresh_);
  std::swap(poly_mode_,go.poly_mode_);
  std::swap(aalines_flag_,go.aalines_flag_);
  std::swap(line_halo_,go.line_halo_);
  std::swap(mat_,go.mat_);
  std::swap(mat_dlist_,go.mat_dlist_);
  std::swap(mat_update_,go.mat_update_);
  std::swap(opacity_,go.opacity_);
  std::swap(smoothf_,go.smoothf_);
  std::swap(outline_flag_,go.outline_flag_);
  std::swap(outline_mode_,go.outline_mode_);
  std::swap(solid_,go.solid_);
  std::swap(solid_color_,go.solid_color_);
  std::swap(clip_flag_,go.clip_flag_);
  std::swap(clip_plane_,go.clip_plane_);
  std::swap(clip_offset_,go.clip_offset_);
  std::swap(c_ops_,go.c_ops_);
  std::swap(labels_,go.labels_);
  std::swap(use_occlusion_,go.use_occlusion_);
}

void GfxObj::RenderGL(RenderPass pass)
{
  LOG_TRACE("object " << GetName() << ": RenderGL()");

  if(pass==0) {
    if(mat_update_) {
      LOG_TRACE("updating material display list");
      if(mat_dlist_==0) {
        mat_dlist_=glGenLists(1);
      }
      glNewList(mat_dlist_,GL_COMPILE);
      mat_.RenderGL();
      glEndList();
      mat_update_=false;
    }
    
    if(rebuild_ || refresh_) {
      PreRenderGL(rebuild_);
      rebuild_=false;
      refresh_=false;
      // there really needs to be a central place
      // where the va attributes are re-applied
      va_.SetOpacity(opacity_);
      va_.SetSolid(solid_);
      va_.SetSolidColor(solid_color_);
    }
  }
  if(IsVisible()) {
    LOG_TRACE("applying local transformation");
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMultMatrix(transform_.GetTransposedMatrix().Data());
    LOG_TRACE("applying material display list");
    glCallList(mat_dlist_);
    LOG_TRACE("calling custom render gl pass " << pass);

    /*
      only STANDARD_RENDER_PASS and GLOW_RENDER_PASS are
      passed down to the custom rendering routines
    */
  
    if(clip_flag_) {
#if OST_SHADER_SUPPORT_ENABLED && GL_VERSION_3_0
      glEnable(GL_CLIP_DISTANCE0);
      GLuint cp = Shader::Instance().GetCurrentProgram();
      if(cp>0) {
        glUniform1i(glGetUniformLocation(cp,"clip_flag"),1);
        glUniform4f(glGetUniformLocation(cp,"clip_plane"),
                   clip_plane_[0],clip_plane_[1],clip_plane_[2],clip_plane_[3]);
      } else {
        GLdouble eq[4];
        for(size_t i=0;i<4;++i) eq[i]=clip_plane_[i];
        glClipPlane(GL_CLIP_PLANE0,eq);
      }
#else
      GLdouble eq[4];
      for(size_t i=0;i<4;++i) eq[i]=clip_plane_[i];
      glClipPlane(GL_CLIP_PLANE0,eq);
#endif
    }

    if(pass==DEPTH_RENDER_PASS) {
      render_depth_only();
    } else if(pass==TRANSPARENT_RENDER_PASS) {
      if(GetOpacity()<1.0) {
        render_depth_only();
        CustomRenderGL(STANDARD_RENDER_PASS);
        if(outline_flag_) {
          va_.SetOutlineMode(outline_mode_);
          CustomRenderGL(pass);
          va_.SetOutlineMode(0);
        }
      }
    } else if(pass==STANDARD_RENDER_PASS) {
      if(GetOpacity()>=1.0) {
        CustomRenderGL(STANDARD_RENDER_PASS);
        if(outline_flag_) {
          va_.SetOutlineMode(outline_mode_);
          CustomRenderGL(pass);
          va_.SetOutlineMode(0);
        }
      }
    } else if(pass==GLOW_RENDER_PASS) {
      CustomRenderGL(GLOW_RENDER_PASS);
    } else if(pass==OVERLAY_RENDER_PASS) {
      LOG_TRACE("drawing labels");
      render_labels();
    }

    if(clip_flag_) {
#if OST_SHADER_SUPPORT_ENABLED && GL_VERSION_3_0
      glDisable(GL_CLIP_DISTANCE0);
#endif
    }

    glPopMatrix();    
  }
}

void GfxObj::InitGL()
{
}

void GfxObj::RenderPov(PovState& pov)
{
  if(IsVisible()) {
    pov.start_obj(GetName(),1.0,1.0,1.0);
    // apply local transformation
    // using transform_
    if(rebuild_ || refresh_) {
      PreRenderGL(true);
    }
    CustomRenderPov(pov);
    pov.end_obj();
  }
}


void GfxObj::Export(Exporter* ex)
{
  if(IsVisible()) {
    ex->NodeStart(GetName(),Exporter::OBJ);
    // in the simplest case, just export va
    if(rebuild_ || refresh_) {
      PreRenderGL(true);
    }
    va_.Export(ex);
    ex->NodeEnd(GetName());
  }
}

void GfxObj::Apply(GfxNodeVisitor& v, GfxNodeVisitor::Stack st)
{
  v.VisitObject(this,st);
}

int GfxObj::GetType() const
{
  return 1;
}

////////////////////////////////////////
// now for the GfxObjBase interface
void GfxObj::SetMatAmb(const Color& c)
{
  mat_.SetAmb(c);
  mat_update_=true;
  Scene::Instance().RequestRedraw();
}

void GfxObj::SetMatDiff(const Color& c)
{
  mat_.SetDiff(c);
  mat_update_=true;
  Scene::Instance().RequestRedraw();
}

void GfxObj::SetMatSpec(const Color& c)
{
  mat_.SetSpec(c);
  mat_update_=true;
  Scene::Instance().RequestRedraw();
}

void GfxObj::SetMatShin(float s)
{
  mat_.SetShin(s);
  mat_update_=true;
  Scene::Instance().RequestRedraw();
}

void GfxObj::SetMatEmm(const Color& c)
{
  mat_.SetEmm(c);
  mat_update_=true;
  Scene::Instance().RequestRedraw();
}

void GfxObj::ContextSwitch()
{
  FlagRebuild();
  GfxNode::ContextSwitch();
}

void GfxObj::SetRenderMode(RenderMode::Type m)
{
  if (render_mode_==m) return;
  render_mode_=m;
  OnRenderModeChange();
  FlagRebuild();
}

RenderMode::Type GfxObj::GetRenderMode() const
{
  return render_mode_;
}

geom::Vec3 GfxObj::GetCenter() const 
{
  return this->GetBoundingBox().GetCenter();
}

void GfxObj::SetLineWidth(float w)
{
  va_.SetLineWidth(w);
  line_width_=std::max(float(0.01),w);
  FlagRefresh();
  Scene::Instance().RenderModeChanged(GetName());
}

void GfxObj::SetPolyMode(unsigned int m)
{
  if(m==poly_mode_) return;
  poly_mode_=std::min((unsigned int)2,m);
  va_.SetPolyMode(poly_mode_);
  FlagRefresh();
}

void GfxObj::SetAALines(bool f)
{
  if(f==aalines_flag_) return;
  va_.SetAALines(f);
  aalines_flag_=f;
  FlagRefresh();
}

void GfxObj::SetLineHalo(float f)
{
  va_.SetLineHalo(f);
  line_halo_=f;
  FlagRefresh();
}

void GfxObj::SetOutline(bool f)
{
  outline_flag_=f;
  FlagRefresh();
  if(f) {
    outline_mode_=std::min(3,std::max(1,outline_mode_));
  }
  Scene::Instance().RequestRedraw();
}

void GfxObj::SetOutlineMode(int m)
{
  outline_mode_=m;
  if(outline_flag_) {
    FlagRefresh();
    Scene::Instance().RequestRedraw();
  }
}

void GfxObj::SetOutlineWidth(float f)
{
  va_.SetOutlineWidth(f);
  Scene::Instance().RequestRedraw();
}

float GfxObj::GetOutlineWidth() const
{
  return va_.GetOutlineWidth();
}

void GfxObj::SetOutlineExpandFactor(float f)
{
  va_.SetOutlineExpandFactor(f);
  Scene::Instance().RequestRedraw();
}

float GfxObj::GetOutlineExpandFactor() const
{
  return va_.GetOutlineExpandFactor();
}

void GfxObj::SetOutlineExpandColor(const Color& c)
{
  va_.SetOutlineExpandColor(c);
  Scene::Instance().RequestRedraw();
}

Color GfxObj::GetOutlineExpandColor() const
{
  return va_.GetOutlineExpandColor();
}

void GfxObj::SetOpacity(float o)
{
  opacity_=o;
  va_.SetOpacity(opacity_);
  FlagRefresh();
  Scene::Instance().RequestRedraw();
}

void GfxObj::SetSolid(bool f)
{
  solid_=f;
  va_.SetSolid(solid_);
  Scene::Instance().RequestRedraw();
}

void GfxObj::SetSolidColor(const Color& c)
{
  solid_color_=c;
  va_.SetSolidColor(solid_color_);
  Scene::Instance().RequestRedraw();
}

void GfxObj::SetClip(bool f)
{
  clip_flag_=f;
  // va_.SetClip(clip_flag_);
  Scene::Instance().RequestRedraw();
}

void GfxObj::SetClipPlane(const geom::Vec4& p)
{
  clip_plane_=p;
  // va_.SetSolidColor(solid_color_);
  Scene::Instance().RequestRedraw();
}

void GfxObj::SetClipOffset(float f)
{
  clip_offset_=f;
  va_.SetClipOffset(clip_offset_);
  Scene::Instance().RequestRedraw();
}

void GfxObj::ColorBy(const mol::EntityView& ev, 
                      const String& prop,
                      const Gradient& g, float minv, float maxv)
{
  LOG_VERBOSE("ColorBy not implemented for this gfx object");
}

void GfxObj::ColorBy(const img::MapHandle& mh,
                      const String& prop,
                      const Gradient& g, float minv, float maxv)
{
  LOG_VERBOSE("ColorBy not implemented for this gfx object");
}

//////////////////////////////////////////////////
// and now for the rest of the GfxObj interface

geom::AlignedCuboid GfxObj::GetBoundingBox(bool use_tf) const
{
  return geom::AlignedCuboid(geom::Vec3(),geom::Vec3());
}

void GfxObj::ProcessLimits(geom::Vec3& minc, geom::Vec3& maxc, 
                           const geom::Transform& tf) const
{
  try {
    geom::AlignedCuboid coord_limits=tf.Apply(this->GetBoundingBox(true));
    // update min/max by transforming all 8 corners of the bounding box and 
    // comparing it against the current min/max
    minc=geom::Min(minc,coord_limits.GetMin());
    maxc=geom::Min(maxc,coord_limits.GetMax());
  } catch(Error& e) {
    // in case the object is empty...
  }
}


void GfxObj::CustomRenderGL(RenderPass pass) {}

void GfxObj::CustomPreRenderGL(bool flag) {}

void GfxObj::CustomRenderPov(PovState& pov) {}

bool GfxObj::OnSelect(const geom::Line3& l, geom::Vec3& result, 
                      float zlim, bool pick_flag)
{
  return false;
}

void GfxObj::OnInput(const InputEvent& e) 
{
  geom::Mat3 rot=gfx::Scene::Instance().GetTransform().GetRot();
  geom::Vec3 cam_trans=gfx::Scene::Instance().GetTransform().GetTrans();  
  float x=e.GetDelta()*M_PI/180.0;
  bool transformed=false;
  if (e.GetCommand()==INPUT_COMMAND_ROTX) {
    transformed=true;
    transform_.SetRot(AxisRotation(geom::Vec3(1,0,0)*rot,x)*transform_.GetRot());
  } else if(e.GetCommand()==INPUT_COMMAND_ROTY) {
    transformed=true;
    transform_.SetRot(AxisRotation(geom::Vec3(0,1,0)*rot,x)*transform_.GetRot());
  } else if(e.GetCommand()==INPUT_COMMAND_ROTZ) {
    transformed=true;
    transform_.SetRot(AxisRotation(geom::Vec3(0,0,1)*rot,x)*transform_.GetRot());
  } else if(e.GetCommand()==INPUT_COMMAND_TRANSX ||
            e.GetCommand()==INPUT_COMMAND_TRANSY) {
 
    double mm[]={1,0,0,0,
                 0,1,0,0,
                 0,0,1,0,
                 transform_.GetTrans()[0]-cam_trans[0],
                 transform_.GetTrans()[1]-cam_trans[1],
                 transform_.GetTrans()[2]-cam_trans[2], 1};
    double pm[16];
    glGetDoublev(GL_PROJECTION_MATRIX,pm);
    GLint vp[4];
    glGetIntegerv(GL_VIEWPORT,vp);
    double wx,wy,wz;
    gluProject(0.0,0.0,0.0,mm,pm,vp,&wx,&wy,&wz);
    double ox,oy,oz;
    gluUnProject(wx+1.0,wy+1.0,wz,mm,pm,vp,&ox,&oy,&oz);

    float fx=ox;
    float fy=oy;   
    geom::Vec3 trans=transform_.GetTrans();
    transformed=true;    
    if (e.GetCommand()==INPUT_COMMAND_TRANSX) {
      trans+=geom::Vec3(-fx*e.GetDelta(), 0.0, 0.0)*rot;
    } else {
      trans+=geom::Vec3(0.0, -fy*e.GetDelta(), 0.0)*rot;
    }
    transform_.SetTrans(trans);
  } else if (e.GetCommand()==INPUT_COMMAND_TRANSZ) {
    float currz=transform_.GetTrans()[2];
    float delta=currz*pow(1.01f,-e.GetDelta())-currz;
    transform_.ApplyZAxisTranslation(delta);
    transformed=true;
  }
  if (transformed) {
    GfxObjP obj=dyn_cast<GfxObj>(this->shared_from_this());
    Scene::Instance().NodeTransformed(obj);
  }
}

const geom::Transform& GfxObj::GetTF() const
{
  return transform_;
}


void GfxObj::SetTF(const geom::Transform& tf)
{
  transform_=tf;
}

void GfxObj::FlagRebuild()
{
  rebuild_=true;
  Scene::Instance().RequestRedraw();
}

void GfxObj::FlagRefresh()
{
  refresh_=true;
  Scene::Instance().RequestRedraw();
}

void GfxObj::SetNormalSmoothFactor(float smoothf)
{
  smoothf_=smoothf;
  FlagRebuild();
}

float GfxObj::GetNormalSmoothFactor() const
{
  return smoothf_;
}

void GfxObj::OnRenderModeChange()
{
  Scene::Instance().ObjectChanged(GetName());
  Scene::Instance().RenderModeChanged(GetName());
}

void GfxObj::PreRenderGL(bool f)
{
  LOG_DEBUG("object " << GetName() << ": PreRenderGL()");
  CustomPreRenderGL(f);
}

void GfxObj::Clear()
{
  va_.Clear();
}

Color GfxObj::Ele2Color(const String& ele)
{
  // TODO 1: allow this table to be set by the user
  // TODO 2: move this table to entity gfx
  if(ele=="C") {return Color(0.83,0.97,0.97);}
  else if(ele=="N") {return Color(0.0,0.0,1.0);}
  else if(ele=="O") {return Color(1.0,0.0,0.0);}
  else if(ele=="S") {return Color(1.0,1.0,0.0);}
  else if(ele=="P") {return Color(1.0,0.0,1.0);}
  else if(ele=="H") {return Color(0.95,0.95,1.0);}
  
  return Color(0.8,0.8,0.8);
}

void GfxObj::AddLabel(const String& s, const geom::Vec3& pos, const Color& col, float psize)
{
  labels_.push_back(TextPrim(s,pos,col,psize));
}

void GfxObj::AddLabel(const String& s, const geom::Vec3& pos, const Color& col)
{
  AddLabel(s,pos,col,1.0);
}

void GfxObj::AddLabel(const String& s, const geom::Vec3& pos, float psize)
{
  AddLabel(s,pos,Color(1,1,1),psize);
}

void GfxObj::AddLabel(const String& s, const geom::Vec3& pos)
{
  AddLabel(s,pos,Color(1,1,1),1.0);
}

void GfxObj::ClearLabels()
{
  labels_.clear();
}


Material GfxObj::GetMaterial() const
{
  return mat_;
}

void GfxObj::SetMaterial(const Material& m)
{
  mat_=m;
}

void GfxObj::GLCleanup()
{
  if(mat_dlist_!=0) {
    glDeleteLists(mat_dlist_,1);
    mat_dlist_=0;
  }
  OnGLCleanup();
}

void GfxObj::OnGLCleanup()
{}

void GfxObj::SmoothVertices(float smoothf)
{
  va_.SmoothVertices(smoothf);
  FlagRefresh();
}

namespace {

/*
commented out to silence compiler warning, as it is not needed anywhere
float normalize(float v, float min_v, float max_v)
{
  return (v-min_v)/(max_v-min_v);
}
*/

}

void GfxObj::render_labels() const
{
  if(labels_.empty()) return;
  for(TextPrimList::const_iterator it=labels_.begin();it!=labels_.end();++it) {
    Scene::Instance().RenderText(*it);
  }
}

void GfxObj::ReapplyColorOps(){
  if(! c_ops_.empty()){
    GfxObjP o=dyn_cast<GfxObj>(shared_from_this());
    for(boost::ptr_vector<gfx::ColorOp>::iterator it=c_ops_.begin();
      it!=c_ops_.end();++it) {
      it->ApplyTo(o);
    }
  }
}

void GfxObj::AppendColorOp(gfx::ColorOp* op){
  c_ops_.push_back(op);
}

void GfxObj::CleanColorOps(){
  c_ops_.release();
}

void GfxObj::Debug(unsigned int flags)
{
  debug_flags_=flags;
  va_.DrawNormals(flags & 0x1);
}

void GfxObj::render_depth_only()
{
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glPushClientAttrib(GL_ALL_ATTRIB_BITS);
  glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
  glDisable(GL_LIGHTING);
  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_BLEND);
  glDisable(GL_FOG);
  glShadeModel(GL_FLAT);
  glDisable(GL_LINE_SMOOTH);
  glDisable(GL_POINT_SMOOTH);
  glDisable(GL_POLYGON_SMOOTH);
  glDisable(GL_DITHER);
  glDisable(GL_NORMALIZE);
  glEnable(GL_DEPTH_TEST);
  CustomRenderGL(STANDARD_RENDER_PASS);
  CustomRenderGL(TRANSPARENT_RENDER_PASS);
  glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);      
  glPopClientAttrib();
  glPopAttrib();
}


}} // ns
