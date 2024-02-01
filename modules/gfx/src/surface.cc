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

#include <map>
#include <limits>

#include <ost/mol/mol.hh>
#include <ost/mol/entity_property_mapper.hh>
#include "surface.hh"
#include "glext_include.hh"
#include "impl/mapped_property.hh"

namespace ost { 

using namespace geom;

namespace gfx {

Surface::Surface(const String& name, const mol::SurfaceHandle& sh):
  GfxObj(name),
  sh_(sh),
  vmap_(),
  recalc_bb_(true)
{
  // TODO replace with def mat for this gfx obj type
  SetMatAmb(Color(0,0,0));
  SetMatDiff(Color(1,1,1));
  SetMatSpec(Color(0.2,0.2,0.2));
  SetMatShin(48);
  Rebuild();
}

void Surface::CustomRenderGL(RenderPass pass)
{
  if(pass==STANDARD_RENDER_PASS) {
    va_.RenderGL();
  }
}

void Surface::CustomRenderPov(PovState& pov)
{
  va_.RenderPov(pov, GetName());
}

void Surface::OnRenderModeChange()
{
  Rebuild();
  GfxObj::OnRenderModeChange();
}

void Surface::Replace(const mol::SurfaceHandle& sh)
{
  sh_=sh;
  recalc_bb_=true;
  Rebuild();
}

void Surface::Rebuild()
{
  if(GetRenderMode()==RenderMode::DEBUG2) {
    vmap_.clear();
    va_.Clear();
    va_.SetMode(0x7);
    va_.SetPolyMode(1);
    va_.SetLighting(false);
    va_.SetCullFace(false);
    va_.SetPointSize(3.0);
    
    std::vector<mol::SurfaceVertexID> svid_list = sh_.GetVertexIDList();
    std::vector<mol::SurfaceTriID> stid_list = sh_.GetTriIDList();
    
    for(std::vector<mol::SurfaceVertexID>::const_iterator it=svid_list.begin();
        it!=svid_list.end();++it) {
      mol::SurfaceVertex vert=sh_.GetVertex(*it);
      Color col(1.0,1.0,1.0,1.0);
      if(vert.type==1) {
        col=Color(1.0,0.0,0.0);
      } else if(vert.type==2) {
        col=Color(0.0,1.0,0.0);
      } else if(vert.type==3) {
        col=Color(0.0,0.0,1.0);
      }
      vmap_[*it] = va_.Add(vert.position,vert.normal,col);
    }
    
    for(std::vector<mol::SurfaceTriID>::const_iterator it=stid_list.begin();
        it!=stid_list.end();++it) {
      mol::SurfaceTri tri=sh_.GetTri(*it);
      va_.AddTri(vmap_[tri.v0],vmap_[tri.v1],vmap_[tri.v2]);
    }
  } else {
    va_.Clear();
    va_.SetMode(0x4);
    va_.SetLighting(true);
    va_.SetTwoSided(true);
    va_.SetCullFace(false);

    if(GetRenderMode()==RenderMode::DEBUG) {
      va_.DrawNormals(true);
    } else {
      va_.DrawNormals(false);
    }
    
    std::vector<mol::SurfaceVertexID> svid_list = sh_.GetVertexIDList();
    std::vector<mol::SurfaceTriID> stid_list = sh_.GetTriIDList();
    
    for(std::vector<mol::SurfaceVertexID>::const_iterator it=svid_list.begin();
        it!=svid_list.end();++it) {
      mol::SurfaceVertex vert=sh_.GetVertex(*it);
      Color col(1.0,1.0,1.0,1.0);
      vmap_[*it] = va_.Add(vert.position,vert.normal,col);
    }
    
    for(std::vector<mol::SurfaceTriID>::const_iterator it=stid_list.begin();
        it!=stid_list.end();++it) {
      mol::SurfaceTri tri=sh_.GetTri(*it);
      va_.AddTri(vmap_[tri.v0],vmap_[tri.v1],vmap_[tri.v2]);
    }
  }
  this->ReapplyColorOps();
  recalc_bb_=true;
}

geom::AlignedCuboid Surface::GetBoundingBox(bool use_tf) const
{
  static geom::AlignedCuboid bb;
  if(recalc_bb_) {
    geom::Vec3 minc(std::numeric_limits<float>::max(),
                    std::numeric_limits<float>::max(),
                    std::numeric_limits<float>::max());
    geom::Vec3 maxc(-std::numeric_limits<float>::max(),
                    -std::numeric_limits<float>::max(),
                    -std::numeric_limits<float>::max());
    
    std::vector<mol::SurfaceVertexID> svid_list = sh_.GetVertexIDList();
    if(! svid_list.empty()) {
      for(std::vector<mol::SurfaceVertexID>::const_iterator it=svid_list.begin();
          it!=svid_list.end();++it) {
        geom::Vec3 pos = sh_.GetVertex(*it).position;
        minc=geom::Min(minc,pos);
        maxc=geom::Max(maxc,pos);
      }
      minc-=geom::Vec3(1.0,1.0,1.0);
      maxc+=geom::Vec3(1.0,1.0,1.0);
      recalc_bb_=false;
      bb=geom::AlignedCuboid(minc,maxc);
    }
  }
  return use_tf ? transform_.Apply(bb) : bb;
}

void Surface::CustomPreRenderGL(bool flag)
{
  va_.FlagRefresh();
  if(flag) {
    //Rebuild();
  } else {
    //RefreshVA(va_);
  }
}

void Surface::SetColor(const Color& col, const String& selection)
{
  UniformColorOp cop = UniformColorOp(selection, col);
  this->Apply(cop);
}

void Surface::ColorBy(const mol::EntityView& ev, 
                      const String& prop,
                      const Gradient& g, float minv, float maxv)
{
  EntityViewColorOp evop = EntityViewColorOp(prop, g,minv,maxv,ev);
  this->Apply(evop);
}

void Surface::ColorBy(const img::MapHandle& mh,
                      const String& prop,
                      const Gradient& g,float minv, float maxv)
{
  MapHandleColorOp mhop = MapHandleColorOp("",prop,g,minv,maxv,mh);
  this->Apply(mhop);
}

namespace {

float normalize(float v, float min_v, float max_v)
{
  return (v-min_v)/(max_v-min_v);
}
}

void Surface::ColorBy(const String& prop, 
                     const Gradient& gradient,
                     float minv,float maxv,
                     mol::Prop::Level level)
{
  GradientLevelColorOp glop = GradientLevelColorOp("",prop, gradient,minv,maxv,level);
  this->Apply(glop);
}

void Surface::ColorBy(const String& prop, 
                     const Color& c1, const Color& c2, 
                     float minv, float maxv,
                     mol::Prop::Level level)
{
  Gradient g;
  g.SetColorAt(0.0,c1);
  g.SetColorAt(1.0,c2);
  this->ColorBy(prop,g,minv,maxv,level);
}

void Surface::ColorBy(const String& prop, 
                     const Gradient& gradient,
                     mol::Prop::Level level)
{
  BasicGradientColorOp bgop = BasicGradientColorOp(prop, gradient,level);
  this->Apply(bgop);
}

void Surface::ColorBy(const String& prop, 
                     const Color& c1, const Color& c2,
                     mol::Prop::Level level)
{
  Gradient g;
  g.SetColorAt(0.0,c1);
  g.SetColorAt(1.0,c2);
  this->ColorBy(prop,g,level);
}

void Surface::Apply(const gfx::UniformColorOp& op, bool store){
  if(store){
    UniformColorOp* op_ptr = new UniformColorOp(op);
    this->AppendColorOp(op_ptr);
  }
  mol::Query q(op.GetSelection());
  gfx::Color col = op.GetColor();
  if(op.GetSelection()=="") {
    for(VMap::const_iterator it=vmap_.begin();it!=vmap_.end();++it) {
      va_.SetColor(it->second,col);
    }
  } else {
    for(VMap::const_iterator it=vmap_.begin();it!=vmap_.end();++it) {
      mol::AtomHandle ah = sh_.GetVertex(it->first).atom;
      if(ah.IsValid()) {
        if(q.IsAtomSelected(ah)) {
          va_.SetColor(it->second,col);
        }
      }
    }
  }
  FlagRefresh();
}

void Surface::Apply(const gfx::BasicGradientColorOp& op, bool store)
{
  if(store){
    BasicGradientColorOp* op_ptr = new BasicGradientColorOp(op);
    this->AppendColorOp(op_ptr);
  }
  const String& prop = op.GetSelection();
  const Gradient& gradient = op.GetGradient();
  mol::Prop::Level level = op.GetLevel();
  // get the min/max values first, requires one lookup round
  float minv=std::numeric_limits<float>::max();
  float maxv=-std::numeric_limits<float>::max();
  mol::EntityPropertyMapper epm(prop, level);
  std::vector<std::pair<VertexID,float> > v2v;
  for(VMap::const_iterator it=vmap_.begin();it!=vmap_.end();++it) {
    mol::AtomHandle ah = sh_.GetVertex(it->first).atom;
    if(ah.IsValid()) {
      float v = epm.Get(ah);
      v2v.push_back(std::make_pair(it->second,v));
      minv=std::min(minv,v);
      maxv=std::max(maxv,v);
    }
  }
  
  // reuse values for speed optimization
  for(std::vector<std::pair<VertexID,float> >::const_iterator it=v2v.begin();it!=v2v.end();++it) {
    va_.SetColor(it->first,gradient.GetColorAt(normalize(it->second,minv,maxv)));
  }
  FlagRefresh();
}

void Surface::Apply(const gfx::GradientLevelColorOp& op, bool store)
{
  if(store){
    GradientLevelColorOp* op_ptr = new GradientLevelColorOp(op);
    this->AppendColorOp(op_ptr);
  }
  mol::Prop::Level level = op.GetLevel();
  const String& prop = op.GetProperty();
  const Gradient& gradient = op.GetGradient();
  float minv = op.GetMinV();
  float maxv = op.GetMaxV();
  // don't bother mapping any surface specific properties, just go
  // for the attached atoms
  mol::EntityPropertyMapper epm(prop, level);
  for(VMap::const_iterator it=vmap_.begin();it!=vmap_.end();++it) {
    mol::AtomHandle ah = sh_.GetVertex(it->first).atom;
    if(ah.IsValid()) {
      va_.SetColor(it->second,gradient.GetColorAt(normalize(epm.Get(ah,minv),minv,maxv)));
    }
  }
  FlagRefresh();
}

void Surface::Apply(const gfx::EntityViewColorOp& op, bool store)
{
  if(store){
    EntityViewColorOp* op_ptr = new EntityViewColorOp(op);
    this->AppendColorOp(op_ptr);
  }
  const mol::EntityView& ev = op.GetEntityView();
  const String& prop = op.GetProperty();
  const Gradient& g = op.GetGradient();
  float minv = op.GetMinV();
  float maxv = op.GetMaxV();
  for(VMap::const_iterator it=vmap_.begin();it!=vmap_.end();++it) {
    va_.SetColor(it->second,impl::MappedProperty(ev,prop,g,minv,maxv,
                                                 sh_.GetVertex(it->first).position));
  }
  FlagRefresh();
}

void Surface::Apply(const gfx::MapHandleColorOp& op, bool store)
{
  if(store){
    MapHandleColorOp* op_ptr = new MapHandleColorOp(op);
    this->AppendColorOp(op_ptr);
  }
  const img::MapHandle& mh = op.GetMapHandle();
  const String& prop = op.GetProperty();
  const Gradient& g = op.GetGradient();
  float minv = op.GetMinV();
  float maxv = op.GetMaxV();
  for(VMap::const_iterator it=vmap_.begin();it!=vmap_.end();++it) {
    va_.SetColor(it->second,impl::MappedProperty(mh,prop,g,minv,maxv,sh_.GetVertex(it->first).position));
  }
  FlagRefresh();
}

void Surface::CleanColorOps()
{
  GfxObj::CleanColorOps();
}

void Surface::ReapplyColorOps()
{
  GfxObj::ReapplyColorOps();
}

ost::mol::SurfaceHandle Surface::GetHandle() const
{
  return this->sh_;
}

}} // ns
