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

#include <ost/mol/atom_handle.hh>

#include "scene.hh"

#include "prim_list.hh"
#include "povray.hh"

namespace ost { namespace gfx {

PrimList::PrimList(const String& name):
  GfxObj(name),
  points_(),
  lines_(),
  spheres_(),
  cyls_(),
  texts_(),
  sphere_detail_(4),
  arc_detail_(4),
  simple_va_(),
  vas_()
{}

void PrimList::Clear()
{
  points_.clear();
  lines_.clear();
  spheres_.clear();
  cyls_.clear();
  texts_.clear();
  vas_.clear();
  Scene::Instance().RequestRedraw();
  this->FlagRebuild();
}

geom::AlignedCuboid PrimList::GetBoundingBox(bool use_tf) const
{
  if(points_.empty() && lines_.empty() && spheres_.empty() && cyls_.empty() && texts_.empty() && vas_.empty()) {
    return geom::AlignedCuboid(geom::Vec3(-1,-1,-1),geom::Vec3(1,1,1));
  }
  geom::Vec3 minc(std::numeric_limits<float>::max(),
                  std::numeric_limits<float>::max(),
                  std::numeric_limits<float>::max());
  geom::Vec3 maxc(-std::numeric_limits<float>::max(),
                  -std::numeric_limits<float>::max(),
                  -std::numeric_limits<float>::max());
  ProcessLimits(minc,maxc,geom::Transform());
  geom::AlignedCuboid bb(minc,maxc);
  return use_tf ? transform_.Apply(bb) : bb;
}

void PrimList::ProcessLimits(geom::Vec3& minc, geom::Vec3& maxc, 
                             const geom::Transform& tf) const
{
  for(SpherePrimList::const_iterator it=points_.begin();it!=points_.end();++it) {
    geom::Vec3 tpos = tf.Apply(it->position);
    minc=geom::Min(minc,tpos);
    maxc=geom::Max(maxc,tpos);
  }
  for(CylinderPrimList::const_iterator it=lines_.begin();it!=lines_.end();++it) {
    geom::Vec3 tpos = tf.Apply(it->start);
    minc=geom::Min(minc,tpos);
    maxc=geom::Max(maxc,tpos);
    tpos = tf.Apply(it->end);
    minc=geom::Min(minc,tpos);
    maxc=geom::Max(maxc,tpos);
  }
  for(SpherePrimList::const_iterator it=spheres_.begin();it!=spheres_.end();++it) {
    geom::Vec3 tpos = tf.Apply(it->position);
    minc=geom::Min(minc,tpos-it->radius);
    maxc=geom::Max(maxc,tpos+it->radius);
  }
  for(CylinderPrimList::const_iterator it=cyls_.begin();it!=cyls_.end();++it) {
    geom::Vec3 tpos = tf.Apply(it->start);
    minc=geom::Min(minc,tpos-it->radius1);
    maxc=geom::Max(maxc,tpos+it->radius1);
    tpos = tf.Apply(it->end);
    minc=geom::Min(minc,tpos-it->radius2);
    maxc=geom::Max(maxc,tpos+it->radius2);
  }
  for(TextPrimList::const_iterator it=texts_.begin();it!=texts_.end();++it) {
    geom::Vec3 tpos = tf.Apply(it->position);
    minc=geom::Min(minc,tpos);
    maxc=geom::Max(maxc,tpos);
  }
  for(std::vector<IndexedVertexArray>::const_iterator it=vas_.begin();it!=vas_.end();++it) {
    geom::AlignedCuboid bb=it->GetBoundingBox();
    minc=geom::Min(minc,bb.GetMin());
    maxc=geom::Max(maxc,bb.GetMax());
  }
  minc-=1.0;
  maxc+=1.0;
}

geom::Vec3 PrimList::GetCenter() const
{
  geom::Vec3 cen;
  size_t sum=0;
  for(SpherePrimList::const_iterator it=points_.begin();it!=points_.end();++it) {
    cen+=it->position;
  }
  sum+=points_.size();
  for(CylinderPrimList::const_iterator it=lines_.begin();it!=lines_.end();++it) {
    cen+=0.5*(it->start+it->end);
  }
  sum+=lines_.size();
  for(SpherePrimList::const_iterator it=spheres_.begin();it!=spheres_.end();++it) {
    cen+=it->position;
  }
  sum+=spheres_.size();
  for(CylinderPrimList::const_iterator it=cyls_.begin();it!=cyls_.end();++it) {
    cen+=0.5*(it->start+it->end);
  }
  sum+=cyls_.size();
  for(TextPrimList::const_iterator it=texts_.begin();it!=texts_.end();++it) {
    cen+=it->position;
  }
  sum+=texts_.size();
  if(sum>0) {
    cen/=static_cast<float>(sum);
  }
  return cen;
}

void PrimList::SetLineWidth(float w)
{
  GfxObj::SetLineWidth(w);
  va_.SetLineWidth(line_width_);
  simple_va_.SetLineWidth(line_width_);
  for(std::vector<IndexedVertexArray>::iterator it=vas_.begin();it!=vas_.end();++it) {
    it->SetLineWidth(line_width_);
  }
}

void PrimList::OnRenderModeChange()
{
  // noop
}

void PrimList::CustomPreRenderGL(bool flag)
{
  prep_va();
  prep_simple_va();
}

void PrimList::CustomRenderGL(RenderPass pass)
{
  if(pass==STANDARD_RENDER_PASS) {
    va_.RenderGL();
    simple_va_.RenderGL();
    render_text();
    for(std::vector<IndexedVertexArray>::iterator it=vas_.begin();it!=vas_.end();++it) {
      it->RenderGL();
    }
  }
}

void PrimList::CustomRenderPov(PovState& pov)
{
  if(points_.empty() && lines_.empty() && spheres_.empty() && cyls_.empty()) return;
  pov.write_merge_or_union(GetName());

  for(SpherePrimList::const_iterator it=points_.begin();it!=points_.end();++it) {
    pov.write_sphere(it->position,0.1,it->color,GetName());
  }
  for(CylinderPrimList::const_iterator it=lines_.begin();it!=lines_.end();++it) {
    pov.write_sphere(it->start,0.1,it->color1,GetName());
    pov.write_sphere(it->end,0.1,it->color2,GetName());
    pov.write_cyl(it->start,it->end,0.1,it->color1,GetName(),true);
  }
  for(SpherePrimList::const_iterator it=spheres_.begin();it!=spheres_.end();++it) {
    pov.write_sphere(it->position,it->radius,it->color,GetName());
  }
  for(CylinderPrimList::const_iterator it=cyls_.begin();it!=cyls_.end();++it) {
    pov.write_sphere(it->start,it->radius1,it->color1,GetName());
    pov.write_sphere(it->end,it->radius2,it->color2,GetName());
    pov.write_cyl(it->start,it->end,it->radius1,it->color1,GetName(),true);
  }
  // TODO text
  pov.inc() << " }\n";
}

void PrimList::AddPoint(const geom::Vec3& p, const Color& col)
{
  points_.push_back(SpherePrim(p, 0.0, col));
  Scene::Instance().RequestRedraw();
  FlagRebuild();
}

void PrimList::AddLine(const geom::Vec3& p1, const geom::Vec3& p2, const Color& col1, const Color& col2)
{
  lines_.push_back(CylinderPrim(p1,p2,0.0,0.0,col1,col2));
  Scene::Instance().RequestRedraw();
  FlagRebuild();
}

void PrimList::AddSphere(const geom::Vec3& c, float r, const Color& col)
{
  spheres_.push_back(SpherePrim(c, r, col));
  Scene::Instance().RequestRedraw();
  FlagRebuild();
}

void PrimList::AddCyl(const geom::Vec3& p1, const geom::Vec3& p2, float r1, float r2, const Color& col1, const Color& col2)
{
  cyls_.push_back(CylinderPrim(p1, p2, r1, r2, col1, col2));
  Scene::Instance().RequestRedraw();
  FlagRebuild();
}

void PrimList::AddText(const std::string& text, const geom::Vec3& pos, const Color& col, float point_size)
{
  texts_.push_back(TextPrim(text,pos,col,point_size));
  Scene::Instance().RequestRedraw();
  FlagRebuild();
}

void PrimList::SetDiameter(float d)
{
  LOG_WARNING("PrimList::SetDiameter is defunct");
}

void PrimList::SetRadius(float r)
{
  LOG_WARNING("PrimList::SetDiameter is defunct");
}

void PrimList::SetSphereDetail(unsigned int d)
{
  sphere_detail_=d;
  Scene::Instance().RequestRedraw();
  FlagRebuild();
}

void PrimList::SetArcDetail(unsigned int d)
{
  arc_detail_=d;
  Scene::Instance().RequestRedraw();
  FlagRebuild();
}

void PrimList::SetColor(const Color& c)
{
  for(SpherePrimList::iterator it=points_.begin();it!=points_.end();++it) {
    it->color=c;
  }
  for(CylinderPrimList::iterator it=lines_.begin();it!=lines_.end();++it) {
    it->color1=c;
    it->color2=c;
  }
  for(SpherePrimList::iterator it=spheres_.begin();it!=spheres_.end();++it) {
    it->color=c;
  }
  for(CylinderPrimList::iterator it=cyls_.begin();it!=cyls_.end();++it) {
    it->color1=c;
    it->color2=c;
  }
  for(TextPrimList::iterator it=texts_.begin();it!=texts_.end();++it) {
    it->color=c;
  }
  Scene::Instance().RequestRedraw();
  FlagRebuild();
}

void PrimList::AddMesh(float* v, float* n, float* c, size_t nv, unsigned int* i, size_t ni)
{
  static float dummy_normal[]={0.0,0.0,1.0};
  static float dummy_color[]={1.0,1.0,1.0,1.0};
  vas_.push_back(IndexedVertexArray());
  IndexedVertexArray& va=vas_.back();
  va.SetLighting(true);
  va.SetTwoSided(true);
  va.SetColorMaterial(true);
  va.SetCullFace(false);
  float* vv=v;
  float* nn=n;
  if(!n) {
    nn=dummy_normal;
    va.SetLighting(false);
  }
  float* cc=c;
  if(!c) {
    cc=dummy_color;
  }
  for(size_t k=0;k<nv;++k) {
    va.Add(geom::Vec3(vv[0],vv[1],vv[2]),
           geom::Vec3(nn[0],nn[1],nn[2]),
           Color(cc[0],cc[1],cc[2],cc[3]));
    vv+=3;
    if(n) nn+=3;
    if(c) cc+=4;
  }
  unsigned int* ii=i;
  for(size_t k=0;k<ni;++k) {
    va.AddTri(ii[0],ii[1],ii[2]);
    ii+=3;
  }
  Scene::Instance().RequestRedraw();
  FlagRebuild();
}

////////////////////////////////
// private methods

void PrimList::prep_simple_va()
{
  simple_va_.Clear();
  simple_va_.SetLighting(false);
  simple_va_.SetCullFace(false);
  simple_va_.SetColorMaterial(false);
  simple_va_.SetMode(0x3);
  simple_va_.SetTwoSided(true);
  simple_va_.SetAALines(GetAALines());
  simple_va_.SetOpacity(GetOpacity());
  simple_va_.SetLineWidth(line_width_);

  for(SpherePrimList::const_iterator it=points_.begin();it!=points_.end();++it) {
    simple_va_.Add(it->position,geom::Vec3(),it->color);
  }

  for(CylinderPrimList::const_iterator it=lines_.begin();it!=lines_.end();++it) {
    VertexID id0 = simple_va_.Add(it->start,geom::Vec3(),it->color1);
    VertexID id1 = simple_va_.Add(it->end,geom::Vec3(),it->color2);
    simple_va_.AddLine(id0,id1);
  }
}

void PrimList::prep_va()
{
  va_.Clear();
  va_.SetLighting(true);
  va_.SetCullFace(true);
  va_.SetColorMaterial(true);
  va_.SetMode(0x4);
  va_.SetLineWidth(line_width_);

  for(SpherePrimList::const_iterator it=spheres_.begin();it!=spheres_.end();++it) {
    va_.AddSphere(SpherePrim(it->position, it->radius, it->color), 
                  GetSphereDetail());
  }

  for(CylinderPrimList::const_iterator it=cyls_.begin();it!=cyls_.end();++it) {
    va_.AddCylinder(CylinderPrim(it->start, it->end, it->radius1, it->radius2, it->color1, it->color2),
                    GetArcDetail(),
                    true);
  }
}

void PrimList::render_text()
{
  for(TextPrimList::const_iterator it=texts_.begin();it!=texts_.end();++it) {
    Scene::Instance().RenderText(*it);
  }
}

}} // ns
