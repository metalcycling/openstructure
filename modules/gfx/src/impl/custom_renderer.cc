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
#include "custom_renderer.hh"

#include <ost/gfx/entity.hh>
#include <ost/gfx/povray.hh>

namespace ost {

namespace gfx {

using namespace impl;

CustomRenderer::CustomRenderer(): options_(new CustomRenderOptions()) 
{
  this->SetName("Ball & Stick");
  this->SetFixedPickRadius(0.2);
}

void CustomRenderer::PrepareRendering()
{
  ConnectRendererBase::PrepareRendering();
  va_.Clear();
  this->PrepareRendering(view_, va_, false);
  sel_va_.Clear();
  if (this->HasSelection()) {
    this->PrepareRendering(sel_view_, sel_va_, true);
    sel_va_.SetLighting(false);
    sel_va_.SetColorMaterial(true);    
  }
}

void CustomRenderer::PrepareRendering(GfxView& view, 
                                      IndexedVertexArray& va, bool is_sel)
{
  float plus=is_sel ? 0.05: 0.0;
  if(options_!=NULL){
    va.SetLighting(true);
    va.SetCullFace(true);
    va.SetColorMaterial(true);
    va.SetMode(0x4);
    // draw all spheres
    for (AtomEntryMap::const_iterator it=view.atom_map.begin();it!=view.atom_map.end();++it) {
      va.AddSphere(SpherePrim(it->second.atom.GetPos(),
                              options_->GetSphereRad()+plus,
                              it->second.color),
                              options_->GetSphereDetail());
    }

    for(BondEntryList::const_iterator it=view.bond_list.begin();it!=view.bond_list.end();++it) {
      const geom::Vec3& p0=it->bond.GetFirst().GetPos();
      const geom::Vec3& p2=it->bond.GetSecond().GetPos();
      geom::Vec3 p1=(p0+p2)*0.5;
      va.AddCylinder(CylinderPrim(p0,p1,options_->GetBondRad()+plus,it->atom1->color),
                     options_->GetArcDetail());
      va.AddCylinder(CylinderPrim(p1,p2,options_->GetBondRad()+plus,it->atom2->color),
                     options_->GetArcDetail());
    }
  }
  sel_state_=0;
  state_=0;
}

void CustomRenderer::RenderPov(PovState& pov, const std::string& name)
{
  if(view_.atom_map.empty() && view_.bond_list.empty()) return;
  pov.write_merge_or_union(name);
  
  for (AtomEntryMap::const_iterator it=view_.atom_map.begin();it!=view_.atom_map.end();++it) {
    pov.write_sphere(it->second.atom.GetPos(),
                     options_->GetSphereRad(),
                     it->second.color,
                     name);
  }

  for(BondEntryList::const_iterator it=view_.bond_list.begin();it!=view_.bond_list.end();++it) {
    const geom::Vec3& p0=it->bond.GetFirst().GetPos();
    const geom::Vec3& p2=it->bond.GetSecond().GetPos();
    geom::Vec3 p1=(p0+p2)*0.5;
    pov.write_cyl(p0,p1,options_->GetBondRad(),it->atom1->color,name,true);
    pov.write_cyl(p1,p2,options_->GetBondRad(),it->atom2->color,name,true);
  }
  pov.inc() << " }\n";
}

bool CustomRenderer::CanSetOptions(RenderOptionsPtr& render_options)
{
 return render_options->GetRenderMode()==RenderMode::CUSTOM;
}

void CustomRenderer::SetOptions(RenderOptionsPtr& render_options)
{
  options_ = boost::static_pointer_cast<CustomRenderOptions>(render_options);
}

RenderOptionsPtr CustomRenderer::GetOptions()
{
 return options_;
}

int CustomRenderer::GetMode()
{
  return 0;
}

CustomRenderer::~CustomRenderer() 
{
}

}} // ns
