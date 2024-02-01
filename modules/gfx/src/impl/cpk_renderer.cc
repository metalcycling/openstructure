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
#include "cpk_renderer.hh"

#include <ost/gfx/gl_helper.hh>

#include <ost/gfx/entity.hh>

#include <ost/gfx/scene.hh>
#if OST_SHADER_SUPPORT_ENABLED
#include <ost/gfx/shader.hh>
#endif //OST_SHADER_SUPPORT_ENABLED
#include <ost/gfx/vertex_array_helper.hh>
#include <ost/gfx/render_mode.hh>
#include <ost/gfx/povray.hh>

namespace ost { namespace gfx { namespace impl {

CPKRenderer::CPKRenderer(): 
  options_(new CPKRenderOptions())
#if OST_SHADER_SUPPORT_ENABLED
  ,fsr_(),sel_fsr_()
#endif  
{
  this->SetName("Spheres");
}

void CPKRenderer::PrepareRendering()
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

void CPKRenderer::PrepareRendering(GfxView& view, IndexedVertexArray& va, bool is_sel)
{
  RGBAColor sel_clr=this->GetSelectionColor();
  float factor=is_sel ? 1.2 : 1.0;
  if(options_!=NULL){
    factor *= options_->GetRadiusMult();
    if(factor>0.0) {
#if OST_SHADER_SUPPORT_ENABLED
      if(options_->GetSphereMode()==1 || options_->GetSphereMode()==2) {
        FastSphereRenderer& fsr = is_sel ? sel_fsr_ : fsr_;
        fsr.Clear();
        for(AtomEntryMap::const_iterator it=view.atom_map.begin();it!=view.atom_map.end();++it) {
          fsr.Add(it->second.atom.GetPos(),
                  is_sel? sel_clr : it->second.color,
                  it->second.vdwr*factor);
        }
      } else {
#endif
        va.SetLighting(true);
        va.SetCullFace(true);
        va.SetColorMaterial(true);
        va.SetMode(0x4);
        
        // draw all spheres
        uint det=options_->GetSphereDetail();
        for(AtomEntryMap::const_iterator it=view.atom_map.begin();it!=view.atom_map.end();++it) {
          va.AddSphere(SpherePrim(it->second.atom.GetPos(),
                                  it->second.vdwr*factor,
                                  is_sel? sel_clr : it->second.color),
                       det);
        }
      }
#if OST_SHADER_SUPPORT_ENABLED
    }
#endif
  }
  sel_state_=0;
  state_=0;
}

void CPKRenderer::Render(RenderPass pass)
{
  if(options_!=NULL){
#if OST_SHADER_SUPPORT_ENABLED
    if(options_->GetSphereMode()==1 || options_->GetSphereMode()==2) {
      if(pass==STANDARD_RENDER_PASS) {
        fsr_.RenderGL();
      } else if(pass==GLOW_RENDER_PASS && this->HasSelection()) {
        sel_fsr_.RenderGL();
      }
      return;
    }
#endif
  }
  // fall back to default if above did not fire
  EntityRenderer::Render(pass);
}

void CPKRenderer::RenderPov(PovState& pov, const std::string& name)
{
  if(view_.atom_map.empty()) return;
  pov.write_merge_or_union(name);
  
  for (AtomEntryMap::const_iterator it=view_.atom_map.begin();it!=view_.atom_map.end();++it) {
    pov.write_sphere(it->second.atom.GetPos(),
                     it->second.vdwr,
                     it->second.color,
                     name);
  }

  pov.inc() << " }\n";
}

bool CPKRenderer::CanSetOptions(RenderOptionsPtr& render_options){
 return render_options.get()->GetRenderMode()==RenderMode::CPK;
}

void CPKRenderer::SetOptions(RenderOptionsPtr& render_options){
  options_ = boost::static_pointer_cast<CPKRenderOptions>(render_options);
}

RenderOptionsPtr CPKRenderer::GetOptions(){
 return options_;
}

}}}
