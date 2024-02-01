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
   Authors: Marco Biasini, Ansgar Philippsen
*/

#include <ost/gfx/entity.hh>
#include <ost/gfx/povray.hh>

#include "trace_renderer.hh"


namespace ost { namespace gfx { namespace impl {

TraceRenderer::TraceRenderer(BackboneTrace* trace): 
  TraceRendererBase(trace, 1), options_(new TraceRenderOptions())
{
  this->SetName("Trace");
}

void TraceRenderer::PrepareRendering() 
{
  TraceRendererBase::PrepareRendering(false);
  va_.Clear();
  this->PrepareRendering(trace_subset_, va_, false);
  sel_va_.Clear();
  if (this->HasSelection()) {
    //this->PrepareRendering(sel_subset_, sel_va_, true);
    this->PrepareRendering(trace_subset_, sel_va_, true);
    sel_va_.SetLighting(false);
  }
}

void TraceRenderer::PrepareRendering(BackboneTrace& trace_subset,
                                     IndexedVertexArray& va, bool is_sel)
{
  const Color& sel_clr=this->GetSelectionColor();
  if(options_!=NULL){
    va.Clear();
    va.SetLighting(true);
    va.SetCullFace(true);
    va.SetColorMaterial(true);
    va.SetTwoSided(false);
    if(is_sel) {
      va.SetOpacity(GetSelectionColor().Alpha());
      for (int node_list=0; node_list<trace_subset.GetListCount(); ++node_list) {
        const NodeEntryList& nl=trace_subset.GetList(node_list);
        for(unsigned int i=0;i<nl.size();++i) {
          mol::AtomHandle ah=nl[i].atom;
          if(sel_.ViewForHandle(ah).IsValid()) {
            geom::Vec3 apos = ah.GetPos();
            va.AddSphere(SpherePrim(apos,
                                    options_->GetTubeRadius()+0.05,
                                    sel_clr),
                         options_->GetArcDetail());
            if(i>0) {
              va.AddCylinder(CylinderPrim(apos+0.5*(nl[i-1].atom.GetPos()-apos),
                                          apos,
                                          options_->GetTubeRadius()+0.05,
                                          sel_clr),
                             options_->GetArcDetail());
            }
            if(i<nl.size()-1) {
              va.AddCylinder(CylinderPrim(apos,
                                          apos+0.5*(nl[i+1].atom.GetPos()-apos),
                                          options_->GetTubeRadius()+0.05,
                                          sel_clr),
                             options_->GetArcDetail());
            }
          }
        }
      }
    } else {
      for (int node_list=0; node_list<trace_subset.GetListCount(); ++node_list) {
        const NodeEntryList& nl=trace_subset.GetList(node_list);
        mol::AtomHandle a1=nl[0].atom;
        va.AddSphere(SpherePrim(a1.GetPos(),
                                options_->GetTubeRadius(),
                                nl[0].color1),
                     options_->GetArcDetail());
        for(unsigned int i=1;i<nl.size();++i) {
          mol::AtomHandle a2=nl[i].atom;
          va.AddSphere(SpherePrim(a2.GetPos(),
                                  options_->GetTubeRadius(),
                                  nl[i].color1),
                       options_->GetArcDetail());
          const geom::Vec3& p0=a1.GetPos();
          const geom::Vec3& p2=a2.GetPos();
          geom::Vec3 p1=(p0+p2)*0.5;
          va.AddCylinder(CylinderPrim(p0,p1,options_->GetTubeRadius(),nl[i-1].color1),
                         options_->GetArcDetail());
          va.AddCylinder(CylinderPrim(p1,p2,options_->GetTubeRadius(),nl[i].color1),
                         options_->GetArcDetail());
          a1=a2;
        }
      }
    }
  }
  sel_state_=0;
  state_=0;
}

void TraceRenderer::RenderPov(PovState& pov, const std::string& name)
{
  if(trace_subset_.GetListCount()==0) return;
  pov.write_merge_or_union(name);

  for (int node_list=0; node_list<trace_subset_.GetListCount(); ++node_list) {
    const NodeEntryList& nl=trace_subset_.GetList(node_list);

    geom::Vec3 p0=nl[0].atom.GetPos();
    float rad0=0.2;
    Color col0=nl[0].color1;
    pov.write_sphere(p0,rad0,col0,name);
    for (unsigned int i=1; i<nl.size();++i) {
      geom::Vec3 p1=nl[i].atom.GetPos();
      float rad1=0.2;
      Color col1=nl[i].color1;
      pov.write_sphere(p1,rad1,col1,name);
      geom::Vec3 pm=(p0+p1)*0.5;
      pov.write_cyl(p0,pm,rad0,col0,name,true);
      pov.write_cyl(pm,p1,rad1,col1,name,true);
      p0=p1;
      rad0=rad1;
      col0=col1;
    }
  }

  pov.inc() << " }\n";
}

bool TraceRenderer::CanSetOptions(RenderOptionsPtr& render_options)
{
  return render_options.get()->GetRenderMode()==RenderMode::TRACE;
}

void TraceRenderer::SetOptions(RenderOptionsPtr& render_options)
{
  options_=boost::static_pointer_cast<TraceRenderOptions>(render_options);
}

RenderOptionsPtr TraceRenderer::GetOptions()
{
 return options_;
}

TraceRenderer::~TraceRenderer() 
{

}

}}}
