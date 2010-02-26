//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2010 by the OpenStructure authors
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

#include "line_trace_renderer.hh"
#include <ost/gfx/entity.hh>

namespace ost { namespace gfx { namespace impl {

LineTraceRenderer::LineTraceRenderer(BackboneTrace& trace): 
  TraceRendererBase(trace, 1), options_(new LineTraceRenderOptions())
{}

void LineTraceRenderer::PrepareRendering() 
{
  TraceRendererBase::PrepareRendering();
  va_.Clear();
  this->PrepareRendering(trace_subset_, va_, false);
  sel_va_.Clear();
  sel_va_.SetOutlineWidth(options_->GetLineWidth()+3.0);
  if (this->HasSelection()) {
    this->PrepareRendering(sel_subset_, sel_va_, true);
    sel_va_.SetLineWidth(options_->GetLineWidth()+4.0);    
  }
}

void LineTraceRenderer::PrepareRendering(TraceSubset& trace_subset,
                                         IndexedVertexArray& va, bool is_sel)
{

  const Color& sel_clr=this->GetSelectionColor();
  if(options_!=NULL){
    va.Clear();
    va.SetLighting(false);
    va.SetCullFace(false);
    va.SetColorMaterial(false);
    va.SetMode(0x2);
    va.SetTwoSided(true);
    va.SetLineWidth(options_->GetLineWidth());
    va.SetPointSize(options_->GetLineWidth());
    va.SetAALines(options_->GetAALines());
    for (int node_list=0; node_list<trace_subset.GetSize(); ++node_list) {
      const NodeListSubset& nl=trace_subset[node_list];

      if (nl.GetSize()==2) {
        VertexID p0, p1;
        if (nl.AtStart()==0) {
          p0=va.Add(nl[0].atom.GetPos(), geom::Vec3(),
                    is_sel ? sel_clr : nl[0].color1);
          p1=va.Add((nl[0].atom.GetPos()+nl[1].atom.GetPos())/2, 
                    geom::Vec3(), is_sel ? sel_clr : nl[1].color1);
        } else {
          p0=va.Add((nl[0].atom.GetPos()+nl[1].atom.GetPos())/2, 
                    geom::Vec3(), is_sel ? sel_clr : nl[0].color1);
          p1=va.Add(nl[1].atom.GetPos(), 
                    geom::Vec3(), is_sel ? sel_clr : nl[1].color1);
        }
        va.AddLine(p0, p1);        
        continue;
      }
      if (nl.GetSize()<3) {
        continue;
      }
      VertexID p0;      
      if (nl.AtStart()==0) {
        p0=va.Add(nl[0].atom.GetPos(), geom::Vec3(),
                  is_sel ? sel_clr : nl[0].color1);
      } else {
        p0=va.Add((nl[0].atom.GetPos()+nl[1].atom.GetPos())/2, 
                  geom::Vec3(), is_sel ? sel_clr : nl[0].color1);
      }
      for (int i=1; i<nl.GetSize()-1;++i) {
        const NodeEntry& entry=nl[i];
        VertexID p1 =va.Add(entry.atom.GetPos(), geom::Vec3(), 
                            is_sel ? sel_clr : entry.color1);
        va.AddLine(p0, p1);
        p0=p1;
      }
      const NodeEntry& entry=nl[nl.GetSize()-1];      
      if (nl.AtEnd()==0) {
        VertexID p1 =va.Add(entry.atom.GetPos(), geom::Vec3(), 
                            is_sel ? sel_clr : entry.color1);
        va.AddLine(p0, p1);                            
      } else {
        geom::Vec3 p=(entry.atom.GetPos()+nl[nl.GetSize()-2].atom.GetPos())*0.5;
        VertexID p1 =va.Add(p, geom::Vec3(), 
                            is_sel ? sel_clr : entry.color1);
        va.AddLine(p0, p1);
      }
    }
  }
  sel_state_=0;
  state_=0;
}

void LineTraceRenderer::Render()
{

}

bool LineTraceRenderer::CanSetOptions(RenderOptionsPtr& render_options)
{
 return render_options.get()->GetRenderMode()==RenderMode::LINE_TRACE;
}

void LineTraceRenderer::SetOptions(RenderOptionsPtr& render_options)
{
  options_=boost::static_pointer_cast<LineTraceRenderOptions>(render_options);
}

RenderOptionsPtr LineTraceRenderer::GetOptions()
{
 return options_;
}

LineTraceRenderer::~LineTraceRenderer() 
{

}

}}}
