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

#include "line_trace_renderer.hh"
#include <ost/gfx/entity.hh>

namespace ost { namespace gfx { namespace impl {

namespace {

void add_atom_and_bond(mol::AtomHandle atom1,mol::AtomHandle atom2, AtomEntryMap& amap, BondEntryList& blist)
{
  AtomEntry ae1(atom1,0.0,0.0,Color(1,1,1));
  amap[atom1.GetHashCode()]=ae1;
  AtomEntry ae2(atom2,0.0,0.0,Color(1,1,1));
  amap[atom2.GetHashCode()]=ae2;
  blist.push_back(BondEntry(mol::BondHandle(),0.0f,
                            &amap[atom1.GetHashCode()],
                            &amap[atom2.GetHashCode()]));
}

}


LineTraceRenderer::LineTraceRenderer(BackboneTrace* trace): 
  TraceRendererBase(trace, 1), 
  options_(new LineTraceRenderOptions()),
  amap_(),
  blist_()
{
  this->SetName("Fast Trace");

  for (int node_list=0; node_list<trace->GetListCount(); ++node_list) {
    const NodeEntryList& nl=trace->GetList(node_list);
    for (unsigned int i=0; i<nl.size()-1;++i) {
      add_atom_and_bond(nl[i].atom,nl[i+1].atom,amap_,blist_);
    }
  }
}

void LineTraceRenderer::PrepareRendering() 
{
  TraceRendererBase::PrepareRendering(false);
  va_.Clear();
  this->PrepareRendering(trace_subset_, va_, false);
  sel_va_.Clear();
  sel_va_.SetOutlineWidth(options_->GetLineWidth()+3.0);
  if (this->HasSelection()) {
    //this->PrepareRendering(sel_subset_, sel_va_, true);
    this->PrepareRendering(trace_subset_, sel_va_, true);
    sel_va_.SetLineWidth(options_->GetLineWidth()+4.0);    
  }
}

void LineTraceRenderer::PrepareRendering(const BackboneTrace& trace_subset,
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
    if(is_sel) {
      for (int node_list=0; node_list<trace_subset.GetListCount(); ++node_list) {
        const NodeEntryList& nl=trace_subset.GetList(node_list);
        if(nl.size()<1) continue;
        for(unsigned int i=0;i<nl.size();++i) {
          const NodeEntry& entry=nl[i];
          if(sel_.ViewForHandle(entry.atom).IsValid()) {
            geom::Vec3 apos = entry.atom.GetPos();
            VertexID p0=va.Add(apos, geom::Vec3(),sel_clr);
            if(i>0) {
              VertexID p1 =va.Add(apos+0.5*(nl[i-1].atom.GetPos()-apos), geom::Vec3(), sel_clr);
              va.AddLine(p0, p1);
            }
            if(i<nl.size()-1) {
              VertexID p1 =va.Add(apos+0.5*(nl[i+1].atom.GetPos()-apos), geom::Vec3(), sel_clr);
              va.AddLine(p0, p1);
            }
          }
        }
      }
    } else {
      for (int node_list=0; node_list<trace_subset.GetListCount(); ++node_list) {
        const NodeEntryList& nl=trace_subset.GetList(node_list);
        
        if (nl.size()<2) continue;
        
        VertexID p0=va.Add(nl[0].atom.GetPos(), geom::Vec3(),
                           nl[0].color1);
        for (unsigned int i=1; i<nl.size()-1;++i) {
          const NodeEntry& entry=nl[i];
          VertexID p1 =va.Add(entry.atom.GetPos(), geom::Vec3(), 
                              entry.color1);
          va.AddLine(p0, p1);
          p0=p1;
#if 0
          VertexID p2 =va.Add(entry.atom.GetPos()+entry.direction, geom::Vec3(), 
                              Color(1,0,0));
          VertexID p3 =va.Add(entry.atom.GetPos()+entry.normal, geom::Vec3(), 
                              Color(0,1,1));
          va.AddLine(p0,p2);
          va.AddLine(p0,p3);
#endif
        }
        const NodeEntry& entry=nl.back();
        VertexID p1 =va.Add(entry.atom.GetPos(), geom::Vec3(), 
                            entry.color1);
        va.AddLine(p0, p1);
      }
    }
  }
  sel_state_=0;
  state_=0;
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
