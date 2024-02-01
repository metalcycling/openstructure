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

#include "sline_renderer.hh"

#include <ost/gfx/entity.hh>

namespace ost { namespace gfx { namespace impl {

SlineRenderer::SlineRenderer(BackboneTrace* trace): 
  TraceRendererBase(trace, 3), 
  options_(new SlineRenderOptions()),
  spline_list_list_(),
  sel_spline_list_list_()
{
  this->SetName("Fast Spline");
}

bool SlineRenderer::CanSetOptions(RenderOptionsPtr& render_options)
{
 return render_options->GetRenderMode()==RenderMode::SLINE;
}

void SlineRenderer::SetOptions(RenderOptionsPtr& render_options)
{
  options_ = boost::static_pointer_cast<SlineRenderOptions>(render_options);
}

RenderOptionsPtr SlineRenderer::GetOptions()
{
 return options_;
}

void SlineRenderer::PrepareRendering()
{
  TraceRendererBase::PrepareRendering(false);
  if(state_>0) {
    prepare_rendering(trace_subset_, va_, spline_list_list_);
    rebuild_spline_obj(va_, spline_list_list_, false);
  }
  if (this->HasSelection() && (state_>0 || sel_state_>0)) {
    Color sel_color=GetSelectionColor();
    rebuild_sel(spline_list_list_,sel_spline_list_list_,sel_color);
    rebuild_spline_obj(sel_va_, sel_spline_list_list_, true);
    va_.SetOpacity(sel_color.Alpha());
  }
  sel_state_=0;
  state_=0;
}

void SlineRenderer::prepare_rendering(const BackboneTrace& subset, 
                                      IndexedVertexArray& va,
                                      SplineEntryListList& spline_list_list)
{
  if(options_==NULL) {
    LOG_DEBUG("SlineRenderer: NULL options, not creating objects");
  }

  int spline_detail=std::max((unsigned int) 1, options_->GetSplineDetail());
  SplineEntryListList tmp_sll;
  for (int node_list=0; node_list<subset.GetListCount(); ++node_list) {
    LOG_DEBUG("SplineRenderer: collecting spline entries for node list " << node_list);
    // first build the spline
    SplineEntryList spl;
    const NodeEntryList& nl=subset.GetList(node_list);
    for (unsigned int i=0; i<nl.size();++i) {
      int type=0;
      const NodeEntry& entry=nl[i];
      SplineEntry ee(entry.atom.GetPos(),entry.direction,
                     entry.normal, entry.rad, 
                     entry.color1, 
                     entry.color2,
                     type, entry.id);
      ee.v1 = entry.v1;
      spl.push_back(ee);
    }
    LOG_DEBUG("SplineRenderer: found " << spl.size() << " entries");
    if(!spl.empty()) {
      tmp_sll.push_back(spl);
    }
  }
  spline_list_list.clear();
  for(SplineEntryListList::const_iterator sit=tmp_sll.begin();sit!=tmp_sll.end();++sit) {
    spline_list_list.push_back(Spline::Generate(*sit,spline_detail));
  }
}

void SlineRenderer::rebuild_spline_obj(IndexedVertexArray& va, const SplineEntryListList& spline_list_list, bool is_sel)
{
  va.Clear();
  va.SetLighting(false);
  va.SetCullFace(false);
  va.SetColorMaterial(false);
  va.SetMode(0x2);
  va.SetTwoSided(true);
  va.SetLineWidth(options_->GetLineWidth() + (is_sel ? 2.0 : 0.0));
  va.SetPointSize(options_->GetLineWidth() + (is_sel ? 2.0 : 0.0));
  va.SetAALines(options_->GetAALines());

  for(SplineEntryListList::const_iterator it=spline_list_list.begin();
      it<spline_list_list.end();++it) {
    const SplineEntryList& slist=*it;

    SplineEntryList::const_iterator sit=slist.begin();
    VertexID p0=va.Add(sit->position, geom::Vec3(),sit->color1);
    for (++sit; sit!=slist.end(); ++sit) {
      VertexID p1 = va.Add(sit->position, geom::Vec3(),sit->color1);
      va.AddLine(p0,p1);
      p0=p1;
    }
  }

}

}}}
