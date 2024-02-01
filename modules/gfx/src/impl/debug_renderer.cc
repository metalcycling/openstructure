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
#include "debug_renderer.hh"

#include <ost/mol/residue_handle.hh>

#include <ost/gfx/entity.hh>

#include <ost/log.hh>
#include <ost/dyn_cast.hh>


namespace ost {

namespace gfx {

using namespace impl;
using namespace mol;

DebugRenderer::DebugRenderer(BackboneTrace* trace): 
  TraceRendererBase(trace, 2), options_(new SlineRenderOptions()) {
  this->SetName("Debug");
}

void DebugRenderer::PrepareRendering() 
{
#if 0  
  va.Clear();
  va.SetLighting(false);
  va.SetCullFace(false);
  va.SetColorMaterial(false);
  va.SetMode(0x2);
  va.SetTwoSided(true);
  va.SetLineWidth(2.0);

  for(NodeEntryListList::const_iterator ll_it=node_list_list_->begin();ll_it!=node_list_list_->end();++ll_it) {

    SplineEntryList spl;
    for(NodeEntryList::const_iterator it=(*ll_it).begin();it!=(*ll_it).end();++it) {
      int type=0;
      ResidueHandle resh = it->atom.GetResidue();
      SecStructure sst=resh.GetSecStructure();
      if(sst.IsHelical()) {
        type=1;
      } else if(sst.IsExtended()) {
        type=2;
      }

      SplineEntry ee(it->atom.GetPos(),it->direction,it->normal,
                     it->rad,it->color1,it->color2,type,it->id);
      ee.v1 = it->v1;
      spl.push_back(ee);
    }

    SplineEntryList sel = Spline::Generate(spl,std::max((unsigned int) 1,options_->GetSplineDetail()));

    SplineEntryList::const_iterator sit = sel.begin();
    geom::Vec3 ap = sit->position;
    VertexID p0 = va.Add(ap, geom::Vec3(),Color());
    VertexID v0 = va.Add(ap+sit->direction,geom::Vec3(),Color(1.0,0.0,0.0));
    va.AddLine(p0,v0);
    v0 = va.Add(ap+sit->normal,geom::Vec3(),Color(0.0,1.0,1.0));
    va.AddLine(p0,v0);
    v0 = va.Add(ap+sit->v0,geom::Vec3(),Color(1.0,0.0,1.0));
    va.AddLine(p0,v0);
    v0 = va.Add(ap+sit->v1,geom::Vec3(),Color(1.0,1.0,0.0));
    va.AddLine(p0,v0);
    ++sit;
    for(;sit!=sel.end();++sit) {
      ap = sit->position;
      VertexID p1 = va.Add(ap, geom::Vec3(),Color());
      va.AddLine(p0,p1);
      VertexID v1 = va.Add(ap+sit->direction,geom::Vec3(),Color(1.0,0.0,0.0));
      va.AddLine(p1,v1);
      //v1 = va.Add(ap+sit->normal,geom::Vec3(),Color(0.0,1.0,1.0));
      //va.AddLine(p1,v1);
      //v1 = va.Add(ap+sit->v0,geom::Vec3(),Color(1.0,0.0,1.0));
      //va.AddLine(p1,v1);
      //v1 = va.Add(ap+sit->v1,geom::Vec3(),Color(1.0,1.0,0.0));
      //va.AddLine(p1,v1);
      geom::Vec3 orth = -geom::Normalize(geom::Cross(sit->normal,sit->direction));
      if(sit->nflip) orth=-orth;
      v1 = va.Add(ap+orth,geom::Vec3(),Color(0.0,1.0,0.0));
      va.AddLine(p1,v1);
      p0=p1;
    }
  }
#endif
}

bool DebugRenderer::CanSetOptions(RenderOptionsPtr& render_options){
 return render_options.get()->GetRenderMode()==RenderMode::SLINE;
}

void DebugRenderer::SetOptions(RenderOptionsPtr& render_options){
  options_ = boost::static_pointer_cast<SlineRenderOptions>(render_options);
}

RenderOptionsPtr DebugRenderer::GetOptions()
{
 return options_;
}

int DebugRenderer::GetMode(){
  return 1;
}

DebugRenderer::~DebugRenderer() {

}

}

}
