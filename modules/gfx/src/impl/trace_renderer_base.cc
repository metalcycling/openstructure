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

#include <ost/gfx/scene.hh>

#include "trace_renderer_base.hh"


namespace ost { namespace gfx { namespace impl {

namespace {

void set_node_entry_color(NodeEntry& e, ColorMask mask,
                                             const Color& c)
{
  if (mask & MAIN_COLOR) e.color1=c;
  if (mask & DETAIL_COLOR) e.color2=c;
}

template <typename T1>
inline void apply_color_op(TraceRendererBase* rend, BackboneTrace& trace_subset, T1 get_col, const ColorOp& op)
{
  rend->UpdateViews();
  ColorMask mask = op.GetMask();
  mol::EntityView view;
  if(op.IsSelectionOnly()){
    mol::Query q(op.GetSelection());
    for (int node_list=0; node_list<trace_subset.GetListCount(); ++node_list) {
      NodeEntryList& nl=trace_subset.GetList(node_list);
      for (unsigned int i=0; i<nl.size();++i) {
        if(nl[i].atom.IsValid()) {
          if (q.IsAtomSelected(nl[i].atom)) {
            std::pair<bool,Color> clr =get_col.ColorOfAtom(nl[i].atom);
            if(clr.first) {
              set_node_entry_color(nl[i],mask,clr.second);
            }
          }
        }
      }
    }
  } else {
    mol::EntityView view = op.GetView();
    for (int node_list=0; node_list<trace_subset.GetListCount(); ++node_list) {
      NodeEntryList& nl=trace_subset.GetList(node_list);
      for (unsigned int i=0; i<nl.size();++i) {
        if(view.ViewForHandle(nl[i].atom)){
          std::pair<bool,Color> clr =get_col.ColorOfAtom(nl[i].atom);
          if(clr.first) {
            set_node_entry_color(nl[i],mask,clr.second);
          }
        }
      }
    }
  }
};


} //ns

TraceRendererBase::TraceRendererBase(BackboneTrace* trace, int n):
  trace_(trace), trace_subset_(), sel_subset_()
{
}

void TraceRendererBase::PrepareRendering(bool twist_hack)
{
  trace_->SetTwistHack(twist_hack);
  trace_subset_.SetTwistHack(twist_hack);
  if(this->HasSelection()) sel_subset_.SetTwistHack(twist_hack);
  if (state_ & DIRTY_VA) {
    trace_->OnUpdatedPositions();
    trace_subset_.OnUpdatedPositions();
    if(this->HasSelection()) sel_subset_.OnUpdatedPositions();
    // don't clear DIRTY_VA flag - derived classed may depend on it
  }
}

void TraceRendererBase::UpdateViews()
{
  if (state_ & DIRTY_VIEW) {
    mol::EntityView view=this->GetEffectiveView();
    trace_subset_=trace_->CreateSubset(view);
    state_&=~DIRTY_VIEW;    
  }
  if (this->HasSelection() && (sel_state_ & DIRTY_VIEW)) {
    sel_subset_ = trace_->CreateSubset(sel_);
    sel_state_&=~DIRTY_VIEW;    
  }
}

geom::AlignedCuboid TraceRendererBase::GetBoundingBox() const
{
  geom::Vec3 mmin(std::numeric_limits<float>::max(),
		  std::numeric_limits<float>::max(),
		  std::numeric_limits<float>::max());
  geom::Vec3 mmax(-std::numeric_limits<float>::max(),
		  -std::numeric_limits<float>::max(),
		  -std::numeric_limits<float>::max());
		  
  assert(!(state_ & DIRTY_VIEW));
  for (int node_list=0; node_list<trace_subset_.GetListCount(); ++node_list) {
    const NodeEntryList& nl=trace_subset_.GetList(node_list);
    for (unsigned int i=0; i<nl.size();++i) {
      const NodeEntry& entry=nl[i];      
      geom::Vec3 p=entry.atom.GetPos();
      mmin=geom::Min(mmin, p);
      mmax=geom::Max(mmax, p);
    }
  }  
  return geom::AlignedCuboid(mmin, mmax);
}

void TraceRendererBase::Apply(const gfx::ByElementColorOp& op)
{
  apply_color_op(this,trace_subset_,ByElementGetCol(),op);
  state_|=DIRTY_VA;
}

void TraceRendererBase::Apply(const gfx::ByChainColorOp& op)
{
  apply_color_op(this,trace_subset_,ByChainGetCol(op),op);
  state_|=DIRTY_VA;
}


void TraceRendererBase::Apply(const gfx::UniformColorOp& op)
{
  apply_color_op(this,trace_subset_,UniformGetCol(op.GetColor()),op);
  state_|=DIRTY_VA;
}

void TraceRendererBase::Apply(const gfx::GradientLevelColorOp& op)
{
  apply_color_op(this,trace_subset_,GradientLevelGetCol(op),op);
  state_|=DIRTY_VA;
}

void TraceRendererBase::Apply(const gfx::EntityViewColorOp& op)
{
  apply_color_op(this,trace_subset_,EntityViewGetCol(op),op);
  state_|=DIRTY_VA;
}

bool TraceRendererBase::HasDataToRender() const
{
  assert(!(state_ & DIRTY_VIEW));
  return this->trace_subset_.GetListCount()>0;
}

void TraceRendererBase::set_node_colors(const Color& col, const mol::Query& q, 
                                        ColorMask mask)
{
  this->UpdateViews();
  for (int node_list=0; node_list<trace_subset_.GetListCount(); ++node_list) {
    NodeEntryList& nl=trace_subset_.GetList(node_list);
    for (unsigned int i=0; i<nl.size();++i) {
      if(q.IsAtomSelected(nl[i].atom)) {
        set_node_entry_color(nl[i],mask,col);
      }
    }
  }
  state_|=DIRTY_VA;  
}

void TraceRendererBase::PickAtom(const geom::Line3& line, Real line_width,
                                 mol::AtomHandle& picked_atom)
{
  static const float max_dist = 2.0;
  float zmin=1.0;
  if (picked_atom.IsValid()) {
    geom::Vec3 pp=Scene::Instance().Project(picked_atom.GetPos());
    if (pp[2]>0.0 && pp[2]<zmin) {
       zmin=pp[2];
    }    
  }

  mol::AtomHandle atom;
  for (int node_list=0; node_list<trace_subset_.GetListCount(); ++node_list) {
    NodeEntryList& nl=trace_subset_.GetList(node_list);
    for (unsigned int i=0; i<nl.size();++i) {
      geom::Vec3 p=nl[i].atom.GetPos();
      float dist = geom::Distance(line, p);
      if(dist<=max_dist) {
        geom::Vec3 pp = Scene::Instance().Project(p);
        if(pp[2]>0.0 && pp[2]<zmin) {
          zmin=pp[2];
          atom=nl[i].atom.GetHandle();
        }
      }
    }
  }
  if (atom.IsValid()) {
    picked_atom=atom;
  }
}


void TraceRendererBase::PickBond(const geom::Line3& line, Real line_width,
                                 mol::BondHandle& picked_bond)
{
  // don't to anything
}

void TraceRendererBase::Apply(const gfx::MapHandleColorOp& op)
{
  apply_color_op(this,trace_subset_,MapHandleGetCol(op),op);
  state_|=DIRTY_VA;
}

void TraceRendererBase::rebuild_sel(const SplineEntryListList& spline_list_list, 
                                    SplineEntryListList& sel_spline_list_list,
                                    const Color& sel_color)
{
  // extract spline segments from list_list that match 
  // (via id) the selection subset
  // first put all ids into a set for fast lookup
  std::set<int> id_set;
  for(int nlc=0;nlc<sel_subset_.GetListCount();++nlc) {
    const NodeEntryList& nelist=sel_subset_.GetList(nlc);
    for(NodeEntryList::const_iterator nit=nelist.begin();nit!=nelist.end();++nit) {
      id_set.insert(nit->id);
    }
  }
  // now find all matching spline segments
  sel_spline_list_list.clear();
  for(SplineEntryListList::const_iterator sit=spline_list_list.begin();sit!=spline_list_list.end();++sit) {
    const SplineEntryList& slist=*sit;
    SplineEntryList nlist;
    unsigned int sc=0;
    while(sc<slist.size()) {
      int curr_id=slist.at(sc).id;
      if(id_set.count(curr_id)>0) {
        // if a match is found, add all until a new id is found
        while(sc<slist.size() &&  slist.at(sc).id==curr_id) {
          nlist.push_back(slist[sc++]);
          // override with the selection color
          nlist.back().color1=sel_color;
          nlist.back().color2=sel_color;
        }
      } else {
        // introduce break
        if(!nlist.empty()) {
          sel_spline_list_list.push_back(nlist);
          nlist.clear();
        }
        // and advance to the next id
        while(sc<slist.size() &&  slist.at(sc).id==curr_id) ++sc;
      }
    }
    if(!nlist.empty()) {
      sel_spline_list_list.push_back(nlist);
      nlist.clear();
    }
  }
}

}}} // ns


