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

#include <ost/mol/mol.hh>

#include "backbone_trace.hh"

/*
  Authors: Ansgar Philippsen, Marco Biasini
 */
 
namespace ost { namespace gfx { namespace impl {

class TraceBuilder: public mol::EntityVisitor {
public:
  TraceBuilder(BackboneTrace* bb_trace):
    backbone_trace_(bb_trace),
    last_residue_(),
    list_()
  {}

  virtual bool VisitChain(const mol::ChainHandle& chain)
  {
    if(last_chain_ && chain!=last_chain_) {
      if(!list_.empty()) {
        backbone_trace_->AddNodeList(list_);
        list_.clear();
      }
      last_chain_=chain;
    }
    return true;
  }

  virtual bool VisitResidue(const mol::ResidueHandle& res)
  {
    if(!mol::InSequence(last_residue_,res)) {
      if(!list_.empty()) {
        backbone_trace_->AddNodeList(list_);
        list_.clear();
      }
    }
    // determine atom to add to list
    mol::AtomHandle ca = res.GetCentralAtom();
    if (ca) {
      NodeEntry entry={ca, GfxObj::Ele2Color(ca.GetProp().element),
                       GfxObj::Ele2Color(ca.GetProp().element),
                       geom::Vec3(0.0,0.0,0.0), // this will be set by the gfx trace obj
                       res.GetCentralNormal(),
                       1.0,
                       geom::Vec3(),geom::Vec3(),geom::Vec3(),
                       false};
      list_.push_back(entry);
    }

    last_residue_=res;

    return false;
  }

  virtual void OnExit()
  {
    if (!list_.empty()) {
      backbone_trace_->AddNodeList(list_);
      list_.clear();
    }
  }
private:
  BackboneTrace*     backbone_trace_;
  mol::ResidueHandle last_residue_;
  mol::ChainHandle   last_chain_;
  NodeEntryList      list_;
};

BackboneTrace::BackboneTrace(const mol::EntityView& ent)
{
  view_=ent;
  TraceBuilder trace(this);
  if (view_) {
    view_.Apply(trace);
  }
}

BackboneTrace::BackboneTrace()
{
  
}

void BackboneTrace::SetView(const mol::EntityView& ent)
{
  view_=ent;
  if (view_) {
    node_list_list_.clear();
    TraceBuilder trace(this);    
    view_.Apply(trace);
  }
}

void BackboneTrace::AddNodeList(const NodeEntryList& l)
{
  if(l.size()>=3) {
    node_list_list_.push_back(l);
    // assign direction and normal vectors for each entry
    // they are composed of the i-1 -> i and i->i+1 directions
    //
    // this same algorithm is used in the spline generation, so
    // perhaps all of this here is not necessary ?!
    //
    NodeEntry* e0=&node_list_list_.back()[0];
    NodeEntry* e1=&node_list_list_.back()[1];
    NodeEntry* e2=&node_list_list_.back()[2];
    geom::Vec3 p0 = e0->atom.GetPos();
    geom::Vec3 p1 = e1->atom.GetPos();
    geom::Vec3 p2 = e2->atom.GetPos();

    e0->direction=geom::Normalize(p1-p0);
    // e0->normal is set afterwards to normal of second one
    // backup residue normal
    e0->v1 = e0->normal;

    //reference normal to avoid twisting
    geom::Vec3 nref=geom::Normalize(geom::Cross(p0-p1,p2-p1));

    // start loop with the second
    unsigned int i=1;
    for(;i<node_list_list_.back().size()-1;++i) {
      geom::Vec3 p10 = p0-p1;
      geom::Vec3 p12 = p2-p1;
      if(p10==-p12 || p10==p12) p12+=geom::Vec3(0.001,0.001,0.001);
      e1->v1=e1->normal;
      // twist avoidance
      if(geom::Dot(e0->v1,e1->v1)<0.0) {
        e1->v1=-e1->v1;
      }
      e1->normal=geom::Normalize(geom::Cross(p10,p12));
      float omega=0.5*acos(geom::Dot(geom::Normalize(p10),geom::Normalize(p12)));
      geom::Vec3 orth=geom::AxisRotation(e1->normal, -omega)*p12;
      e1->direction=geom::Normalize(geom::Cross(e1->normal,orth));

      // align normals to avoid twisting
      //if(geom::Dot(e1->normal,nref)<0.0) e1->normal=-e1->normal;
      //nref=e1->normal;
      // skip over shift for the last iteration
      if(i==node_list_list_.back().size()-2) break;
      // shift to i+1 for next iteration
      e0=&node_list_list_.back()[i];
      e1=&node_list_list_.back()[i+1];
      e2=&node_list_list_.back()[i+2];
      p0 = e0->atom.GetPos();
      p1 = e1->atom.GetPos();
      p2 = e2->atom.GetPos();
    }
    // finish remaining values
    // i is at size-2 due to break statement above
    node_list_list_.back()[0].normal=node_list_list_.back()[1].normal;
    node_list_list_.back()[i+1].direction=geom::Normalize(p2-p1);
    node_list_list_.back()[i+1].v1=node_list_list_.back()[i+1].normal;
    if (geom::Dot(node_list_list_.back()[i].v1,
                  node_list_list_.back()[i+1].v1)<0.0) {
      node_list_list_.back()[i+1].v1=-node_list_list_.back()[i+1].v1;
    }
    node_list_list_.back()[i+1].normal=node_list_list_.back()[i].normal;
  }
}

int BackboneTrace::GetListCount() const
{
  return node_list_list_.size();
}

const NodeEntryList& BackboneTrace::GetList(int index) const
{
  return node_list_list_[index];
}

NodeEntryList& BackboneTrace::GetList(int index)
{
  return node_list_list_[index];
}

void BackboneTrace::Rebuild()
{
  TraceBuilder trace(this);
  view_.Apply(trace);  
}

NodeListSubset::NodeListSubset(BackboneTrace& trace, int index):
  trace_(trace), list_index_(index), at_start_(0), at_end_(0)
{
}

TraceSubset::TraceSubset(BackboneTrace& trace, 
                         const mol::EntityView& view, int n):
  trace_(trace), overshoot_(n)
{
  this->Update(view);
}

TraceSubset::TraceSubset(BackboneTrace& trace, int n):
  trace_(trace), overshoot_(n)
{
  
}


void TraceSubset::NodeListEnd(NodeListSubset& nl, int c, int s)
{
  if (nl.GetSize()>0) {
    int n=std::min(s-(nl.indices_.back()+1), overshoot_);
    nl.at_end_=n;    
    while (n>0) {
      nl.indices_.push_back(nl.indices_.back()+1);
      --n;
    }
  }
}

void TraceSubset::NodeListStart(NodeListSubset& nl, int c)
{
  if (nl.GetSize()==0 && c>0) {
    for (int i=std::max(0,c-overshoot_); i<c; ++i) {
      nl.indices_.push_back(i);
    }
    nl.at_start_=nl.indices_.size();
  }
}


void TraceSubset::Update(const mol::EntityView& view)
{

  lists_.clear();
  if (!view.IsValid()) {
    return;
  }
  for (int i=0; i<trace_.GetListCount(); ++i) {
    const NodeEntryList& l=trace_.GetList(i); 
    int c=0;
    NodeEntryList::const_iterator j=l.begin(),e=l.end();  
    NodeListSubset curr(trace_, i);
    while (j!=e) {  
      while (j!=e && view.FindAtom(j->atom)) {
        if (curr.indices_.empty()) {
          this->NodeListStart(curr, c);
        }
        curr.indices_.push_back(c);
        ++j; ++c;
      }
      if (curr.GetSize()) {
        this->NodeListEnd(curr, curr.indices_.back(), l.size());
        lists_.push_back(curr);
        curr=NodeListSubset(trace_, i);
      } else {
        ++c; ++j;        
      }
    }
  }
}

NodeListSubset& NodeListSubset::operator=(const NodeListSubset& rhs)
{
  trace_=rhs.trace_;
  list_index_=rhs.list_index_;
  indices_=rhs.indices_;
  return *this;
}

TraceSubset& TraceSubset::operator=(const TraceSubset& rhs)
{
  trace_=rhs.trace_;
  lists_=rhs.lists_;
  overshoot_=rhs.overshoot_;
  return *this;
}

}}}
