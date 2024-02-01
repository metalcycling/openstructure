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
#include <cassert>

#include <ost/mol/entity_handle.hh>
#include <ost/mol/entity_view.hh>
#include <ost/mol/atom_view.hh>

#include "surface_impl.hh"

namespace ost { namespace mol { namespace impl {

SurfaceImpl::SurfaceImpl():
  vertex_list_(),
  face_list_(),
  vertex_organizer_(5.0)
{}

void SurfaceImpl::Attach(const EntityHandle& eh, Real cutoff)
{
  for(VertexList::iterator vit=vertex_list_.begin();vit!=vertex_list_.end();++vit) {
    AtomHandleList ahl = eh.FindWithin(vit->position,cutoff);
    if(ahl.empty()) continue;
    AtomHandleList::const_iterator ait=ahl.begin();
    Real best_dist2 = geom::Length2(ait->GetPos()-vit->position);
    AtomHandle best_atom = *ait;
    ++ait;
    for(;ait!=ahl.end();++ait) {
      Real dist2 = geom::Length2(ait->GetPos()-vit->position);
      if(dist2<best_dist2) {
        best_dist2=dist2;
        best_atom = *ait;
      }
    }
    vit->atom=best_atom;
  }
}

void SurfaceImpl::Attach(const EntityView& ev, Real cutoff)
{
  for(VertexList::iterator vit=vertex_list_.begin();vit!=vertex_list_.end();++vit) {
    AtomViewList avl = ev.FindWithin(vit->position,cutoff);
    if(avl.empty()) continue;
    AtomViewList::const_iterator ait=avl.begin();
    Real best_dist2 = geom::Length2(ait->GetPos()-vit->position);
    AtomView best_atom = *ait;
    ++ait;
    for(;ait!=avl.end();++ait) {
      Real dist2 = geom::Length2(ait->GetPos()-vit->position);
      if(dist2<best_dist2) {
        best_dist2=dist2;
        best_atom = *ait;
      }
    }
    vit->atom=best_atom.GetHandle();
  }
}

SurfaceVertexList SurfaceImpl::FindWithin(const geom::Vec3& pos, Real dist)
{
  SurfaceVertexList vlist;
  SpatialSurfaceVertexOrganizer::ItemList vplist = vertex_organizer_.FindWithin(pos,dist);
  for(SpatialSurfaceVertexOrganizer::ItemList::const_iterator it=vplist.begin();
      it!=vplist.end();++it) {
    vlist.push_back(vertex_list_[*it]);
  }
  return vlist;
}

SurfaceVertexID SurfaceImpl::AddVertex(const SurfaceVertex& v)
{
  vertex_list_.push_back(v);
  vertex_id_list_.push_back(vertex_list_.size()-1);
  vertex_organizer_.Add(vertex_id_list_.back(), v.position);
  return vertex_id_list_.back();
}

SurfaceTriID SurfaceImpl::AddTri(SurfaceVertexID v0, SurfaceVertexID v1, SurfaceVertexID v2)
{
  SurfaceTri f={v0,v1,v2};
  face_list_.push_back(f);
  tri_id_list_.push_back(face_list_.size()-1);
  return tri_id_list_.back();
}

SurfaceVertex& SurfaceImpl::Vertex(SurfaceVertexID id)
{
  return vertex_list_.at(id);
}

const SurfaceVertex& SurfaceImpl::Vertex(SurfaceVertexID id) const
{
  return vertex_list_.at(id);
}

SurfaceTri& SurfaceImpl::Tri(SurfaceTriID id)
{
  return face_list_.at(id);
}

const SurfaceTri& SurfaceImpl::Tri(SurfaceTriID id) const
{
  return face_list_.at(id);
}

std::vector<SurfaceVertexID> SurfaceImpl::GetVertexIDList() const
{
  return vertex_id_list_;
}

std::vector<SurfaceTriID> SurfaceImpl::GetTriIDList() const
{
  return tri_id_list_;
}

void SurfaceImpl::Invert()
{
  for(VertexList::iterator vit=vertex_list_.begin();vit!=vertex_list_.end();++vit) {
    vit->normal=-vit->normal;
  }
  
}

}}} // ns
