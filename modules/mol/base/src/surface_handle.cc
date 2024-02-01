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
#include <ost/mol/surface_handle.hh>
#include <ost/mol/impl/surface_impl.hh>

namespace ost { namespace mol {

// non-managed version
SurfaceHandle CreateSurface()
{
  return SurfaceHandle(impl::SurfaceImplP(new impl::SurfaceImpl()));
}


SurfaceHandle::SurfaceHandle():
  impl_()
{}

SurfaceHandle::SurfaceHandle(const impl::SurfaceImplP& p):
  impl_(p)
{}

SurfaceVertexList SurfaceHandle::FindWithin(const geom::Vec3& pos, Real dist) const
{
  return impl_->FindWithin(pos, dist);
}

void SurfaceHandle::Attach(const EntityHandle& eh, Real cutoff)
{
  impl_->Attach(eh,cutoff);
}

void SurfaceHandle::Attach(const EntityView& ev, Real cutoff)
{
  impl_->Attach(ev,cutoff);
}

SurfaceVertexID SurfaceHandle::AddVertex(const SurfaceVertex& v)
{
  return impl_->AddVertex(v);
}

SurfaceVertex SurfaceHandle::GetVertex(SurfaceVertexID id) const
{
  return impl_->Vertex(id);
}

SurfaceTriID SurfaceHandle::AddTri(SurfaceVertexID v0,
           SurfaceVertexID v1,
           SurfaceVertexID v2)
{
  return impl_->AddTri(v0,v1,v2);
}

SurfaceTri SurfaceHandle::GetTri(SurfaceTriID id) const
{
  return impl_->Tri(id);
}

std::vector<SurfaceVertexID> SurfaceHandle::GetVertexIDList() const
{
  return impl_->GetVertexIDList();
}

std::vector<SurfaceTriID> SurfaceHandle::GetTriIDList() const
{
  return impl_->GetTriIDList();
}

void SurfaceHandle::Invert()
{
  impl_->Invert();
}

std::ostream& operator<<(std::ostream& os, const SurfaceHandle& surf)
{
  if (surf.IsValid()) {
    os << "valid surface";
  } else {
    os << "invalid surface";
  }
  return os;
}

}} // ns
