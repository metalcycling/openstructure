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
#ifndef OST_SURFACE_IMPL_HH
#define OST_SURFACE_IMPL_HH

#include <vector>

#include <ost/mol/surface_prop.hh>
#include <ost/mol/spatial_organizer.hh>

#include "surface_impl_fw.hh"

namespace ost { namespace mol { namespace impl {

/// \internal
typedef SpatialOrganizer<SurfaceVertexID> SpatialSurfaceVertexOrganizer;

/// \internal
class DLLEXPORT_OST_MOL SurfaceImpl {

  typedef std::vector<SurfaceVertex> VertexList;
  typedef std::vector<SurfaceTri> FaceList;

public:
  SurfaceImpl();

  SurfaceVertexID AddVertex(const SurfaceVertex& v);
  SurfaceTriID AddTri(SurfaceVertexID v0, SurfaceVertexID v1, SurfaceVertexID v2);

  void Attach(const EntityHandle& eh, Real cutoff);
  void Attach(const EntityView& ev, Real cutoff);

  SurfaceVertex& Vertex(SurfaceVertexID id);
  const SurfaceVertex& Vertex(SurfaceVertexID id) const;

  SurfaceTri& Tri(SurfaceTriID id);
  const SurfaceTri& Tri(SurfaceTriID id) const;

  std::vector<SurfaceVertexID> GetVertexIDList() const;
  std::vector<SurfaceVertexID> GetTriIDList() const;

  SurfaceVertexList FindWithin(const geom::Vec3& pos, Real distance);

  void Invert();
  
private:
  VertexList vertex_list_;
  FaceList face_list_;

  SpatialSurfaceVertexOrganizer vertex_organizer_;

  std::vector<SurfaceVertexID> vertex_id_list_;
  std::vector<SurfaceTriID> tri_id_list_;
};

}}} // ns

#endif
