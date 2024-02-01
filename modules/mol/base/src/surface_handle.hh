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
#ifndef OST_SURFACE_HANDLE_HL
#define OST_SURFACE_HANDLE_HL

#include <vector>

#include <ost/mol/surface_prop.hh>
#include <ost/mol/impl/surface_impl_fw.hh>

namespace ost { namespace mol {

class SurfaceHandle;
class EntityView;

//  Create a non-managed surface
DLLEXPORT_OST_MOL SurfaceHandle CreateSurface();

typedef std::vector<SurfaceTriID> SurfaceTriIDList;
typedef std::vector<SurfaceVertexID> SurfaceVertexIDList;

// molecular surface
class DLLEXPORT_OST_MOL SurfaceHandle {
public:
  // creates invalid handle
  SurfaceHandle();
  
  // internally used ctor
  SurfaceHandle(const impl::SurfaceImplP& p);

  // associate each vertex with the closest atom from the given entity,
  // but only if an atom is actually within 'cutoff' distance
  void Attach(const EntityHandle& eh, Real cutoff);

  // dito for view
  void Attach(const EntityView& ev, Real cutoff);

  SurfaceVertexList FindWithin(const geom::Vec3&, Real dist) const;

  // returns unique id
  SurfaceVertexID AddVertex(const SurfaceVertex& v);

  // retrieve vertex based on id
  SurfaceVertex GetVertex(SurfaceVertexID id) const;

  // assigns face based on id
  SurfaceTriID AddTri(SurfaceVertexID v0, SurfaceVertexID v1, SurfaceVertexID v2);

  SurfaceTri GetTri(SurfaceTriID id) const;

  SurfaceVertexIDList GetVertexIDList() const;
  SurfaceTriIDList GetTriIDList() const;
  
  // flip normals
  void Invert();

  bool IsValid() const {return static_cast<bool>(impl_);}

  bool operator==(const SurfaceHandle& ref) const { return impl_==ref.impl_; }

  bool operator!=(const SurfaceHandle& ref) const { return !this->operator==(ref); }

private:
  impl::SurfaceImplP impl_;

};

DLLEXPORT_OST_MOL std::ostream& operator<<(std::ostream& os,
                                           const SurfaceHandle& surf);

typedef std::vector<SurfaceHandle> SurfaceHandleList;

}} // ns

#endif
