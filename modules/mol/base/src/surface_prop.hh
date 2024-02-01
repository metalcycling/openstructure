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
#ifndef OST_SURFACE_PROP_HH
#define OST_SURFACE_PROP_HH

#include <vector>
#include <boost/shared_ptr.hpp>

#include <ost/geom/geom.hh>

#include "atom_handle.hh"

namespace ost { namespace mol {

// molecular surface vertex, POD
struct SurfaceVertex {
  SurfaceVertex() {}
  SurfaceVertex(const geom::Vec3& pos,
                const geom::Vec3& norm,
                int t,
                const AtomHandle& a=AtomHandle()):
    position(pos),
    normal(norm),
    type(t),
    atom(a)
  {}

  bool operator==(const SurfaceVertex& rhs) const {
    return position==rhs.position && normal==rhs.normal && atom==rhs.atom;
  }
  bool operator!=(const SurfaceVertex& rhs) const {
    return !this->operator==(rhs);
  }
  geom::Vec3 position;
  geom::Vec3 normal;
  int type;
  AtomHandle atom; // attached atom
};

typedef unsigned int SurfaceVertexID;

// molecular surface face, POD
struct SurfaceTri {
  SurfaceVertexID v0,v1,v2;
};

typedef unsigned int SurfaceTriID;

typedef std::vector<SurfaceVertex> SurfaceVertexList;

}}

#endif
