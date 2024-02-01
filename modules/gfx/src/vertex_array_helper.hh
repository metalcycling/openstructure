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
#ifndef OST_GFX_VERTEX_ARRAY_HELPER_HH
#define OST_GFX_VERTEX_ARRAY_HELPER_HH

#include <vector>
#include <ost/geom/geom.hh>

namespace ost { namespace gfx { 

static const unsigned int VA_SPHERE_MAX_DETAIL=12;
static const unsigned int VA_ICO_SPHERE_MAX_DETAIL=4;
static const unsigned int VA_CYL_MAX_DETAIL=12;

namespace detail {

struct PrebuildSphereEntry {
  std::vector<geom::Vec3> vlist;
  std::vector<unsigned int> ilist;
};

PrebuildSphereEntry& GetPrebuildSphere(unsigned int level);

std::vector<geom::Vec3> GetPrebuildCyl(unsigned int level);

std::vector<geom::Vec3> GetPrebuildIcoSphere(unsigned int level);

}}} //

#endif
