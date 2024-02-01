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
#ifndef OST_MAP_ISO_GEN_HH
#define OST_MAP_ISO_GEN_HH

/*
  Author: Ansgar Philippsen
*/

#include <vector>
#include <ost/geom/geom.hh>
#include <ost/img/map.hh>

#include <ost/gfx/map_iso_prop.hh>


namespace ost { namespace gfx {

namespace map_iso {

// this must match the signature in map_iso_spec.hh
// \internal
typedef void (*AddLF)(IndexedVertexArray& va,
                      unsigned int vertex_id[12]);
  
typedef std::vector<AddLF> AddLFList;

} // map_iso ns

namespace impl {
  
// \internal
class DLLEXPORT_OST_GFX IsosurfaceGenerator {
public:
  // singleton
  static IsosurfaceGenerator& Instance();

  IsosurfaceGenerator();

  void Generate(const img::MapHandle& mh, IndexedVertexArray& va, float level);

  void AddLF(uint code, IndexedVertexArray& va,unsigned int vertex_id[12]); 

private:
  map_iso::AddLFList add_lf_func_;
};

}}}

#endif
