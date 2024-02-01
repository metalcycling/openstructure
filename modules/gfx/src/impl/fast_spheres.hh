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
#ifndef OST_GFX_IMPL_FAST_SPHERES_HH
#define OST_GFX_IMPL_FAST_SPHERES_HH

#include <ost/geom/vec3.hh>
#include <ost/gfx/color.hh>

/*
  Author: Ansgar Philippsen
*/

#include <vector>

namespace ost { namespace gfx { namespace impl {

  class DLLEXPORT_OST_GFX FastSphereRenderer {
  public:
    struct VData
    {
      float pos[3];
      float col[4];
      float rad;
    };

    typedef std::vector<VData> DataList;

    FastSphereRenderer(size_t reserve=0);
    FastSphereRenderer(const DataList&);

    void Clear();
    void Add(float pos[3], float col[4], float rad);
    void Add(const geom::Vec3& pos, const Color& col, float rad);

    void RenderGL();

    const DataList& Data() const {return data_;}

  private:
    DataList data_;
  };
  
}}} // ns

#endif
