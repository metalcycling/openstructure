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
#ifndef OST_GFX_IMPL_MAPPED_PROPERTY_HH
#define OST_GFX_IMPL_MAPPED_PROPERTY_HH

/*
  Author: Marco Biasini
 */

#include <ost/mol/entity_view.hh>
#include <ost/img/map.hh>

#include <ost/gfx/module_config.hh>
#include <ost/gfx/gradient.hh>

namespace ost { namespace gfx { namespace impl {

inline float Normalize(float v, float min_v, float max_v)
{
  return (v-min_v)/(max_v-min_v);
}

inline float Clamp(float v, float min_v, float max_v)
{
  return std::max(min_v,std::min(max_v,v));
}
  
float DLLEXPORT_OST_GFX MappedProperty(const mol::EntityView& ev, 
                                       const String& prop,
                                       const geom::Vec3& pos);
Color DLLEXPORT_OST_GFX MappedProperty(const mol::EntityView& ev, 
                                       const String& prop,
                                       const Gradient& g, float minv, float maxv,
                                       const geom::Vec3& pos);                      
float DLLEXPORT_OST_GFX MappedProperty(const img::MapHandle& mh, 
                                       const String& prop,
                                       const geom::Vec3& pos);
Color DLLEXPORT_OST_GFX MappedProperty(const img::MapHandle& mh, 
                                       const String& prop,
                                       const Gradient& g, float minv, float maxv,
                                       const geom::Vec3& pos);

}}}

#endif
