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
#ifndef OST_GFX_MAP_SLAB_HH
#define OST_GFX_MAP_SLAB_HH

/*
  Author: Ansgar Philippsen
*/

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

#include <ost/geom/geom.hh>

#include <ost/img/map.hh>

#include "gfx_object.hh"
#include "glext_include.hh"

namespace ost { namespace gfx {

class MapSlab;
typedef boost::shared_ptr<MapSlab> MapSlabP;

/// \brief Render slice of \ref img::ImageHandle "3d image"
/// 
/// A slab is a planar slice through a \ref img::ImageHandle "3D map". Its 
/// orientation is defined by a \ref geom::Plane "plane". The value range and 
/// coloring of the slab can be adjusted with the ColorBy() methods.
class DLLEXPORT_OST_GFX MapSlab: public GfxObj {
  /// \internal
  struct PolyListEntry{
    PolyListEntry() {}
    PolyListEntry(const geom::Vec3& pp, const geom::Vec2& tt): p(pp), t(tt) {}
    geom::Vec3 p;
    geom::Vec3 t;
  };
  typedef std::vector<PolyListEntry> PolyList;

public:
  /// \brief create new slab
  /// 
  /// \param name is the name of the object
  /// \param mh is a 3D image
  /// \param p is the plane that is intersected with \a mh
  MapSlab(const String& name, const img::MapHandle& mh, 
          const geom::Plane& p);

  virtual geom::AlignedCuboid GetBoundingBox(bool use_tf=true) const;

  virtual geom::Vec3 GetCenter() const;
  virtual void CustomRenderGL(RenderPass pass);
  
  /// \brief get plane of slab
  virtual void CustomRenderPov(PovState& pov);

  geom::Plane GetPlane() const;
  /// \brief set slab plane
  /// 
  /// Will force an update of the slab.
  void SetPlane(const geom::Plane& p);
  
  /// \name adjust coloring and value range
  //@{
  /// \brief color by gradient, using explicit min/max values
  void ColorBy(const Gradient& g, float minv, float maxv);
  /// \brief color by gradient, using min/max of slab
  void ColorBy(const Gradient& g); 
  /// \brief color by interpolating between two colors with explicit min 
  ///     and max values
  void ColorBy(const Color& c1, const Color& c2, float minv, float maxv);
  /// \brief color by interpolating between two colors with min/max of slab
  void ColorBy(const Color& c1, const Color& c2);
  //@}
protected:
  virtual void CustomPreRenderGL(bool flag);

  void Rebuild();

private:
  img::MapHandle mh_;
  geom::Plane plane_;
  uint usize_;
  uint vsize_;
  GLuint tex_id_;
  boost::shared_array<float> tex_data_;
  bool rebuild_;

  bool valid_;
  geom::Vec3 c0_,c1_,c2_,c3_;
  PolyList poly_list_;

  float minv_,maxv_;
  Gradient grad_;

  bool calc_poly();
  void rebuild_tex();
  void render_slab();
};

}}

#endif
