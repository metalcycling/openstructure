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
#ifndef OST_GFX_SURFACE_HH
#define OST_GFX_SURFACE_HH

/*
  graphical representation of all kinds of surfaces

  Author: Ansgar Philippsen
*/

#include <boost/shared_ptr.hpp>

#include <ost/config.hh>
#include <ost/mol/surface.hh>
#include "gfx_object.hh"
#include "vertex_array.hh"

#include <ost/gfx/color_ops/uniform_color_op.hh>
#include <ost/gfx/color_ops/basic_gradient_color_op.hh>
#include <ost/gfx/color_ops/gradient_level_color_op.hh>
#include <ost/gfx/color_ops/entity_view_color_op.hh>
#include <ost/gfx/color_ops/map_handle_color_op.hh>

namespace ost { namespace gfx {

class Surface;
typedef boost::shared_ptr<Surface> SurfaceP;

class DLLEXPORT_OST_GFX Surface: public GfxObj {
  typedef std::map<mol::SurfaceVertexID, VertexID> VMap;

public:
  Surface(const String& name, const mol::SurfaceHandle& sh);
  
  virtual void OnRenderModeChange();

  virtual void CustomRenderGL(RenderPass pass);
  virtual void CustomRenderPov(PovState& pov);
  virtual geom::AlignedCuboid GetBoundingBox(bool use_global=false) const;

  mol::SurfaceHandle GetHandle() const;

  void Rebuild();
  void Replace(const mol::SurfaceHandle& sh);

  /// \brief set color for selection
  void SetColor(const Color& col, const String& selection=String(""));

  // GfxObj property interface
  virtual void ColorBy(const mol::EntityView& ev, 
                       const String& prop,
                       const Gradient& g, float minv, float maxv);                    
  // GfxObj property interface
  virtual void ColorBy(const img::MapHandle& mh,
                       const String& prop,
                       const Gradient& g,float minv, float maxv);

  // map property to color gradient from minv to maxv
  void ColorBy(const String& prop, 
               const Gradient& gradient,
               float minv,float maxv,
               mol::Prop::Level hint=mol::Prop::UNSPECIFIED);

  // convenience
  void ColorBy(const String& prop, 
               const Gradient& gradient,
               mol::Prop::Level hint=mol::Prop::UNSPECIFIED);

  // convenience
  void ColorBy(const String& prop, 
               const Color& c1, const Color& c2, 
               float min, float max,
               mol::Prop::Level hint=mol::Prop::UNSPECIFIED);

  // convenience
  void ColorBy(const String& prop, 
               const Color& c1, const Color& c2,
               mol::Prop::Level hint=mol::Prop::UNSPECIFIED);


  void Apply(const gfx::UniformColorOp& op, bool store=true);
  void Apply(const gfx::BasicGradientColorOp& op, bool store=true);
  void Apply(const gfx::GradientLevelColorOp& op, bool store=true);
  void Apply(const gfx::EntityViewColorOp& op, bool store=true);
  void Apply(const gfx::MapHandleColorOp& op, bool store=true);

  void CleanColorOps();
  void ReapplyColorOps();

protected:
  virtual void CustomPreRenderGL(bool flag);

private:

  mol::SurfaceHandle sh_;
  VMap vmap_;
  
  mutable bool recalc_bb_;

  boost::ptr_vector<gfx::ColorOp> c_ops_;
};

}} // ns

#endif
