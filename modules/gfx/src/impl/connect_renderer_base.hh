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
#ifndef OST_GFX_IMPL_CONNECT_RENDERER_HH
#define OST_GFX_IMPL_CONNECT_RENDERER_HH

#include <ost/gfx/module_config.hh>
#include <ost/gfx/impl/entity_detail.hh>
#include <ost/gfx/impl/entity_renderer.hh>
/*
   Author: Marco Biasini, Ansgar Philippsen
 */
namespace ost { namespace gfx { namespace impl {

/// \internal
/// \brief base-class for entity renderers at full connectivity detail like the
///       \ref "simple" SimpleRenderer and \ref "balls & sticks" CustomRenderer
///       renderers.
class DLLEXPORT_OST_GFX ConnectRendererBase : public EntityRenderer {
public:
  ConnectRendererBase();
  
  virtual void PrepareRendering();
  
  virtual geom::AlignedCuboid GetBoundingBox() const;
  
  virtual void Apply(const gfx::ByElementColorOp& op);
  virtual void Apply(const gfx::ByChainColorOp& op);
  virtual void Apply(const gfx::UniformColorOp& op);
  virtual void Apply(const gfx::GradientLevelColorOp& op);
  virtual void Apply(const gfx::EntityViewColorOp& op);

  virtual void PickAtom(const geom::Line3& line, Real line_width,
                        mol::AtomHandle& picked_atom);
  virtual void Apply(const gfx::MapHandleColorOp& op);

  virtual void PickBond(const geom::Line3& line, Real line_width,
                        mol::BondHandle& picked_bond);
  virtual void UpdateViews();  
protected:
  void SetFixedPickRadius(float radius) { pick_radius_=radius;}    
  bool HasFixedPickRadius() const { return pick_radius_!=0.0; }
  float GetFixedPickRadius() const { return pick_radius_; }

protected:
  float     pick_radius_;
  GfxView  view_;
  GfxView  sel_view_;
};
 
}}}

#endif
