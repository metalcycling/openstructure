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
#ifndef OST_GFX_IMPL_TRACE_RENDERER_BASE_HH
#define OST_GFX_IMPL_TRACE_RENDERER_BASE_HH

/*
  Authors: Marco Biasini, Ansgar Philippsen
*/

#include <ost/gfx/module_config.hh>
#include <ost/gfx/impl/entity_renderer.hh>
#include <ost/gfx/impl/backbone_trace.hh>

namespace ost { namespace gfx { namespace impl {

/// \internal
/// \brief base class for trace-based entity renders like 
///     \ref "cartoon renderer" CartoonRenderer, \ref "trace" TraceRenderer,
///     \ref "line trace" LineTraceRenderer and 
///     \ref "smooth trace" SlineRenderer.
/// 
/// All trace-based entity renderers share a  common instance of BackboneTrace
/// (which is held by gfx::Entity); it encapulates a smoothed C-alpha trace.
class DLLEXPORT_OST_GFX TraceRendererBase : public EntityRenderer {
public:
  TraceRendererBase(BackboneTrace* trace, int overshoot);
  
  virtual geom::AlignedCuboid GetBoundingBox() const;  
  
  virtual void PrepareRendering(bool twist_hack);
  virtual bool HasDataToRender() const;
  virtual void UpdateViews();  
  virtual void Apply(const gfx::ByElementColorOp& op);
  virtual void Apply(const gfx::ByChainColorOp& op);
  virtual void Apply(const gfx::UniformColorOp& op);
  virtual void Apply(const gfx::GradientLevelColorOp& op);
  virtual void Apply(const gfx::EntityViewColorOp& op);
  virtual void Apply(const gfx::MapHandleColorOp& op);
  virtual void PickAtom(const geom::Line3& line, Real line_width,
                        mol::AtomHandle& picked_atom);

  virtual void PickBond(const geom::Line3& line, Real line_width,
                        mol::BondHandle& picked_bond);
protected:
  
  void set_node_colors(const Color& c, const mol::Query& q, ColorMask mask);

  void rebuild_sel(const SplineEntryListList& spline_list_list, 
                   SplineEntryListList& sel_spline_list_list,
                   const Color& sel_color);

  BackboneTrace* trace_;
  BackboneTrace trace_subset_;
  BackboneTrace sel_subset_;
};

}}}

#endif
