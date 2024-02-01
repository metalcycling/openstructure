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
#ifndef OST_GFX_IMPL_TRACE_RENDERER_HH
#define OST_GFX_IMPL_TRACE_RENDERER_HH

/*
   Authors: Marco Biasini, Ansgar Philippsen
*/

#include <ost/gfx/impl/backbone_trace.hh>
#include <ost/gfx/impl/entity_detail.hh>
#include <ost/gfx/impl/trace_renderer_base.hh>

#include <ost/gfx/entity_fw.hh>
#include <ost/gfx/render_options/trace_render_options.hh>

namespace ost { namespace gfx { namespace impl {

using namespace impl;

/// \brief render backbone trace of entity
/// \internal
class DLLEXPORT_OST_GFX TraceRenderer: public TraceRendererBase {
public:
  TraceRenderer(BackboneTrace* trace);

  virtual void PrepareRendering();
  virtual void PrepareRendering(BackboneTrace& trace_subset, 
                                IndexedVertexArray& va, bool is_sel);

  virtual void RenderPov(PovState& pov, const std::string& name);

  virtual bool CanSetOptions(RenderOptionsPtr& render_options);
  virtual void SetOptions(RenderOptionsPtr& render_options);
  virtual RenderOptionsPtr GetOptions();

  virtual ~TraceRenderer();
private:
  TraceRenderOptionsPtr  options_;
};

}}}

#endif // OST_GFX_TRACE_RENDERER_HH
