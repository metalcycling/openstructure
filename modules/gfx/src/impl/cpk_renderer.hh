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
#ifndef OST_GFX_CPK_RENDERER_HH
#define OST_GFX_CPK_RENDERER_HH

#include <ost/geom/vec3.hh>

#include <ost/gfx/impl/entity_renderer.hh>
#include <ost/gfx/entity_fw.hh>
#include <ost/gfx/impl/connect_renderer_base.hh>
#include <ost/gfx/render_options/cpk_render_options.hh>

#if OST_SHADER_SUPPORT_ENABLED
#include "fast_spheres.hh"
#endif

namespace ost { namespace gfx { namespace impl {

/// \internal
class DLLEXPORT_OST_GFX CPKRenderer: public ConnectRendererBase {
public:
  CPKRenderer();

  virtual void PrepareRendering();
  
  virtual void RenderPov(PovState& pov, const std::string& name);

  virtual void Render(RenderPass pass);

  virtual bool CanSetOptions(RenderOptionsPtr& render_options);
  virtual void SetOptions(RenderOptionsPtr& render_options);
  virtual RenderOptionsPtr GetOptions();

private:
  void PrepareRendering(GfxView& view, IndexedVertexArray& va, bool is_sel);

  CPKRenderOptionsPtr options_;
#if OST_SHADER_SUPPORT_ENABLED
  FastSphereRenderer fsr_,sel_fsr_;
#endif
};

}}}

#endif // OST_GFX_CPK_RENDERER_HH
