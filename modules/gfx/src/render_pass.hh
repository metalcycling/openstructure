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
#ifndef OST_GFX_RENDER_PASS_HH
#define OST_GFX_RENDER_PASS_HH

/*
   Author: Marco Biasini
 */
namespace ost { namespace gfx {

typedef enum {
  /// \brief shaded and textured renderpass, intended to be shown directly 
  ///    to the user
  STANDARD_RENDER_PASS,
  /// \brief only renders "glowing" parts objects
  GLOW_RENDER_PASS,
  /// \brief renders objects with solid color as used for shadow map
  ///     generation
  DEPTH_RENDER_PASS,
  /// \ brief for transparent objects
  /// the standard render pass is not supposed to render transparent
  /// objects, as they should get drawn after the rest of the scene
  TRANSPARENT_RENDER_PASS,
  OVERLAY_RENDER_PASS
} RenderPass;
 
}}

#endif
