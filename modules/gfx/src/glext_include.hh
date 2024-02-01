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
#ifndef OST_GFX_GLEXT_INCLUDE_HH
#define OST_GFX_GLEXT_INCLUDE_HH

/*
  wrapper for GLEW on platforms that do not define the extension symbols in 
  gl.h
  
  Author: Marco Biasini
 */

#include <ost/config.hh>

#if OST_SHADER_SUPPORT_ENABLED
#  if OST_MODULE==OST_GFX
#    define GLEW_BUILD
#  endif
#  if !defined(__APPLE__)
#   include <ost/gfx/GL/glew.h>
#  endif
#  if defined(_MSC_VER)
#    include <ost/gfx/GL/wglew.h>
#  endif
#endif

#include <ost/gfx/gl_include.hh>

#if defined(__APPLE__)
// On all macOS version we support, OpenGL 2.0 is available, so it's safe to 
// hardcode the value here...
//same for windows vista and above, XP only has 1.1
#define OST_GL_VERSION_2_0 1
#else
#ifdef GLEW_VERSION_2_0
#define OST_GL_VERSION_2_0 1
#endif
#endif
#endif

