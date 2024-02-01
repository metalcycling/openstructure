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
#ifndef OST_GFX_OFFSCREEN_BUFFER_HH
#define OST_GFX_OFFSCREEN_BUFFER_HH



/*
 THE OFFSCREEN BUFFER IS DEPRECATED AND CURRENTLY NOT USED
 Offscreen rendering for image export is now performed in a QT managed
 offscreen buffer (e.g. gui/src/gl_canvas.hh). Reason for that is maintenance,
 as there are platform dependent buffer implementations.   
*/



/*
  Wraps GL offscreen rendering per platform

  boilerplate header, includes platform dependent stuff
*/

#include <vector>

namespace ost { namespace gfx {

class OffscreenBufferFormat
{
public:
  OffscreenBufferFormat(): cbits(8),abits(8),dbits(12),accum(false),multisample(false),samples(0) {}

  unsigned int cbits,abits,dbits; // color, alpha, depth bits
  bool accum;
  bool multisample;
  unsigned int samples;
};

}} // ns

/* 
   instead of creating an abstract base class and
   making runtime polymorphic classes for each
   platform, we do a bit more typing and copy
   the minimal interface to each of the platform
   specific OffscreenBuffer implementations

   OffscreenBuffer interface:

   OffscreenBuffer(unsigned int width, unsigned int height, const OffscreenBufferFormat& f, bool shared=true);
   bool Resize(unsigned int w, unsigned int h);
   bool MakeActive();
   bool IsActive();
   bool IsValid();
*/

#if defined(__linux__)
#if OST_MESA_SUPPORT_ENABLED
#include "impl/mesa_offscreen_buffer.hh"
#else
#include "impl/glx_offscreen_buffer.hh"
#endif
#elif defined(__APPLE__)
#include "impl/cgl_offscreen_buffer.hh"
#elif defined(_MSC_VER)
#include "impl/wgl_offscreen_buffer.hh"
#else
#error platform not found for offscreen rendering
#endif
 
#ifdef Complex
#  undef Complex
#endif

#endif
