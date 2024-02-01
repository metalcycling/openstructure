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

/*
  Author: Marco Biasini
*/

#include <ost/log.hh>
#include <ost/gfx/offscreen_buffer.hh>
#include <ost/gfx/glext_include.hh>
#include <ost/gfx/scene.hh>

namespace ost { namespace gfx {

OffscreenBuffer::OffscreenBuffer(unsigned int width, unsigned int height, 
                                 const OffscreenBufferFormat& f, bool shared):
  width_(width), height_(height), valid_(false), active_(false)
{
  CGLPixelFormatAttribute attributes[]={
    kCGLPFAPBuffer,
    kCGLPFAColorSize, CGLPixelFormatAttribute(f.cbits),
    kCGLPFAAlphaSize, CGLPixelFormatAttribute(f.abits),
    kCGLPFADepthSize, CGLPixelFormatAttribute(f.dbits),
    CGLPixelFormatAttribute(0)
  };
  GLint npix=0;
  CGLError err=CGLChoosePixelFormat(attributes, &pix_format_, &npix);
  if(err) {
    LOG_ERROR("error creating offscreen rendering context. "
               "CGLChoosePixFormat failed:" << CGLErrorString(err));
    return;
  }
  // if a context exists, share resources such as shader programs and display
  // lists.
  if (shared) {
    err=CGLCreateContext(pix_format_, CGLGetCurrentContext(), &context_);    
  } else {
    err=CGLCreateContext(pix_format_, NULL, &context_);
  }
  if(err) {
    LOG_ERROR("error creating offscreen rendering context. "
               "CGLCreateContext failed" << CGLErrorString(err));
    return;
  }
  err=CGLCreatePBuffer(width, height, GL_TEXTURE_RECTANGLE_EXT, GL_RGBA, 0, 
                       &pbuffer_);
  if (err) {
    LOG_ERROR("error creating offscreen rendering context. "
               "CGLCreatePBuffer failed: " << CGLErrorString(err));
    return;
  }
  GLint screen=0;
  assert(CGLGetVirtualScreen(context_, &screen)==0);  
  err=CGLSetPBuffer(context_, pbuffer_, 0, 0, screen);
  if (err) {
    LOG_ERROR("error creating offscreen rendering context. "
               "CGLSetPBuffer failed: " << CGLErrorString(err));
    return;
  }

  valid_=true;
}

bool OffscreenBuffer::Resize(unsigned int width, unsigned int height)
{
  CGLPBufferObj      new_pbuffer;  
  CGLError err=CGLCreatePBuffer(width, height, GL_TEXTURE_RECTANGLE_EXT,
                                GL_RGBA, 0,  &new_pbuffer);
  if (err) {
    LOG_ERROR("error resizing offscreen rendering context: "
               "CGLCreatePBuffer failed: " << CGLErrorString(err));
    return false;
  }
  GLint screen=0;

  // TODO: this assert seems dangerous and unnecessary
  assert(CGLGetVirtualScreen(context_, &screen)==0);  
  //

  err=CGLSetPBuffer(context_, new_pbuffer, 0, 0, screen);
  if (err) {
    LOG_ERROR("error resizing offscreen rendering context. "
               "CGLSetPBuffer failed: " << CGLErrorString(err));
    return false;
  }
  CGLDestroyPBuffer(pbuffer_);
  pbuffer_=new_pbuffer;

  width_=width;
  height_=height;
  return true;
}

bool OffscreenBuffer::MakeActive()
{
  if(active_) return true;

  if (CGLError err=CGLSetCurrentContext(context_)) {
    LOG_ERROR("error switching to offscreen rendering context. "
	       "CGLSetCurrentContext failed: " << CGLErrorString(err));
    return false;
  }
  return true;
}

}} // ns
