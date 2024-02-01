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
#ifndef OST_GFX_TEXTURE_HH
#define OST_GFX_TEXTURE_HH

/*
  texture

  Author: Ansgar Philippsen
*/

#include <boost/shared_array.hpp>

#include "module_config.hh"
#include "glext_include.hh"

namespace ost { namespace gfx {

struct Bitmap;

class Texture
{
public:
  Texture():
    w_(0),
    h_(0),
    d_()
  {}

  Texture(GLint w, GLint h):
    w_(w),
    h_(h),
    d_(new float[4*w*h])
  {}

  Texture(const Bitmap& b);

  bool IsValid() const {return static_cast<bool>(d_);}

  float* data() {return &d_[0];}
  
  GLint width() const {return w_;}
  GLint height() const {return h_;}
  GLint format() const {return GL_RGBA;}
  GLint type() const {return GL_FLOAT;}
  GLint border() const {return 0;}

private:
  GLint w_,h_;
  boost::shared_array<float> d_;
};

}} // ns

#endif
