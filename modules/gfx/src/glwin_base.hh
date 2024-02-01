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
#ifndef OST_GFX_GLWIN_BASE_HH
#define OST_GFX_GLWIN_BASE_HH

/*
  abstract base class for use in the GUI

  Author: Ansgar Philippsen
*/

#include <ost/gfx/module_config.hh>
#include <ost/message.hh>

namespace ost { namespace gfx {

class DLLEXPORT_OST_GFX GLWinBase {
public:
  virtual ~GLWinBase() {}

  virtual void MakeActive() {}

  virtual void DoRefresh() = 0;

  virtual void StatusMessage(const String& m) = 0;

  virtual bool HasStereo() const = 0;

  virtual bool HasMultisample() const = 0;

  // The GLWin is supposed to setup the OpenGL resources. To avoid implementing
  // our own offscreen buffers, image export from OpenGL also happens here
  virtual void Export(const String& fname, unsigned int width,
                      unsigned int height, bool transparent) {
  	throw ost::Error("Image export not implemented for this GLWin");
  }

  virtual void Export(const String& fname, unsigned int width,
                      unsigned int height, int max_samples, 
                      bool transparent) {
  	throw ost::Error("Image export with multi-sampling not implemented for "
  		             "this GLWin");
  }

  virtual void Export(const String& fname, bool transparent) {
  	throw ost::Error("Image export not implemented for this GLWin");
  }

};

}}

#endif
