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
#ifndef OST_GFX_BITMAP_EXPORT_HH
#define OST_GFX_BITMAP_EXPORT_HH

#include <boost/shared_array.hpp>

/*
  bitmap export for raw 2D image data

  Author: Ansgar Philippsen
*/

#include <ost/gfx/module_config.hh>

namespace ost { namespace gfx {

// very rudimentary bitmap support
// TODO: gl tex mapping association
struct Bitmap
{
  /*
   channels:
     1 = greyscale
     2 = grey + alpha
     3 = rgb
     4 = rgb + alpha
   */
  unsigned int channels;
  unsigned int width, height;
  boost::shared_array<unsigned char> data;
};

void DLLEXPORT_OST_GFX ExportBitmap(const String& fname, std::string ext, unsigned int width, unsigned int height,unsigned char* data);

Bitmap DLLEXPORT_OST_GFX ImportBitmap(const String& fname, std::string ext="");


}} // ns


#endif
