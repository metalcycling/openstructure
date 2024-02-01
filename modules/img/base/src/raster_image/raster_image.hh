//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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
  generic raster image

  Author: Ansgar Philippsen, Andreas Schenk
*/

#ifndef IMG_RASTER_IMAGE_H
#define IMG_RASTER_IMAGE_H

#include <ost/base.hh>
#include <ost/img/vecmat.hh>
#include <ost/img/normalizer_impl.hh>

namespace ost { namespace img {

// fw declarations
class Data;
class Point;

//! Generic raster image
/*!
  RasterImage serves as a conversion
  class between img and external libs, 
  such as wxWidgets. The underlying
  representation is always RGB mode,
  saved consecutively in 3-byte triplets.
*/

/*
  this should be recoded to accept the scale (as is) and the
  (fractional) pixel value of the central pixel. This allows
  the pixel range of the target image to be calculated, and
  the user of the ri (e.g. the viewer) can easily translate 
  and zoom without knowledge of the internal offsets, since:

  center_new_x = center_old_x + 2^zoom * mouse_x
  center_new_y = center_old_y + 2^zoom * mouse_y
*/

class DLLEXPORT_OST_IMG_BASE RasterImage {
  typedef unsigned char uchar;
 public:
  //! coloring mode
  enum Mode {GREY=1, PHASECOLOR, SIGNCOLOR};

  //! Single pixel spec
  struct Pixel {
    uchar r,g,b;
    //! default constructor
    Pixel(): r(0), g(0), b(0) {}
    //! copy-constructor
    Pixel(const Pixel& p) {r=p.r; g=p.g; b=p.b;}
    //! initialize with three uchar
    Pixel(uchar rr,uchar gg,uchar bb): r(rr), g(gg), b(bb) {}
    //! initialize with array
    Pixel(uchar rgb[3]) {r=rgb[0]; g=rgb[1]; b=rgb[2];}
  };

  //! Initialization must specify size
  RasterImage(unsigned int width, unsigned int height);

  // Initialization with size and already allocated memory area
  RasterImage(unsigned int width, unsigned int height, uchar *dptr);

  ~RasterImage();

  unsigned int GetWidth() const {return width_;}
  unsigned int GetHeight() const {return height_;}

  //! Set pixel at specified point (no boundary check is performed!)
  void SetPixel(unsigned int x, unsigned int y, const Pixel& p);
  //! Get pixel at specified point (no boundary check is performed!)
  Pixel GetPixel(unsigned int x, unsigned int y) const;

  //! Get raw data pointer (consecutive unsigned char triplets)
  uchar* GetDataPtr() const {return (uchar *)data_;}

  //! Release data pointer
  /*!
    Internally allocated storage will no longer be freed upon
    RasterImage deletion
  */
  uchar* ReleaseDataPtr() {external_data_=true; return (uchar *)data_;}

  //! fill with image data
  void Fill(const Data& d, int logscale, const Vec3& offset, int z, const NormalizerPtr& norm,
            Mode mode=GREY, bool fast_low_mag_flag=true, bool fast_high_mag_flag=true);
  void Fill(const Data& d, int logscale, const Vec3& offset, int z, const NormalizerPtr& norm,
            Mode mode, bool fast_low_mag_flag, bool fast_high_mag_flag,int x1,int y1,int x2,int y2);


 private:
  unsigned int width_,height_;
  Pixel *data_;
  bool external_data_;
 
  int wh2i(int w, int h) const {return h*width_+w;}
};

}} // namespace img

#endif
