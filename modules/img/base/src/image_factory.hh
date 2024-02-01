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
  factory methods to create images, either from
  scratch or based on some other data instance

  Author: Ansgar Philippsen
*/

#ifndef IMG_IMAGE_FACTORY_HH
#define IMG_IMAGE_FACTORY_HH

#include <ost/img/module_config.hh>
#include "extent.hh"
#include "point.hh"
#include "data_types.hh"
#include "size.hh"

namespace ost { namespace img {

class Data;
class ImageHandle;

//! Create image from scratch
/*!
  A handle to a newly created handle is returned. The
  image is defined by an Extent, the type (\c REAL (default) 
  or \c COMPLEX), and the domain (\c SPATIAL (default), 
  \c FREQUENCY or \c HALF_FREQUENCY ).
*/
ImageHandle DLLEXPORT_OST_IMG_BASE CreateImage(const Extent& ext=Extent(),
                                                DataType type=REAL,
                                                DataDomain dom=SPATIAL);

//! convenience variants to create images
ImageHandle DLLEXPORT_OST_IMG_BASE CreateImage(const Size& s,
                                                DataType type=REAL,
                                                DataDomain dom=SPATIAL);

ImageHandle DLLEXPORT_OST_IMG_BASE CreateImage(const Size& s, const Point& o,
                                                DataType type=REAL,
                                                DataDomain dom=SPATIAL);

ImageHandle DLLEXPORT_OST_IMG_BASE CreateImage(const Point& p1, const Point& p2,
                                                DataType type=REAL,
                                                DataDomain dom=SPATIAL);

ImageHandle DLLEXPORT_OST_IMG_BASE CreateImage(const Point& p1, const Size& s,
                                                DataType type=REAL,
                                                DataDomain dom=SPATIAL);

//! Create image based on some data object (ie other image or function)
/*
  The returned image will have the extent of the given Data
*/
ImageHandle DLLEXPORT_OST_IMG_BASE GenerateImage(const Data& d);

//! Create image based on some data object (ie other image or function)
/*!
  The returned image will have the extent given by second parameter
*/
ImageHandle DLLEXPORT_OST_IMG_BASE GenerateImage(const Data& d, const Extent& e);

}} // namespace

#endif
