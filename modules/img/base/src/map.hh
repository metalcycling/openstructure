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
#ifndef OST_MAP_HH
#define OST_MAP_HH

#include <ost/img/image.hh>

namespace ost { namespace img {

typedef ImageHandle MapHandle;

ImageHandle DLLEXPORT_OST_IMG_BASE CreateMap(const Extent& ext=Extent(),
                        DataType type=REAL,
                        DataDomain dom=SPATIAL);

//! convenience variant to create images
ImageHandle DLLEXPORT_OST_IMG_BASE CreateMap(const Size& s, const Point& o,
                                  DataType type=REAL,
                                  DataDomain dom=SPATIAL);

ImageHandle DLLEXPORT_OST_IMG_BASE CreateMap(const Point& p1, const Point& p2,
                                  DataType type=REAL,
                                  DataDomain dom=SPATIAL);

ImageHandle DLLEXPORT_OST_IMG_BASE CreateMap(const Point& p1, const Size& s,
                                  DataType type=REAL,
                                  DataDomain dom=SPATIAL);

}}

#endif
