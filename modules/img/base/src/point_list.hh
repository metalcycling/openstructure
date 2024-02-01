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
  Author: Andreas Schenk
*/

#ifndef POINT_LIST_HH_
#define POINT_LIST_HH_

#include "point.hh"
#include "peak.hh"

namespace ost { namespace img {


//! list of points
/*!
  convenience definition
*/

class DLLEXPORT_OST_IMG_BASE PointList: public std::vector<Point>
{
public:
  PointList();
  PointList(PointList::size_type n);
  PointList(PointList::size_type n, const Point& p);
  PointList( PointList::const_iterator first, PointList::const_iterator last);
  explicit PointList(const PeakList& pl);
};


}} // namespace img

#endif /*POINT_LIST_HH_*/
