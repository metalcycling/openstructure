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

#include "point_list.hh"

namespace ost { namespace img {

PointList::PointList():
 std::vector<Point>()
{
}
PointList::PointList(PointList::size_type n):
 std::vector<Point>(n)
{
}
PointList::PointList(PointList::size_type n, const Point& p):
 std::vector<Point>(n,p)
{
}
PointList::PointList( PointList::const_iterator first, PointList::const_iterator last):
 std::vector<Point>(first,last)
{
}

PointList::PointList(const PeakList& plist)
{
  for(PeakList::const_iterator it=plist.begin();it!=plist.end();++it) {
    push_back(Point(*it));
  }
}


}} // namespace img
