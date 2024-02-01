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
  Author: Ansgar Philippsen
*/

#ifndef IMG_ALG_PEAK_LIST_HH
#define IMG_ALG_PEAK_LIST_HH

#include <vector>

#include <ost/base.hh>
#include "point.hh"

namespace ost { namespace img { 

class DLLEXPORT_OST_IMG_BASE Peak: public Point {
public:
  Peak(int x=0, int y=0, int z=0) : Point(x,y,z), v_() {}
  Peak(const Point& p): Point(p), v_() {}
  Peak(const Point& p, Real v): Point(p), v_(v) {}
  
  Real GetValue() const {return v_;}
private:
  Real v_;
};

typedef std::vector<Peak> PeakList;

//PointList PeakList2PointList(const PeakList&);

DLLEXPORT_OST_IMG_BASE PeakList ImportPeakList(const String& filename);
DLLEXPORT_OST_IMG_BASE void ExportPeakList(const PeakList& l, const String& filename);


}} // ns

#endif
