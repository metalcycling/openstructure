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
  boost.python wrapper for img base library

  Author: Ansgar Philippsen
*/

#include <boost/python.hpp>
using namespace boost::python;

#include <ost/base.hh>
#include <ost/message.hh>

#include <ost/img/algorithm.hh>

#include <ost/img/map.hh>


// externally located global functions
void export_Data();
void export_DataAlgorithm();
void export_Extent();
void export_ImageHandle();
void export_ImageList();
void export_ConstImageHandle();
void export_Peak();
void export_Point();
void export_PointList();
void export_Size();
void export_RasterImage();
void export_Units();
void export_Mask();
void export_Phase();
void export_Progress();
void export_Map();

// actual module definitions
BOOST_PYTHON_MODULE(_ost_img)
{
  export_Data();
  export_DataAlgorithm();
  export_Extent();
  export_ImageHandle();
  export_ImageList();
  export_ConstImageHandle();
  export_Point();
  export_Peak();
  export_PointList();
  export_Phase();
  export_RasterImage();
  export_Size();
  export_Mask();
  export_Progress();
  export_Map();

}
