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
#include <boost/python.hpp>
using namespace boost::python;

#include <ost/img/map.hh>

namespace {

using namespace ost::img;
using namespace ost;

inline MapHandle c0() {return CreateMap();}

inline MapHandle c1a(const Extent& e)
{
  return CreateMap(e);
}

inline MapHandle c2a(const Extent& e, DataType t)
{
  return CreateMap(e,t);
}

inline MapHandle c3a(const Extent& e, DataType t, DataDomain d)
{
  return CreateMap(e,t,d);
}

#define M_CREATE_MAP(NAME1,NAME2,NAME3,PARAM1,PARAM2) \
inline MapHandle NAME1 (const PARAM1 & p1, const PARAM2 & p2) \
{ \
  return CreateMap(p1,p2); \
} \
inline MapHandle NAME2 (const PARAM1 & p1, const PARAM2 & p2,DataType t) \
{ \
  return CreateMap(p1,p2,t); \
} \
inline MapHandle NAME3 (const PARAM1 & p1, const PARAM2 & p2,DataType t, DataDomain d) \
{ \
  return CreateMap(p1,p2,t,d); \
}

M_CREATE_MAP(c1b,c2b,c3b,Size,Point);
M_CREATE_MAP(c1c,c2c,c3c,Point,Point);
M_CREATE_MAP(c1d,c2d,c3d,Point,Size);

#undef M_CREATE_IMAGE
}


void export_Map()
{

  using namespace ost::img;
  def("CreateMap",c0);
  def("CreateMap",c1a);
  def("CreateMap",c2a);
  def("CreateMap",c3a);
  def("CreateMap",c1b);
  def("CreateMap",c2b);
  def("CreateMap",c3b);
  def("CreateMap",c1c);
  def("CreateMap",c2c);
  def("CreateMap",c3c);
  def("CreateMap",c1d);
  def("CreateMap",c2d);
  def("CreateMap",c3d);
}
