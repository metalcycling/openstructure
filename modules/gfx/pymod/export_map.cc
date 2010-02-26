//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2010 by the OpenStructure authors
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

#include <ost/gfx/map_iso.hh>
#include <ost/gfx/map_slab.hh>
using namespace ost;
using namespace ost::gfx;

namespace {

void ms_color_by_01(MapSlab *s, const Gradient& g, float minv, float maxv)
{
  s->ColorBy(g,minv,maxv);
}

void ms_color_by_02(MapSlab *s, const Gradient& g)
{
  s->ColorBy(g);
}

void ms_color_by_03(MapSlab *s, const Color& c1, const Color& c2, float minv, float maxv)
{
  s->ColorBy(c1,c2,minv,maxv);
}

void ms_color_by_04(MapSlab *s, const Color& c1, const Color& c2)
{
  s->ColorBy(c1,c2);
}

} // anon ns

void export_Map()
{
  class_<MapIso, bases<GfxObj>, boost::shared_ptr<MapIso>,
         boost::noncopyable>("MapIso", init<const String&, const ::iplt::MapHandle&, float, optional<uint> >())
    .def("SetLevel",&MapIso::SetLevel)
    .def("GetLevel",&MapIso::GetLevel)
    .def("GetMean", &MapIso::GetMean)
    .def("Rebuild", &MapIso::Rebuild)
    .def("SetNSF",&MapIso::SetNSF)
    .def("SetDebugOctree", &MapIso::SetDebugOctree)    
  ;

  class_<MapSlab, bases<GfxObj>, boost::shared_ptr<MapSlab>,
         boost::noncopyable>("MapSlab",init<const String&,const ::iplt::MapHandle&, const geom::Plane>())
    .def("SetPlane",&MapSlab::SetPlane)
    .def("GetPlane",&MapSlab::GetPlane)
    .def("ColorBy", ms_color_by_01)
    .def("ColorBy", ms_color_by_02)
    .def("ColorBy", ms_color_by_03)
    .def("ColorBy", ms_color_by_04)
  ;
}
