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
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
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

bool get_gds() {return MapIso::global_downsampling_flag;}
void set_gds(bool f) {MapIso::global_downsampling_flag=f;}

} // anon ns

void export_Map()
{

  enum_<MapIsoType>("MapIsoType")
    .value("ORIGINAL_MAP", ORIGINAL_MAP)
    .value("DOWNSAMPLED_MAP", DOWNSAMPLED_MAP)
  ;

  class_<MapIso, bases<GfxObj>, boost::shared_ptr<MapIso>,
         boost::noncopyable>("MapIso", init<const String&, const ::img::MapHandle&, float, optional<uint> >())
    .def("SetLevel",&MapIso::SetLevel)
    .def("GetLevel",&MapIso::GetLevel)
    .def("GetMinLevel",&MapIso::GetMinLevel)
    .def("GetMaxLevel",&MapIso::GetMaxLevel)
    .def("GetMean", &MapIso::GetMean)
    .def("GetHistogram",&MapIso::GetHistogram)
    .def("SetHistogramBinCount",&MapIso::SetHistogramBinCount)
    .def("GetHistogramBinCount",&MapIso::GetHistogramBinCount)

    .def("GetMap", &MapIso::GetMap,return_value_policy<reference_existing_object>())
    .def("GetOriginalMap", &MapIso::GetOriginalMap,return_value_policy<reference_existing_object>())
    .def("GetDownsampledMap", &MapIso::GetDownsampledMap,return_value_policy<reference_existing_object>())
    .def("ShowDownsampledMap", &MapIso::ShowDownsampledMap)
    .def("ShowOriginalMap", &MapIso::ShowOriginalMap)
    .def("IsDownsampledMapAvailable", &MapIso::IsDownsampledMapAvailable)
    .def("GetShownMapType", &MapIso::GetShownMapType)
    .def("MakeOctreeDirty", &MapIso::MakeOctreeDirty)
    .def("IfOctreeDirty", &MapIso::IfOctreeDirty)
    .def("Rebuild", &MapIso::Rebuild)
    .def("SetNSF",&MapIso::SetNSF)
    .def("SetColor", &MapIso::SetColor)
    .def("GetColor", &MapIso::GetColor, return_value_policy<copy_const_reference>())
    .def("SetDebugOctree", &MapIso::SetDebugOctree)    
    .add_static_property("global_downsampling_flag",get_gds,set_gds)
  ;

  class_<MapSlab, bases<GfxObj>, boost::shared_ptr<MapSlab>,
         boost::noncopyable>("MapSlab",init<const String&,const ::img::MapHandle&, const geom::Plane>())
    .def("SetPlane",&MapSlab::SetPlane)
    .def("GetPlane",&MapSlab::GetPlane)
    .def("ColorBy", ms_color_by_01)
    .def("ColorBy", ms_color_by_02)
    .def("ColorBy", ms_color_by_03)
    .def("ColorBy", ms_color_by_04)
  ;
}
