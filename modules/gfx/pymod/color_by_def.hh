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
#ifndef OST_GFX_PYMOD_COLOR_BY_DEF_HH
#define OST_GFX_PYMOD_COLOR_BY_DEF_HH

#include <ost/config.hh>
#include <ost/gfx/gfx_object.hh>
#include <ost/img/alg/stat.hh>


using namespace ost;
using namespace ost::gfx;

namespace ost_gfx {

inline void color_by_e0(GfxObj* go,
                const mol::EntityView& ev, 
                const String& prop,
                const Gradient& g, float minv, float maxv)
{
  go->ColorBy(ev,prop,g,minv,maxv);
}

inline void color_by_e1(GfxObj* go,
			const mol::EntityView& ev, 
			const String& prop,
			const Color& c1, const Color& c2, float minv, float maxv)
{
  Gradient g;
  g.SetColorAt(0.0,c1);
  g.SetColorAt(1.0,c2);
  go->ColorBy(ev,prop,g,minv,maxv);
}

inline void color_by_e2(GfxObj* go,
			const mol::EntityHandle& eh, 
			const String& prop,
			const Gradient& g, float minv, float maxv)
{
  go->ColorBy(eh.CreateFullView(),prop,g,minv,maxv);
}

inline void color_by_e3(GfxObj* go,
			const mol::EntityHandle& eh, 
			const String& prop,
			const Color& c1, const Color& c2, float minv, float maxv)
{
  Gradient g;
  g.SetColorAt(0.0,c1);
  g.SetColorAt(1.0,c2);
  go->ColorBy(eh.CreateFullView(),prop,g,minv,maxv);
}

inline void color_by_e4(GfxObj* go,
			const mol::EntityView& ev, 
			const String& prop,
			const Gradient& g)
{
  std::pair<float,float> minmax = ev.GetMinMax(prop);
  go->ColorBy(ev,prop,g,minmax.first, minmax.second);
}

inline void color_by_e5(GfxObj* go,
			const mol::EntityView& ev, 
			const String& prop,
			const Color& c1, const Color& c2)
{
  Gradient g;
  g.SetColorAt(0.0,c1);
  g.SetColorAt(1.0,c2);
  std::pair<float,float> minmax = ev.GetMinMax(prop);
  go->ColorBy(ev,prop,g,minmax.first, minmax.second);
}

inline void color_by_e6(GfxObj* go,
			const mol::EntityHandle& eh, 
			const String& prop,
			const Gradient& g)
{
  color_by_e4(go,eh.CreateFullView(),prop,g);
}

inline void color_by_e7(GfxObj* go,
			const mol::EntityHandle& eh, 
			const String& prop,
			const Color& c1, const Color& c2)
{
  color_by_e5(go,eh.CreateFullView(),prop,c1,c2);
}

inline void color_by_m0(GfxObj* go,
			const ::img::MapHandle& mh, 
			const String& prop,
			const Gradient& g, float minv, float maxv)
{
  go->ColorBy(mh,prop,g,minv,maxv);
}

inline void color_by_m1(GfxObj* go,
			  const ::img::MapHandle& mh, 
			  const String& prop,
			  const Gradient& g)
{
  ost::img::alg::Stat stat;
  mh.Apply(stat);
  float min = static_cast<float>(stat.GetMinimum());
  float max = static_cast<float>(stat.GetMaximum());
  std::pair<float,float> minmax = std::make_pair(min,max);
  go->ColorBy(mh,prop,g,minmax.first, minmax.second);
}

inline void color_by_m2(GfxObj* go,
			const ::img::MapHandle& mh, 
			const String& prop,
			const Color& c1, const Color& c2)
{
  Gradient g;
  g.SetColorAt(0.0,c1);
  g.SetColorAt(1.0,c2);
  ost::img::alg::Stat stat;
  mh.Apply(stat);
  float min = static_cast<float>(stat.GetMinimum());
  float max = static_cast<float>(stat.GetMaximum());
  std::pair<float,float> minmax = std::make_pair(min,max);
  go->ColorBy(mh,prop,g,minmax.first, minmax.second);
}

inline void color_by_m3(GfxObj* go,
			const ::img::MapHandle& mh, 
			const String& prop,
			const Color& c1, const Color& c2, float minv, float maxv)
{
  Gradient g;
  g.SetColorAt(0.0,c1);
  g.SetColorAt(1.0,c2);
  go->ColorBy(mh,prop,g,minv, maxv);
}

}

#define COLOR_BY_DEF_STD() \
    .def("ColorBy",ost_gfx::color_by_e0)\
    .def("ColorBy",ost_gfx::color_by_e1)\
    .def("ColorBy",ost_gfx::color_by_e2)\
    .def("ColorBy",ost_gfx::color_by_e3)\
    .def("ColorBy",ost_gfx::color_by_e4)\
    .def("ColorBy",ost_gfx::color_by_e5)\
    .def("ColorBy",ost_gfx::color_by_e6)\
    .def("ColorBy",ost_gfx::color_by_e7)

#  define COLOR_BY_DEF_MAP() \
    .def("ColorBy",ost_gfx::color_by_m0)\
    .def("ColorBy",ost_gfx::color_by_m1)\
    .def("ColorBy",ost_gfx::color_by_m2)\
    .def("ColorBy",ost_gfx::color_by_m3)

#define COLOR_BY_DEF() \
  COLOR_BY_DEF_STD() \
  COLOR_BY_DEF_MAP()
#endif

