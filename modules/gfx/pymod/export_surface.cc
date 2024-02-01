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
using namespace boost::python;

#include <ost/gfx/surface.hh>

using namespace ost;
using namespace ost::gfx;

#include "color_by_def.hh"

namespace {

void color_by_01(Surface* s,
                 const String& prop, 
                 const Gradient& gradient,
                 float minv,float maxv,
                 mol::Prop::Level hint)
{
  s->ColorBy(prop,gradient,minv,maxv,hint);
}

void color_by_02(Surface* s,
                 const String& prop, 
                 const Gradient& gradient,
                 float minv,float maxv)
{
  s->ColorBy(prop,gradient,minv,maxv);
}

void color_by_03(Surface* s,
                 const String& prop, 
                 const Gradient& gradient,
                 mol::Prop::Level hint)
{
  s->ColorBy(prop,gradient,hint);
}

void color_by_04(Surface* s,
                 const String& prop, 
                 const Gradient& gradient)
{
  s->ColorBy(prop,gradient);
}

void color_by_05(Surface* s,
                 const String& prop, 
                 const Color& c1, const Color& c2,
                 float minv,float maxv,
                 mol::Prop::Level hint)
{
  s->ColorBy(prop,c1,c2,minv,maxv,hint);
}

void color_by_06(Surface* s,
                 const String& prop, 
                 const Color& c1, const Color& c2,
                 float minv,float maxv)
{
  s->ColorBy(prop,c1,c2,minv,maxv);
}

void color_by_07(Surface* s,
                 const String& prop, 
                 const Color& c1, const Color& c2,
                 mol::Prop::Level hint)
{
  s->ColorBy(prop,c1,c2,hint);
}

void color_by_08(Surface* s,
                 const String& prop, 
                 const Color& c1, const Color& c2)
{
  s->ColorBy(prop,c1,c2);
}

void set_color1(Surface* s, const Color& c) 
{
  s->SetColor(c);
}

void set_color2(Surface* s, const Color& c, const String& t) 
{
  s->SetColor(c,t);
}

} // anon ns

void export_Surface()
{
  class_<Surface, boost::shared_ptr<Surface>, bases<GfxObj>, boost::noncopyable>("Surface", init<const String&, const mol::SurfaceHandle&>())
    .def("Replace",&Surface::Replace)
    .def("ColorBy", color_by_01)
    .def("ColorBy", color_by_02)
    .def("ColorBy", color_by_03)
    .def("ColorBy", color_by_04)
    .def("ColorBy", color_by_05)
    .def("ColorBy", color_by_06)
    .def("ColorBy", color_by_07)
    .def("ColorBy", color_by_08)
    COLOR_BY_DEF()
    .def("SetColor", set_color1)
    .def("SetColor", set_color2)
    .def("CleanColorOps", &Surface::CleanColorOps)
    .def("ReapplyColorOps", &Surface::ReapplyColorOps)
    .def("GetHandle", &Surface::GetHandle)
    .add_property("handle", &Surface::GetHandle)
    ;

}
