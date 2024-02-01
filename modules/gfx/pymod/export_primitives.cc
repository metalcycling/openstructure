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

#include <ost/gfx/primitives.hh>
using namespace ost;
using namespace ost::gfx;

void export_primitives()
{
  class_<Primitive, bases<GfxNode>, boost::noncopyable>("Primitive", no_init)
    .def("HasOutline", &Primitive::HasOutline)
    .def("HasFill", &Primitive::HasFill)    
    .def("SetFill", &Primitive::SetFill)
    .def("SetFillColor", &Primitive::SetFillColor)    
    .def("SetOutlineColor", &Primitive::SetOutlineColor)
    .def("SetOutline", &Primitive::SetOutline)            
    .def("SetLineWidth", &Primitive::SetLineWidth)
    .def("GetLineWidth", &Primitive::GetLineWidth)
    .def("GetFillColor", &Primitive::GetFillColor, 
         return_value_policy<copy_const_reference>())
    .def("GetOutlineColor", &Primitive::GetOutlineColor,
         return_value_policy<copy_const_reference>())
  ;
  class_<Cuboid, bases<Primitive>, 
         boost::noncopyable>("Cuboid", init<const String&, const geom::Cuboid>())
  ; 
  class_<Quad, bases<Primitive>, 
       boost::noncopyable>("Quad", init<const String&, const geom::Vec3&,
                                        const geom::Vec3&, const geom::Vec3&,
                                        const geom::Vec3&>())
  ;  
}
