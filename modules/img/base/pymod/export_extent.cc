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

#include <boost/python.hpp>
#include "wrap_extent_iterator.hh"

using namespace boost::python;


void export_Extent()
{
  using namespace ost::img;

  bool    (Extent::*contains1)(const Point&) const             = &Extent::Contains;
  bool    (Extent::*contains2)(const Extent&) const             = &Extent::Contains;

  // extent wrapper
  class_<Extent>("Extent", init<>() )
    .def(init<const Extent&>(args("extent")) )
    .def(init<const Point&, const Point&>(args("start","end")) )
    .def(init<const Size&>(args("size")) )
    .def(init<const Size&, const Point&>(args("size","origin")) )
    .def(init<const Point&, const Size&>(args("start","size")) )
    .def("GetStart",&Extent::GetStart,
	 return_value_policy<return_by_value>() )
    .def("GetEnd",&Extent::GetEnd,
	 return_value_policy<return_by_value>() )
    .def("SetStart",&Extent::SetStart, args("start"))   
    .def("SetEnd",&Extent::SetEnd, args("end"))
    .def("GetWidth",&Extent::GetWidth)
    .def("GetHeight",&Extent::GetHeight)
    .def("GetDepth",&Extent::GetDepth)
    .def("GetVolume",&Extent::GetVolume)
    .def("Contains",contains1, args("point"))
    .def("Contains",contains2, args("point"))
    .def("GetCenter",&Extent::GetCenter)
    .def("GetSize",&Extent::GetSize,
	 return_value_policy<return_by_value>() )
    .def("GetDim",&Extent::GetDim)
    .def("WrapAround",&Extent::WrapAround, args("point"))
    .def("Mirror",&Extent::Mirror)
    .def("Shift",&Extent::Shift)
    .def("__iter__",&WrapExtentIterator::Create2)
    .def(self_ns::str(self))
    .def(self == self)
    .def(self != self)
    ;

  // extent iterator wrapper
  class_<WrapExtentIterator>("WrapExtentIterator",no_init)
    .def("__iter__",&WrapExtentIterator::Iter,
	 return_internal_reference<>() )
    .def("next",&WrapExtentIterator::Next)
    .def("__next__",&WrapExtentIterator::Next)
  ;

  class_<ExtentIterator>("ExtentIterator", init<const Extent&,optional<unsigned int> >() )
    .def("__iter__",&WrapExtentIterator::Create)
    ;

  def("Overlap",Overlap);
  def("HasOverlap",HasOverlap);

  implicitly_convertible<Size,Extent>();
}
