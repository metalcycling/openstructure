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
#include <ost/img/size.hh>
#include <ost/img/point.hh>
using namespace boost::python;

namespace ost { namespace img {

int size_len(const Size&) {return 3;}
int size_getitem(const Size& s, int i) {return s[i];}
void size_setitem(Size& s, int i, int v) {s[i]=v;}

}}

void export_Size()
{
  using namespace ost::img;
  // size wrapper
  class_<Size>("Size", init<optional<unsigned int, unsigned int, unsigned int> >(args("width","height","depth")) )
    .def(init<const Point&,const Point&>(args("point1","point2")))
    .def("GetHalf",&Size::GetHalf) 
    .def("GetWidth",&Size::GetWidth)
    .def("GetHeight",&Size::GetHeight)
    .def("GetDepth",&Size::GetDepth)
    .def("GetSlab",&Size::GetSlab)
    .def("GetVolume",&Size::GetVolume)
    .def("GetDim",&Size::GetDim)
    .def("__len__",size_len)
    .def("__getitem__",size_getitem)
    .def("__setitem__",size_setitem)
    .def(self += self)
    .def(self + self)
    .def(self -= self)
    .def(self - self)
    .def(self_ns::str(self))
    ;
}
