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
  Authors: Andreas Schenk, Ansgar Philippsen
*/

#include <boost/python.hpp>
using namespace boost::python;

#include <ost/img/mask.hh>
#include <ost/config.hh>
#if(OST_INFO_ENABLED)
#include <ost/img/mask_info_convert.hh>
#endif

namespace ost { namespace img {

MaskPtr (*mask1)(const Polygon2&)      =&Mask;
MaskPtr  (*mask2)(const Extent&)       =&Mask;
MaskPtr  (*mask3)(const Circle2&)       =&Mask;
MaskPtr  (*mask4)(const Sphere&)       =&Mask;


bool isinside_wrapv3(MaskPtr* mptr,const Vec2& v)
{
  return (*mptr)->IsInside(v); 
}
bool isoutside_wrapv3(MaskPtr* mptr,const Vec2& v)
{
  return (*mptr)->IsOutside(v); 
}
bool isinside_wrapv2(MaskPtr* mptr,const Vec2& v)
{
  return (*mptr)->IsInside(v);
}
bool isoutside_wrapv2(MaskPtr* mptr,const Vec2& v)
{
  return (*mptr)->IsOutside(v);
}

void shift_wrapv2(MaskPtr* mptr,const Vec2& v)
{
  (*mptr)->Shift(v);
}
void shift_wrapv3(MaskPtr* mptr,const Vec2& v)
{
  (*mptr)->Shift(v); 
}
void expand_wrap(MaskPtr* mptr,Real d)
{
  (*mptr)->Expand(d); 
}
void scale_wrap(MaskPtr* mptr,Real d)
{
  (*mptr)->Scale(d); 
}

MaskPtr clone_wrap(MaskPtr* mptr)
{
  return (*mptr)->Clone();
}

}}//ns

void export_Mask()
{
  using namespace ost::img;

  class_<MaskPtr>("MaskPtr", no_init)
    .def(self ^ self)
    .def(self & self)
    .def(self | self)
    .def(~self)
    .def("Clone",clone_wrap)
    .def("IsInside",isinside_wrapv3)
    .def("IsOutside",isoutside_wrapv3)
    .def("IsInside",isinside_wrapv2)
    .def("IsOutside",isoutside_wrapv2)
    .def("Shift",shift_wrapv3)
    .def("Shift",shift_wrapv2)
    .def("Expand",expand_wrap)
    .def("Scale",scale_wrap)
    ;
  
  def("Mask",mask1);
  def("Mask",mask2);
  def("Mask",mask3);
  def("Mask",mask4);
  
#if(OST_INFO_ENABLED)
  def("InfoToMask",InfoToMask);
  def("MaskToInfo",MaskToInfo);
#endif
  
}
