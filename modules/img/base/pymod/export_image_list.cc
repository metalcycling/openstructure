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
  Author: Andreas Schenk
*/

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;

#include <ost/img/algorithm.hh>
#include <ost/img/image.hh>
#include <ost/img/image_list.hh>

using namespace ost::img;
using namespace ost;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(gallery_overloads, GetGallery, 0, 2)



void export_ImageList()
{
  void (ImageList::*apply_const)(NonModAlgorithm&) const = &ImageList::Apply;
  void (ImageList::*apply_const_ip)(NonModAlgorithm&) const = &ImageList::ApplyIP;
  void (ImageList::*apply_modip_ip)(ModIPAlgorithm&) = &ImageList::ApplyIP;
  void (ImageList::*apply_cmodip_ip)(const ConstModIPAlgorithm&) = &ImageList::ApplyIP;
  ImageList (ImageList::*apply_modip)(ModIPAlgorithm&) const = &ImageList::Apply;
  ImageList (ImageList::*apply_cmodip)(const ConstModIPAlgorithm&) const = &ImageList::Apply;
  void (ImageList::*apply_modop_ip)(ModOPAlgorithm&) = &ImageList::ApplyIP;
  void (ImageList::*apply_cmodop_ip)(const ConstModOPAlgorithm&) = &ImageList::ApplyIP;
  ImageList (ImageList::*apply_modop)(ModOPAlgorithm&) const = &ImageList::Apply;
  ImageList (ImageList::*apply_cmodop)(const ConstModOPAlgorithm&) const = &ImageList::Apply;

/*
  Important:
  Overloaded binary operators using Real need to be wrapped after 
  corresponding operators using a complex, otherwise python uses the 
  wrong version of the operator. (under Linux)
  Ex.: image-=2.0 is interpreted as image-=complex(2.0) which results in image+=2.0
*/

  class_<ImageList >("ImageList", init<>() )
    .def(vector_indexing_suite<ImageList >())
    .def("GetGallery",&ImageList::GetGallery,gallery_overloads())
    .def("GetImageStack",&ImageList::GetImageStack)
    .def("Apply",apply_const)
    .def("ApplyIP",apply_const_ip)
    .def("ApplyIP",apply_modip_ip)
    .def("ApplyIP",apply_cmodip_ip)
    .def("Apply",apply_modip)
    .def("Apply",apply_cmodip)
    .def("ApplyIP",apply_modop_ip)
    .def("ApplyIP",apply_cmodop_ip)
    .def("Apply",apply_modop)
    .def("Apply",apply_cmodop)

    .def(self += ImageHandle())
    .def(self -= ImageHandle())
    .def(self *= ImageHandle())
    .def(self /= ImageHandle())
    .def(self += Complex())
    .def(self -= Complex())
    .def(self *= Complex())
    .def(self /= Complex())
    .def(self += Real())
    .def(self -= Real())
    .def(self *= Real())
    .def(self /= Real())
    ;
}

