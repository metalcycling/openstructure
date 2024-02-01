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
using namespace boost::python;

#include <ost/message.hh>
#include <ost/img/algorithm.hh>
#include <ost/img/image.hh>
#include "wrap_extent_iterator.hh"

using namespace ost::img;
using namespace ost;

namespace {

inline ImageHandle c0() {return CreateImage();}

inline ImageHandle c1a(const Extent& e) 
{
  return CreateImage(e);
}

inline ImageHandle c2a(const Extent& e, DataType t) 
{
  return CreateImage(e,t);
}

inline ImageHandle c3a(const Extent& e, DataType t, DataDomain d) 
{
  return CreateImage(e,t,d);
}

#define M_CREATE_IMAGE(NAME1,NAME2,NAME3,PARAM1,PARAM2) \
inline ImageHandle NAME1 (const PARAM1 & p1, const PARAM2 & p2) \
{ \
  return CreateImage(p1,p2); \
} \
inline ImageHandle NAME2 (const PARAM1 & p1, const PARAM2 & p2,DataType t) \
{ \
  return CreateImage(p1,p2,t); \
} \
inline ImageHandle NAME3 (const PARAM1 & p1, const PARAM2 & p2,DataType t, DataDomain d) \
{ \
  return CreateImage(p1,p2,t,d); \
}

inline ImageHandle create_im_int(int width, int height, int depth)
{
  return CreateImage(img::Size(width, height, depth));
}

M_CREATE_IMAGE(c1b,c2b,c3b,Size,Point);
M_CREATE_IMAGE(c1c,c2c,c3c,Point,Point);
M_CREATE_IMAGE(c1d,c2d,c3d,Point,Size);

#undef M_CREATE_IMAGE

inline ImageHandle g1a(const Data& d) 
{
  return GenerateImage(d);
}

inline ImageHandle g2a(const Data& d,const Extent& e) 
{
  return GenerateImage(d,e);
}

inline ImageHandle ih_copy1(ImageHandle& h) {return h.Copy();}
inline ImageHandle ih_copy2(ImageHandle& h, bool cc) {return h.Copy(cc);}

img::Point tuple_to_point(size_t dim, tuple t)
{
  size_t tuple_len=len(t);
  if (tuple_len!=dim) {
    throw ost::Error("length of tuple and image dimension must match");
  }
  img::Point p;
  for (size_t i=0; i<tuple_len; ++i) {
    p[i]=extract<int>(t[i]);
  }
  return p;
}

object ih_get_item(ConstImageHandle ih, img::Point p)
{
  if (ih.GetType()==img::REAL) {
    return object(ih.GetReal(p));
  } else {
    return object(ih.GetComplex(p));
  }  
}
object ih_get_item_a(ConstImageHandle ih, tuple t)
{
  img::Point p=tuple_to_point(ih.GetSize().GetDim(), t);
  return ih_get_item(ih, p);
}

object ih_get_item_c(ConstImageHandle ih, int x)
{
  if (ih.GetSize().GetDim()!=1) {
    throw ost::Error("Can't address point of multi-dimensional image with scalar");
  }
  return ih_get_item(ih, img::Point(x, 0, 0));
}
object ih_get_item_b(ConstImageHandle ih, img::Point p)
{
  return ih_get_item(ih, p);
}

void ih_set_item(ImageHandle ih, img::Point p, object value)
{
  if (ih.GetType()==img::REAL) {
    ih.SetReal(p, extract<Real>(value));
  } else {
    ih.SetComplex(p, extract<Complex>(value));
  }
}

void ih_set_item_a(ImageHandle ih, tuple t, object value)
{
  img::Point p=tuple_to_point(ih.GetSize().GetDim(), t);
  ih_set_item(ih, p, value);
}

void ih_set_item_b(ImageHandle ih, img::Point p, object value)
{
  ih_set_item(ih, p, value);
}

void ih_set_item_c(ImageHandle ih, int x, object value)
{
  if (ih.GetSize().GetDim()!=1) {
    throw ost::Error("Can't address point of multi-dimensional image with scalar");
  }  
  ih_set_item(ih, img::Point(x, 0, 0), value);
}


} // anon ns

void export_ImageHandle()
{
  void (ImageHandle::*apply_const)(NonModAlgorithm&) const = &ImageHandle::Apply;
  void (ImageHandle::*apply_const_ip)(NonModAlgorithm&) const = &ImageHandle::ApplyIP;
  void (ImageHandle::*apply_modip_ip)(ModIPAlgorithm&) = &ImageHandle::ApplyIP;
  void (ImageHandle::*apply_cmodip_ip)(const ConstModIPAlgorithm&) = &ImageHandle::ApplyIP;
  ImageHandle (ImageHandle::*apply_modip)(ModIPAlgorithm&) const = &ImageHandle::Apply;
  ImageHandle (ImageHandle::*apply_cmodip)(const ConstModIPAlgorithm&) const = &ImageHandle::Apply;
  void (ImageHandle::*apply_modop_ip)(ModOPAlgorithm&) = &ImageHandle::ApplyIP;
  void (ImageHandle::*apply_cmodop_ip)(const ConstModOPAlgorithm&) = &ImageHandle::ApplyIP;
  ImageHandle (ImageHandle::*apply_modop)(ModOPAlgorithm&) const = &ImageHandle::Apply;
  ImageHandle (ImageHandle::*apply_cmodop)(const ConstModOPAlgorithm&) const = &ImageHandle::Apply;

/*
  Important:
  Overloaded binary operators using Real need to be wrapped after 
  corresponding operators using a complex, otherwise python uses the 
  wrong version of the operator. (under Linux)
  Ex.: image-=2.0 is interpreted as image-=complex(2.0) which results in image+=2.0
*/

  class_<ImageHandle, bases<Data> >("ImageHandle", no_init )
    .def("Copy",ih_copy1) 
    .def("Copy",ih_copy2,args("content_flag")) 
    .def("Extract",&ImageHandle::Extract,args("extent"))
    .def("Paste",&ImageHandle::Paste,args("data"))
    .def("Set",&ImageHandle::Set,args("image"))
    .def("Swap",&ImageHandle::Swap, args("other"))
    .def("Reset",&ImageHandle::Reset, args("ext","type","dom","enc"))
    .def("IsValid",&ImageHandle::IsValid)
    .def("MemSize",&ImageHandle::MemSize)
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
    .def("SetReal",&ImageHandle::SetReal, args("point","value"))
    .def("SetComplex",&ImageHandle::SetComplex, args("point","value"))
    .def("SetSpatialOrigin",&ImageHandle::SetSpatialOrigin)
    .def("GetSpatialOrigin",&ImageHandle::GetSpatialOrigin)
    .add_property("spatial_origin",&ImageHandle::GetSpatialOrigin,&ImageHandle::SetSpatialOrigin)
    .def("CenterSpatialOrigin",&ImageHandle::CenterSpatialOrigin)
    .def("SetAbsoluteOrigin",&ImageHandle::SetAbsoluteOrigin)
    .def("GetAbsoluteOrigin",&ImageHandle::GetAbsoluteOrigin)
    .add_property("absolute_origin",&ImageHandle::GetAbsoluteOrigin,&ImageHandle::SetAbsoluteOrigin)
    .def("IndexToCoord",&ImageHandle::IndexToCoord)
    .def("CoordToIndex",&ImageHandle::CoordToIndex)
    .def("FractionalIndexToCoord",&ImageHandle::FractionalIndexToCoord)
    .def("__iter__",&WrapExtentIterator::Create3)
    .def(self == self)
    .def(self != self)
    .def("__getitem__", ih_get_item_a)
    .def("__getitem__", ih_get_item_b)
    .def("__getitem__", ih_get_item_c)
    .def("__setitem__", ih_set_item_a)
    .def("__setitem__", ih_set_item_b)
    .def("__setitem__", ih_set_item_c)
    .def(self += self)
    .def(self + self)
    .def(self -= self)
    .def(self - self)
    .def(self *= self)
    .def(self * self)
    .def(self /= self)
    .def(self / self)
    .def(self += Complex())
    .def(self + Complex())
    .def(Complex() + self)
    .def(self -= Complex())
    .def(self - Complex())
    .def(Complex() - self)
    .def(self *= Complex())
    .def(self * Complex())
    .def(Complex() * self)
    .def(self /= Complex())
    .def(self / Complex())
    .def(self += Real())
    .def(self + Real())
    .def(Real() + self)
    .def(self -= Real())
    .def(self - Real())
    .def(Real() - self)
    .def(self *= Real())
    .def(self * Real())
    .def(Real() * self)
    .def(self /= Real())
    .def(self / Real())
    ;

  def("CreateImage",c0);
  def("CreateImage", create_im_int, (arg("width"), 
      arg("height")=0, arg("depth")=0));
  def("CreateImage",c1a);
  def("CreateImage",c2a);
  def("CreateImage",c3a);
  def("CreateImage",c1b);
  def("CreateImage",c2b);
  def("CreateImage",c3b);
  def("CreateImage",c1c);
  def("CreateImage",c2c);
  def("CreateImage",c3c);
  def("CreateImage",c1d);
  def("CreateImage",c2d);
  def("CreateImage",c3d);
  def("GenerateImage",g1a);
  def("GenerateImage",g2a);
}

void export_ConstImageHandle()
{
  void (ConstImageHandle::*apply_const)(NonModAlgorithm&) const = &ConstImageHandle::Apply;
  ImageHandle (ConstImageHandle::*apply_modip)(ModIPAlgorithm&) const = &ConstImageHandle::Apply;
  ImageHandle (ConstImageHandle::*apply_modop)(ModOPAlgorithm&) const = &ConstImageHandle::Apply;
  ImageHandle (ConstImageHandle::*apply_cmodip)(const ConstModIPAlgorithm&) const = &ConstImageHandle::Apply;
  ImageHandle (ConstImageHandle::*apply_cmodop)(const ConstModOPAlgorithm&) const = &ConstImageHandle::Apply;

  class_<ConstImageHandle, bases <ConstData> >("ConstImageHandle", init<const ImageHandle&>() )
    .def("IsValid",&ConstImageHandle::IsValid)
    .def("MemSize",&ConstImageHandle::MemSize)
    .def("Apply",apply_const)
    .def("Apply",apply_modip)
    .def("Apply",apply_modop)
    .def("Apply",apply_cmodip)
    .def("Apply",apply_cmodop)
    .def("GetSpatialOrigin",&ConstImageHandle::GetSpatialOrigin)
    .add_property("spatial_origin",&ImageHandle::GetSpatialOrigin)
    .def("GetAbsoluteOrigin",&ConstImageHandle::GetAbsoluteOrigin)
    .add_property("absolute_origin",&ImageHandle::GetAbsoluteOrigin)
    .def("IndexToCoord",&ConstImageHandle::IndexToCoord)
    .def("CoordToIndex",&ConstImageHandle::CoordToIndex)
    .def("FractionalIndexToCoord",&ConstImageHandle::FractionalIndexToCoord)
    .def("__iter__",&WrapExtentIterator::Create3)
    .def(self == self)
    .def(self != self)
    .def(self + Real())
    .def(Real() + self)
    .def(self - Real())
    .def(Real() - self)
    .def(self * Real())
    .def(Real() * self)
    .def(self / Real())
    .def("__getitem__", ih_get_item_a)
    .def("__getitem__", ih_get_item_b)
    .def("__getitem__", ih_get_item_c)
    .def(self + Complex())
    .def(Complex() + self)
    .def(self - Complex())
    .def(Complex() - self)
    .def(self * Complex())
    .def(Complex() * self)
    .def(self / Complex())
     ;

  implicitly_convertible<ImageHandle,ConstImageHandle>();
}
