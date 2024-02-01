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
  image factory

  Author: Ansgar Philippsen
*/

#include <ost/message.hh>

#include "image_factory.hh"

#include "extent.hh"
#include "data.hh"
#include "image_handle.hh"
#include "image_state.hh"

namespace ost { namespace img {

ImageHandle DoCreateImage(const Extent& ext,DataType type,DataDomain dom)
{
  return ImageHandle::Create(ext,type,dom);
}


ImageHandle CreateImage(const Extent& ext, DataType type, DataDomain dom)
{
  return DoCreateImage(ext,type,dom);
}

ImageHandle CreateImage(const Size& s, DataType type, DataDomain dom)
{
  return DoCreateImage(Extent(s),type,dom);
}

ImageHandle CreateImage(const Size& s, const Point& o,DataType type,DataDomain dom)
{
  return CreateImage(Extent(s,o),type,dom);
}

ImageHandle CreateImage(const Point& p1, const Point& p2,DataType type,DataDomain dom)
{
  return CreateImage(Extent(p1,p2),type,dom);
}

ImageHandle CreateImage(const Point& p1, const Size& s,DataType type,DataDomain dom)
{
  return CreateImage(Extent(p1,s),type,dom);
}

ImageHandle GenerateImage(const Data& d)
{
  return GenerateImage(d,d.GetExtent());
}

ImageHandle GenerateImage(const Data& d, const Extent& e)
{
  ImageHandle h=CreateImage(e,d.GetType(),d.GetDomain());

  if(h.GetType()==REAL) {
    for(ExtentIterator it(h.GetExtent());!it.AtEnd();++it) {
      h.SetReal(it,d.GetReal(it));
    }
  } else {
    for(ExtentIterator it(h.GetExtent());!it.AtEnd();++it) {
      h.SetComplex(it,d.GetComplex(it));
    }
  }
  return h;
}

}};
