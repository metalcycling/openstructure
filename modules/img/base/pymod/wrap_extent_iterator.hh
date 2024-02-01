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
  basic image handle, manages shared instances of an image

  offers the Image interface

  Author: Ansgar Philippsen
*/


#ifndef WRAP_EXTENT_ITERATOR_HH_
#define WRAP_EXTENT_ITERATOR_HH_


#include <ost/img/extent.hh>
#include <ost/img/point.hh>
#include <ost/img/extent_iterator.hh>
#include <ost/img/image.hh>
#include <boost/python.hpp>
using namespace boost::python;

namespace ost { namespace img {

class WrapExtentIterator {
public:
  WrapExtentIterator(const ExtentIterator& it): it_(it) {}

  ExtentIterator& Iter() {return it_;}

  Point Next() {
    if(it_.AtEnd()) {
      PyErr_SetObject(PyExc_StopIteration, Py_None);
      throw_error_already_set();
    }

    Point p=it_;
    it_.Inc();
    return p;
  }

  static WrapExtentIterator Create(ExtentIterator& it) 
  {
    return WrapExtentIterator(it);
  }

  static WrapExtentIterator Create2(const Extent& ext) 
  {
    return WrapExtentIterator(ExtentIterator(ext));
  }
  static WrapExtentIterator Create3(const ImageHandle& ih) 
  {
    return WrapExtentIterator(ih.GetIterator());
  }
private:
  ExtentIterator it_;
};

}} // namespace

#endif /*WRAP_EXTENT_ITERATOR_HH_*/
