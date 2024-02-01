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

#include <ost/gfx/bitmap_io.hh>
using namespace ost;
using namespace ost::gfx;

namespace {
  unsigned int get_width(const Bitmap& bm) {return bm.width;}
  void set_width(Bitmap& bm, unsigned int w) {bm.width=w;}
  unsigned int get_height(const Bitmap& bm) {return bm.height;}
  void set_height(Bitmap& bm, unsigned int w) {bm.height=w;}

  Bitmap import1(const std::string& n) {return ImportBitmap(n);}
  Bitmap import2(const std::string& n, std::string e) {return ImportBitmap(n,e);}
}

void export_bitmap()
{
  class_<Bitmap>("Bitmap",init<>())
    .add_property("width",get_width,set_width)
    .add_property("height",get_height,set_height)
    ;

  def("ExportBitmap",ExportBitmap);
  def("ImportBitmap",import1);
  def("ImportBitmap",import2);
}
