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
  ImageFormat

  Author: Andreas Schenk
*/

#ifndef IMG_IMAGE_FORMAT_H
#define IMG_IMAGE_FORMAT_H

#include <ost/img/data_types.hh>
#include <ost/base.hh>
#include <ost/io/convert.hh>
#include <ost/io/module_config.hh>

namespace ost { namespace io {

class DLLEXPORT_OST_IO ImageFormatBase
{
 protected:

  ImageFormatBase():
    format_string_("")
  {};

  ImageFormatBase(const String& formatstring):
    format_string_(formatstring)
  {};

 public:

  Real GetMaximum() const { return 1.0;}
  Real GetMinimum() const { return 0.0;}


  template <typename T>
  const T& As() const
  {
    return reinterpret_cast<const T&>(*this);
  }
  String GetFormatString() const {return format_string_;}

 private:

  String format_string_;

};

enum Format
{
  OST_BIT8_FORMAT,OST_BIT16_FORMAT,OST_BIT32_FORMAT,OST_FLOAT_FORMAT,OST_DOUBLE_FORMAT,OST_DEFAULT_FORMAT
};

enum Subformat {
  MRC_OLD_FORMAT,MRC_NEW_FORMAT,MRC_AUTO_FORMAT
};



class DLLEXPORT_OST_IO UndefinedImageFormat : public ImageFormatBase
{
 public:
  UndefinedImageFormat():
    ImageFormatBase(FORMAT_STRING)
  {}

  static String FORMAT_STRING;

};


}} // namespaces

#endif
