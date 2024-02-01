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

#ifndef OST_IO_CONVERT_H
#define OST_IO_CONVERT_H

// endian macros deprecated in Boost 1.68
#include <boost/version.hpp>
#if BOOST_VERSION<106800
#include <boost/detail/endian.hpp>
// recent Boost versions use boost/predef/other/endian.h anyway
#ifndef BOOST_ENDIAN_BIG_BYTE
#ifdef BOOST_BIG_ENDIAN
#define BOOST_ENDIAN_BIG_BYTE 1
#else
#define BOOST_ENDIAN_BIG_BYTE 0
#endif
#endif
#else
// see https://www.boost.org/doc/libs/1_68_0/doc/html/predef/reference.html#predef.reference.other_macros.boost_endian
// defines BOOST_ENDIAN_BIG_BYTE and BOOST_ENDIAN_LITTLE_BYTE
// -> set to 1 or 0 respectively of endianness of machine
#include <boost/predef/other/endian.h>
#endif

namespace ost { namespace io {

enum Endianess {
  OST_BIG_ENDIAN,OST_LITTLE_ENDIAN,OST_VAX_DATA,
#if BOOST_ENDIAN_BIG_BYTE
  OST_LOCAL_ENDIAN = OST_BIG_ENDIAN
#else
  OST_LOCAL_ENDIAN = OST_LITTLE_ENDIAN
#endif
};



template<int CONVERSIONTYPE, class VALUETYPE>
class Convert{
public:
  static VALUETYPE From(const VALUETYPE& value);
  static VALUETYPE To(const VALUETYPE& value);
  static void FromIP(VALUETYPE* value);
  static void ToIP(VALUETYPE* value);
};


}  // io ns



} // ost ns


#endif


