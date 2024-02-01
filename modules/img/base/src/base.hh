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
  Some basic types - since this file is included by
  pretty much every header, it should be as slim as
  possible.

  Author: Ansgar Philippsen
*/

#ifndef IMG_BASE_H
#define IMG_BASE_H

#include <string>
#include <cmath>
#include <math.h>
#include <complex>

#include "dllexport.hh"
#include <ost/base.hh>

#ifdef _MSC_VER
#pragma warning(disable:4251)
#pragma warning(disable:4275)
#pragma warning(disable:4244)
#pragma warning(disable:4231)
#endif

#ifdef IRIX
using  std::cos;
using  std::sin;
#endif

namespace ost { namespace img {

  // String is not always predefined
#ifndef String
  typedef std::basic_string<char> String;
#else
  using String
#endif


#ifdef Complex
  // undef Complex from xorg X.h to avoid clash with our definition
  #undef Complex
#endif
  typedef std::complex<Real> Complex;
  typedef unsigned short Word;

  enum DataType {REAL=0, COMPLEX, WORD, DEFAULTTYPE}; //! underlying data type
  enum DataDomain {SPATIAL=0, FREQUENCY, HALF_FREQUENCY, DEFAULTDOMAIN}; //! data domain

 #ifdef IGNORE
 # undef IGNORE
 #endif
  enum DivZeroMethod {
    SET_ZERO=0,
    DIV_EPS=1,
    INTERPOLATE
  };

  enum DomainToColorMode {MODULUS, PHASECOLOR};


#ifndef uchar
  typedef unsigned char uchar;
#endif

#ifndef ushort
  typedef unsigned short int ushort;
#endif

#ifndef uint
  typedef unsigned int uint;
#endif

}} // namespace

#endif
