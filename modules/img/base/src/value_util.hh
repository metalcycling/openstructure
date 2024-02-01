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
  templated value utilities

  Author: Ansgar Philippsen
*/

#ifndef IMG_VALUE_UTIL_H
#define IMG_VALUE_UTIL_H

#include <time.h>
#include <boost/random.hpp>
#include <ost/img/data_types.hh>

#include <ost/img/module_config.hh>

namespace ost { namespace img {

// randomization
// declaration
template <typename T>
DLLEXPORT_OST_IMG_BASE T Random();

namespace {
  boost::mt19937 RandomGenerator(time(NULL));
  boost::uniform_01<boost::mt19937> UniformRandom(RandomGenerator);
}

// specialization
template <>
inline
Complex Random<Complex>() {
  Real r=UniformRandom();
  Real p=UniformRandom()*2.0*M_PI;
  return Complex(r*cos(p),r*sin(p));
}

template <>
inline
Word Random<Word>() 
{
  Real r=UniformRandom();
  return static_cast<Word>(r*65536.0);
}
  
// definition of generic, ie non-complex case
template <typename T>
inline
T Random() {
  return UniformRandom();
}
 
typedef boost::variate_generator<boost::mt19937&, boost::uniform_int<> > UniformIntGenerator;

inline
UniformIntGenerator GetUniformIntGenerator(int min, int max)
{
  boost::uniform_int<> dist(min,max);
  return UniformIntGenerator(RandomGenerator, dist);
}


// value to value conversion
  
// declaration
template <typename V, typename R>
DLLEXPORT_OST_IMG_BASE R Val2Val(const V& v);
  
// specialization
template <>
inline
Real Val2Val<Complex, Real>(const Complex& c) {return std::abs(c);}

template <>
inline
Word Val2Val<Complex, Word>(const Complex& c) {return static_cast<Word>(std::abs(c));}
  
// generic case
template <typename V, typename R>
inline
R Val2Val(const V& v) {return static_cast<R>(v);}


template <typename T>
DLLEXPORT_OST_IMG_BASE DataType Val2Type();

template<>
inline
DataType Val2Type<Complex>() {return COMPLEX;}

template <>
inline
DataType Val2Type<Real>() {return REAL;}

template <>
inline
DataType Val2Type<Word>() {return WORD;}


template <typename T>
DLLEXPORT_OST_IMG_BASE String Val2String();

template<>
inline
String Val2String<Complex>() {return "COMPLEX";}

template <>
inline
String Val2String<Real>() {return "REAL";}

template <>
inline
String Val2String<Word>() {return "WORD";}


}} // namespace

#endif

