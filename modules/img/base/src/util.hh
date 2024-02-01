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
  various little things

  Author: Ansgar Philippsen
*/

#ifndef IMG_UTIL_H
#define IMG_UTIL_H

#include <cmath>

#include <ost/img/module_config.hh>

namespace ost { namespace img { namespace util {

// byte swapping stuff

template <typename T> void swap_buf(T* p, int n);

inline
void swap_short(short int *p, short int n) {swap_buf<short int>(p,n);}
inline
void swap_ushort(ushort *p, short int n) {swap_buf<unsigned short int>(p,n);}
inline
void swap_int(int *p, int n) {swap_buf<int>(p,n);}
inline
void swap_uint(uint *p, int n) {swap_buf<uint>(p,n);}
inline
void swap_float(float *p, int n) {swap_buf<float>(p,n);}
inline
void swap_double(Real *p, int n) {swap_buf<Real>(p,n);}
inline
void swap_long(long *p, int n) {swap_buf<long>(p,n);}

}}} // namespace img

// definitions
#include "util.cc"

#endif


