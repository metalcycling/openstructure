//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2003-2010 by the IPLT authors
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
#ifndef OST_IO_SWAP_UTIL_HH
#define OST_IO_SWAP_UTIL_HH


/*
  byte swapping utils, originally copied from IPLT

  Author: Ansgar Philippsen
*/

#include <cmath>

namespace ost { namespace io {

// byte swapping stuff

template <typename T> void swap_buf(T* p, int n);

inline
void swap_short(short int *p, short int n) {swap_buf<short int>(p,n);}
inline
void swap_ushort(unsigned short int *p, short int n) {swap_buf<unsigned short int>(p,n);}
inline
void swap_int(int *p, int n) {swap_buf<int>(p,n);}
inline
void swap_uint(uint *p, int n) {swap_buf<uint>(p,n);}
inline
void swap_float(float *p, int n) {swap_buf<float>(p,n);}
inline
void swap_double(double *p, int n) {swap_buf<double>(p,n);}
inline
void swap_long(long *p, int n) {swap_buf<long>(p,n);}

//definitions
template <typename T>
void swap_buf(T* p, int n)
{
  unsigned char *a,tmp;

  for(int k=0;k<n;k++) {
    a=reinterpret_cast<unsigned char *>(&p[k]);

    if(sizeof(T)==2) {
      tmp=a[0];
      a[0]=a[1];
      a[1]=tmp;
    } else if(sizeof(T)==4) {
      tmp=a[0];
      a[0]=a[3];
      a[3]=tmp;
      tmp=a[1];
      a[1]=a[2];
      a[2]=tmp;
    } else if(sizeof(T)==8) {
      tmp = a[0];
      a[0] = a[7];
      a[7] = tmp;
      tmp = a[1];
      a[1] = a[6];
      a[6] = tmp;
      tmp = a[2];
      a[2] = a[5];
      a[5] =tmp;
      tmp = a[3];
      a[3] = a[4];
      a[4] = tmp;
    }
  }
}

}} // namespace

#endif


