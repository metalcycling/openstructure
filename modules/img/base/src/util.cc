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

#ifndef IMG_UTIL_H
#error "util.cc must be included from util.hh"
#endif

namespace ost { namespace img { namespace util {

template <typename T> 
void swap_buf(T* p, int n)
{
  uchar *a,tmp;
  
  for(int k=0;k<n;k++) {
    a=reinterpret_cast<uchar *>(&p[k]);
    
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

}}} // namespace
