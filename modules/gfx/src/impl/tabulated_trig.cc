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

#include <ost/gfx/impl/tabulated_trig.hh>

/*
   Author: Marco Biasini
 */
namespace ost { namespace gfx { namespace impl {
    

class CosTable {
public:
  const static int TAB_ENTRIES=360;
  CosTable() 
  {
    float delta=2.0*M_PI/TAB_ENTRIES;
    for (int i=0; i<TAB_ENTRIES; ++i) {
      tab_[i]=cos(delta*static_cast<float>(i));
    }
  }
  float Get(float angle) const 
  {
    if (angle<0.0) {
      angle=-angle;
    }
    float bounded=fmod(angle, float(2.0*M_PI));
    float index=TAB_ENTRIES*bounded/(2.0*M_PI);
    int bef=static_cast<int>(index);
    float t=index-bef;
    return tab_[bef]*(1.0f-t)+tab_[(bef+1) % TAB_ENTRIES]*t;
  }
  float tab_[TAB_ENTRIES];      
};

float TabCos(float angle)
{
  return cos(angle);
  /*
  static CosTable table;
  return table.Get(angle);*/
}

float TabSin(float angle)
{
  return sin(angle);
  /*
  return TabCos(angle-M_PI/2.0);
  */
}


}}}
