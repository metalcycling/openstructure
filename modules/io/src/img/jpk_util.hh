//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
//
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
#ifndef OST_IO_JPK_UTIL_HH
#define OST_IO_JPK_UTIL_HH

namespace ost { namespace io { namespace jpk_tags {

const static unsigned int SLOT_COUNT=0x8080;
const static unsigned int DEFAULT_SLOT=0x8081;
const static unsigned int GRID_REFLECT=0x8045;
const static unsigned int GRID_U_LENGTH=0x8042;
const static unsigned int GRID_V_LENGTH=0x8043;


inline unsigned int ScalingVar1(unsigned int n) 
{ 
  return (0x80A4+(n)*0x30);
}

inline unsigned int ScalingVar2(unsigned int n)     
{ 
  return (0x80A5+(n)* 0x30); 
}

inline unsigned int ScalingMultiply(unsigned int n) 
{ 
  return ScalingVar1(n); 
}

inline unsigned int SlotName(unsigned n)
{ 
  return (0x8090+(n)*0x30); 
}

inline unsigned int ScalingOffset(unsigned int n)   
{ 
  return ScalingVar2(n); 
}
inline unsigned int EncoderUnit(unsigned int n)
{ 
  return (0x80A2+(n)*0x30); 
}



}}}

#endif
