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
#ifndef OST_MOL_SECONDARY_STRUCTURE_HH
#define OST_MOL_SECONDARY_STRUCTURE_HH

#include <ost/mol/module_config.hh>

namespace ost { namespace mol {

/// \brief Secondary structure types as defined by DSSP. 
/// For convenience, the enum values match the characters used in the DSP file 
/// format.
struct DLLEXPORT_OST_MOL SecStructure {
  typedef enum {
    ALPHA_HELIX     ='H',
    PI_HELIX        ='I',
    THREE_TEN_HELIX ='G',
    TURN            ='T',
    EXTENDED        ='E',
    BETA_BRIDGE     ='B',
    BEND            ='S',
    COIL            ='C'
  } Type;
  explicit SecStructure(char type) : type_(Type(type)) {}
  explicit SecStructure(Type type) : type_(type) {}  
  SecStructure() : type_(COIL) { }
  bool operator==(const Type& rhs) const 
  { 
    return type_==rhs; 
  }
  
  bool operator==(const SecStructure& rhs) const 
  { 
    return type_==rhs.type_; 
  }
  
  bool operator!=(const Type& rhs) const
  { 
    return !this->operator==(rhs); 
  }
  bool operator!=(const SecStructure& rhs) const 
  { 
    return !this->operator==(rhs); 
  }
  
  bool IsHelical() const 
  { 
    return type_==ALPHA_HELIX || type_==PI_HELIX || type_==THREE_TEN_HELIX;
  }
  
  bool IsExtended() const
  {
    return type_==EXTENDED || type_==BETA_BRIDGE;
  }
  
  bool IsCoil() const 
  {
    return type_==COIL || type_==BEND || type_==TURN;
  }
  operator char() const { return type_; }
private:
  Type type_;
};

}} // ns
#endif
