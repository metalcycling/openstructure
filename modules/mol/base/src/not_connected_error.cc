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
#include "not_connected_error.hh"

namespace ost { namespace mol {
  
NotConnectedError::NotConnectedError(AtomHandle one, 
                                     AtomHandle two)
  : Error(one.GetQualifiedName()+" and "+
                       two.GetQualifiedName()+" are not connected"), 
    one_(one), two_(two) {
}

std::pair<AtomHandle, AtomHandle> NotConnectedError::GetAtoms() const {
  return std::make_pair(one_, two_);
}

}} // ns
