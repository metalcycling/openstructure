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
#ifndef OST_NOT_CONNECTED_ERROR_HH
#define OST_NOT_CONNECTED_ERROR_HH

#include <ost/message.hh>
#include <ost/mol/atom_handle.hh>

namespace ost { namespace mol {
  
/// \brief Error to indicate that two atoms that are supposed to be connected
///        are not.
class DLLEXPORT NotConnectedError : public Error {
public:
  NotConnectedError(AtomHandle one, AtomHandle two);
  virtual ~NotConnectedError() throw() {}
  /// \brief get the two unconnected atoms
  std::pair<AtomHandle, AtomHandle> GetAtoms() const;
private:
  AtomHandle one_;
  AtomHandle two_;
};  

}} //ns

#endif
