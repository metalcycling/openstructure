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
#ifndef OST_INVALID_HANDLE_HH
#define OST_INVALID_HANDLE_HH

#include <ost/message.hh>

#include <ost/module_config.hh>
#include <ost/message.hh>

namespace ost {

/// \brief Signals access of member functions of invalid handles
class DLLEXPORT InvalidHandle : public Error {
public:
  InvalidHandle()
    : Error("Can not access invalid handle or view") 
  {}
}; 

template <typename H>
inline void CheckHandleValidity(const H& handle)
{
  if (!handle.IsValid())
    throw InvalidHandle();
}

}

#endif

