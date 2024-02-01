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
#ifndef OST_BASE_DYN_CAST_HH
#define OST_BASE_DYN_CAST_HH

#include <boost/shared_ptr.hpp>

namespace ost {

/// \brief a convenient shortcut for the painfully long 
///   boost::dynamic_pointer_cast
template <typename T, typename F>
inline boost::shared_ptr<T> dyn_cast(F f)
{
  return boost::dynamic_pointer_cast<T>(f);
}

}

#endif
