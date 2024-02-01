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
#ifndef OST_PLATFORM_HH
#define OST_PLATFORM_HH

#include <ost/module_config.hh>
#include <ost/base.hh>

/*
  Author: Marco Biasini
 */
namespace ost { 

/// \brief set path prefix
void DLLEXPORT_OST_BASE SetPrefixPath(const String& prefix);

/// \brief get the path prefix
String DLLEXPORT_OST_BASE GetPrefixPath();

String DLLEXPORT_OST_BASE GetSharedDataPath();

}

#endif
