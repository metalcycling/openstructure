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
#ifndef GEOM_EXCEPTION_HH
#define GEOM_EXCEPTION_HH

#include <stdexcept>
#include <string>

#include <ost/geom/module_config.hh>

namespace geom {

/*
  The rationale behind this base class is
  to povide a common exception for all
  geom exceptions, so that implementing
  code may choose to catch exceptions from
  the geom library undifferentiatedly
*/
class DLLEXPORT GeomException: public std::runtime_error {
public:
  GeomException(const String& m):
    std::runtime_error(m) {}
};

class DLLEXPORT DivideByZeroException: public GeomException {
public:
  DivideByZeroException(const String& m="Divide By Zero"):
    GeomException(m) {}
};

class DLLEXPORT OutOfRangeException: public GeomException {
public:
  OutOfRangeException(const String& m="Out Of Range"):
    GeomException(m) {}
};

} // namespace geom

#endif
