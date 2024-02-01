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
/*
  Author: Valerio Mariani
*/

#ifndef OST_BOOST_FILESYSTEM_HELPER_HH
#define OST_BOOST_FILESYSTEM_HELPER_HH

#include <boost/filesystem/path.hpp>
#include <boost/version.hpp>

namespace {

inline
String BFPathToString(const boost::filesystem::path& path)
{
#if BOOST_FILESYSTEM_VERSION==3 || BOOST_VERSION<103400
  return path.string();
#else
  return path.file_string();
#endif
}

inline String BFPathStem(const boost::filesystem::path& path) {
#if BOOST_FILESYSTEM_VERSION<103400
  String name = BFPathToString(path);
  size_t n = name.rfind('.');
  return name.substr(0, n);
#else
  return path.stem();
#endif
}

}



#endif // OST_BOOST_FILESYSTEM_HELPER
