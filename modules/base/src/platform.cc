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
#include <cstdlib>
#include <boost/filesystem/path.hpp>
#include <ost/platform.hh>
using boost::filesystem::path;
namespace ost {

static String path_prefix="";
void SetPrefixPath(const String& prefix)
{
  path_prefix=prefix;
}

String GetPrefixPath()
{
  //assert(path_prefix!="");
  if(path_prefix=="") {
    char* ost_root=getenv("OST_ROOT");
    if(!ost_root) {
      throw std::runtime_error("missing PrefixPath, and no fallback OST_ROOT env var");
    }
    path_prefix=std::string(ost_root);
  }
  return path_prefix;
}

String GetSharedDataPath()
{
  return (path(GetPrefixPath()) / "share" / "openstructure").string();
}

}
