//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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
  info path

  Author: Ansgar Philippsen
*/

#include <vector>
#include <iostream>

#include <ost/base.hh>

#include <ost/info/module_config.hh>


namespace ost { namespace info {

//! encapsulates path to info groups and items
/*!
  path elements in the String are separated by the typical path
  separator '/', where a path that starts with a slash denotes
  an absolute path from the root node, while a path without leading
  slash denotes a local path.
*/
class DLLEXPORT_OST_INFO InfoPath {
friend InfoPath operator+(const InfoPath& p1, const InfoPath& p2);
public:
  //! initialize with String
  InfoPath(const String& p);
  InfoPath(const char* cp);
  
  /// \brief retrieve list of path elements
  std::vector<String> GetList() const {return list_;}

  /// \brief returns true if this path is relative
  bool IsRelative() const {return relative_flag_;}
  
  /// \brief return last path element
  String Leaf() const { return list_.empty() ? "" : list_.back(); }
  
  /// \brief return path with last element stripped away
  InfoPath Strip() const;

private:
  InfoPath(const std::vector<String> lst,bool rf);

  std::vector<String> list_;
  bool relative_flag_;

  void init(const String& path);
};

DLLEXPORT_OST_INFO std::ostream& operator<<(std::ostream&, const InfoPath& p);
#ifndef _MSC_VER
DLLEXPORT_OST_INFO InfoPath operator+(const InfoPath& p1, const InfoPath& p2);
#endif

}} // ns
