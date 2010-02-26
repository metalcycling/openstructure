//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2010 by the OpenStructure authors
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

#include <ost/string_ref.hh>

namespace ost {
  
std::pair<bool, int> StringRef::to_int() const
{
  int n=0;
  bool empty=true;
  int sig=1;
  // skip empty characters... at beginning
  for (const char* c=begin_; c!=end_; ++c) {
    if (*c=='-' && empty) {
      empty=false;
      sig=-1;
      continue;
    }
    if (*c>='0' && *c<='9') {
      empty=false;
      n=n*10+int(*c-'0');
      continue;
    }
    return std::make_pair(false, 0);
  }
  if (empty) {
    return std::make_pair(false, 0);
  }
  return std::make_pair(true, sig*n);
}

std::pair<bool, float> StringRef::to_float() const
{
  float n=0;
  bool empty=true;
  int sig=1;
  bool after_dot=false;
  float factor=0.1;
  for (const char* c=begin_; c!=end_; ++c) {
    if (*c=='-' && empty) {
      empty=false;
      sig=-1;
      continue;
    }
    if (*c=='.') {
      if (after_dot==true) {
        return std::make_pair(false, 0.0f);
      }
      after_dot=true;
      continue;
    }
    if (*c>='0' && *c<='9') {
      empty=false;
      if (after_dot==true) {
        n+=factor*int(*c-'0');
        factor*=0.1;
      } else {
        n=n*10+int(*c-'0');
      }
      continue;
    }
    return std::make_pair(false, 0.0f);
  }
  if (empty) {
    return std::make_pair(false, 0.0f);
  }
  return std::make_pair(true, sig*n);
}

std::ostream& operator<<(std::ostream& stream, const StringRef& strref)
{
  if (strref.empty()) {
    return stream;
  }
  stream.write(strref.begin(), strref.size());
  return stream;
}

}
