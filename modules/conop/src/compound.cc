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
#include "compound.hh"

namespace ost { namespace conop {


int Compound::GetAtomSpecIndex(const String& name) const {
  AtomSpecList::const_iterator i;
  // BZDNG-261: first search all primary atom names before falling back to
  // alternative names. There are some files where alternative atom names are 
  // used as primary names for other atoms
  for (i=atom_specs_.begin(); i!=atom_specs_.end(); ++i) {
    if ((*i).name==name)
      return std::distance(atom_specs_.begin(), i);
  } 
  for (i=atom_specs_.begin(); i!=atom_specs_.end(); ++i) {
    if ((*i).alt_name==name)
      return std::distance(atom_specs_.begin(), i);
  }
 return -1;
}
String Date::ToString() const
{
  std::stringstream ss;
  ss << year << "-";
  ss.fill('0');
  ss.width(2);
  ss << month << "-";
  ss.fill('0');
  ss.width(2);
  ss << day;
  return ss.str();
}
}}
