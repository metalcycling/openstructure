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
  Author: Ansgar Philippsen
*/

#include <map>

#include <ost/base.hh>
#include "module_config.hh"
#include "info_group.hh"

namespace ost { namespace info {

class DLLEXPORT_OST_INFO InfoSet {
public:
  InfoSet(const String& name="");

  // retrieve name of this set
  String GetName() const;

  // regex match of given String to all entries
  bool Match(const String& s) const;

  // add a new entry
  void Add(const String& entry);

  // clear all entries
  void Clear();

private:
  String name_;
  std::vector<String> entry_list_;
};

typedef std::map<String,InfoSet> InfoSetMap;

/*
  for each <set> in the given group, iterates
  over the <entry/> tags and assembles an
  info set; the map of all sets is then returned
*/
DLLEXPORT InfoSetMap ExtractSets(const InfoGroup&);

}} // ns
