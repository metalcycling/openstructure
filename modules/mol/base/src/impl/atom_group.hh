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
#ifndef OST_IMPL_ATOM_GROUP_HH
#define OST_IMPL_ATOM_GROUP_HH

#include <ost/mol/module_config.hh>
#include <ost/mol/impl/atom_impl_fw.hh>
#include <vector>
#include <map>
#include <ost/geom/vec3.hh>

namespace ost { namespace mol {namespace impl {

struct AtomGroupEntry {
  AtomGroupEntry() {}
  AtomGroupEntry(const AtomImplPtr& a, const geom::Vec3& p,
                 Real o, Real b)
    : atom(a), pos(p), occ(o), b_factor(b) {}
  impl::AtomImplW   atom;
  geom::Vec3        pos;
  Real              occ;
  Real              b_factor;
};

typedef std::vector<AtomGroupEntry> AtomGroupEntryList;

struct AtomGroup {
  String             name;
  AtomGroupEntryList atoms;  
};

typedef std::map<String, AtomGroup> AtomEntryGroups;

}}} // ns

#endif
