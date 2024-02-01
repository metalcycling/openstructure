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
#ifndef OST_BOND_TABLE_HH
#define OST_BOND_TABLE_HH

#include <map>

#include <ost/mol/impl/connector_impl_fw.hh>
#include <ost/mol/bond_handle.hh>

namespace ost { namespace mol {

template <typename ATOM>
struct BondTableEntry {
  BondTableEntry(const BondHandle& b, const ATOM& atom)
    : bond(b) {
      atoms[0]=atom;
  }
  BondHandle bond;  
  ATOM       atoms[2];
  bool IsComplete() const {
    return atoms[1].IsValid();
  }
};

template <typename ATOM>
struct BondTable {
  typedef BondTableEntry<ATOM> EntryType;
  typedef typename std::map<impl::ConnectorImpl*, EntryType>    MapType;
  MapType bonds;
  void Update(const BondHandle& bond, const ATOM& a) {
    if(!bond.IsValid()) return;
    typename MapType::iterator i=bonds.find(bond.Impl().get());
    if (i!=bonds.end()) {
      i->second.atoms[1]=a;
    } else {
      bonds.insert(std::make_pair(bond.Impl().get(), EntryType(bond, a)));
    }
  }
  
};

}}
#endif
