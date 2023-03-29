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
#include "atom_handle.hh"
#include <ost/mol/impl/atom_impl.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/mol/impl/connector_impl.hh>
#include "bond_handle.hh"
#include "entity_handle.hh"

namespace ost { namespace  mol {

AtomHandle::AtomHandle()
{}

AtomHandle::AtomHandle(const impl::AtomImplPtr& impl):
  AtomBase(impl) 
{}

EntityHandle AtomHandle::GetEntity() const 
{
  this->CheckValidity();
  return EntityHandle(Impl()->GetEntity());
}
  
ResidueHandle AtomHandle::GetResidue() const 
{
  this->CheckValidity();
  return ResidueHandle(Impl()->GetResidue());
}

int AtomHandle::GetBondCount() const 
{
  this->CheckValidity();  
  return Impl()->GetSecondaryConnectors().size()+
         (Impl()->GetPrimaryConnector() ? 1 : 0);
}

BondHandleList AtomHandle::GetBondList() const 
{
  this->CheckValidity();  
  BondHandleList l;
  l.reserve(Impl()->GetSecondaryConnectors().size());
  impl::ConnectorImplList::const_iterator i;
  for (i=Impl()->GetSecondaryConnectors().begin();
       i!=Impl()->GetSecondaryConnectors().end(); ++i) {
    l.push_back(*i);
  }
  if (Impl()->GetPrimaryConnector())
    l.push_back(Impl()->GetPrimaryConnector());
  return l;
}

void AtomHandle::Apply(EntityVisitor& v) 
{
  this->CheckValidity();  
  Impl()->Apply(v);
}

bool AtomHandle::operator==(const AtomHandle& ref) const 
{
  return Impl()==ref.Impl();
}

bool AtomHandle::operator!=(const AtomHandle& ref) const 
{
  return Impl()!=ref.Impl();
}

BondHandle AtomHandle::FindBondToAtom(const AtomHandle& handle) const 
{
  this->CheckValidity();
  return BondHandle(impl::GetConnector(handle.Impl(), Impl()));
}

AtomHandleList AtomHandle::GetBondPartners() const 
{
  this->CheckValidity();
  AtomHandleList partners;
  impl::ConnectorImplList::const_iterator i;  
  for (i=Impl()->GetSecondaryConnectors().begin();
       i!=Impl()->GetSecondaryConnectors().end(); ++i) {
    if (AtomHandle((*i)->GetFirst())!=*this) {
      partners.push_back(AtomHandle((*i)->GetFirst()));
    } else {
      partners.push_back(AtomHandle((*i)->GetSecond()));
    }
  }
  if (Impl()->GetPrimaryConnector()) {
    if (AtomHandle(Impl()->GetPrimaryConnector()->GetFirst())!=*this) {
      partners.push_back(AtomHandle(Impl()->GetPrimaryConnector()->GetFirst()));
    } else {
      partners.push_back(AtomHandle(Impl()->GetPrimaryConnector()->GetSecond()));
    }
  }
  return partners;
}

AtomHandle AtomHandle::GetHandle() const
{
  return *this;
}

unsigned long AtomHandle::GetHashCode() const
{
  this->CheckValidity();  
  return reinterpret_cast<unsigned long>(Impl().get());
}

}} // ns


