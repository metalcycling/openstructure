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
#include "entity_base.hh"
#include <ost/mol/impl/entity_impl.hh>
#include <ost/invalid_handle.hh>

namespace ost { namespace  mol {

EntityBase::EntityBase()
{ }

EntityBase::EntityBase(const impl::EntityImplPtr& impl):
  impl_(impl) {
}

impl::EntityImplPtr& EntityBase::Impl() { 
  return impl_; 
}

const impl::EntityImplPtr& EntityBase::Impl() const {
  return impl_;
}
GenericPropContainerImpl* EntityBase::GpImpl()
{
  return Impl().get();
}

const GenericPropContainerImpl* EntityBase::GpImpl() const
{
  return Impl().get();
}

void EntityBase::CheckValidity() const {
  if (!impl_)
    throw InvalidHandle();
}

const String& EntityBase::GetName() const {
  this->CheckValidity();
  return impl_->GetName();
}

void EntityBase::SetName(const String& ent_name) {
  this->CheckValidity();
  impl_->SetName(ent_name);
}

DLLEXPORT_OST_MOL std::ostream& operator<<(std::ostream& os,   
                                           const EntityBase& ent) {
  if (ent.Impl()) {
    os << "entity";
  } else {
    os << "invalid";
  }

  return os;
}
}} //ns

