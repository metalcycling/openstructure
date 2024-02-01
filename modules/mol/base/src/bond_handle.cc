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
#include "bond_handle.hh"
#include <ost/mol/impl/connector_impl.hh>
#include "entity_visitor.hh"
#include <ost/mol/impl/atom_impl.hh>
#include <ost/invalid_handle.hh>

namespace ost {  namespace mol {

BondHandle::BondHandle():
  impl_()
{}

BondHandle::BondHandle(const impl::ConnectorImplP& im):
  impl_(im)
{}
GenericPropContainerImpl* BondHandle::GpImpl()
{
  return impl_.get();
}

const GenericPropContainerImpl* BondHandle::GpImpl() const
{
  return impl_.get();
}

BondHandle::operator bool() const 
{
  return impl_.get()!=0;
}

bool BondHandle::IsValid() const {
  return static_cast<bool>(impl_);
}

AtomHandle BondHandle::GetFirst() const
{
  this->CheckValidity();  
  AtomHandle nrvo(impl_->GetFirst());
  return nrvo;
}

AtomHandle BondHandle::GetSecond() const
{
  this->CheckValidity();  
  AtomHandle nrvo(impl_->GetSecond());
  return nrvo;
}

AtomHandle BondHandle::GetOther(const AtomHandle& a) const
{
  this->CheckValidity();
  if(a.IsValid()) {
    if(a.Impl()==impl_->GetFirst()) {
      return GetSecond();
    } else {
      return GetFirst();
    }
  }
  AtomHandle nrvo(impl_->GetFirst());
  return nrvo;
}

geom::Vec3 BondHandle::GetPos() const
{
  this->CheckValidity();  
  return impl_->GetPos();
}

geom::Vec3 BondHandle::GetOriginalPos() const
{
  this->CheckValidity();
  return impl_->GetOriginalPos();
}

Real BondHandle::GetLength() const
{
  this->CheckValidity();  
  return impl_->GetLength();
}

unsigned char BondHandle::GetBondOrder() const
{
  this->CheckValidity();
  return impl_->GetBondOrder();
}

void BondHandle::SetBondOrder(unsigned char bo)
{
  this->CheckValidity();
  impl_->SetBondOrder(bo);
}

void BondHandle::Apply(EntityVisitor& v)
{
  this->CheckValidity();
  v.VisitBond(*this);
}

void BondHandle::Apply(EntityViewVisitor& v)
{
  this->CheckValidity();
  v.VisitBond(*this);
}

unsigned long BondHandle::GetHashCode() const
{
  this->CheckValidity();
  return reinterpret_cast<unsigned long>(impl_.get());
}

bool BondExists(const AtomHandle& a, const AtomHandle& b) {
  return impl::ConnectorExists(a.Impl(), b.Impl());
}

void BondHandle::CheckValidity() const {
  if (!impl_)
    throw InvalidHandle();
}

bool BondHandle::operator==(const BondHandle& rhs) const {
  return impl_==rhs.impl_;
}

bool BondHandle::operator!=(const BondHandle& rhs) const {
  return impl_!=rhs.impl_;
}

std::ostream& operator<<(std::ostream& os, const BondHandle& bond)
{
   os << bond.GetFirst() << " - " << bond.GetSecond();
   return os;
}

}} // ns
