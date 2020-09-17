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
#include "chain_base.hh"
#include <ost/mol/impl/chain_impl.hh>
#include <ost/invalid_handle.hh>

namespace ost { namespace mol {

ChainBase::ChainBase()
{ }

ChainBase::ChainBase(const impl::ChainImplPtr& impl): 
  impl_(impl) 
{}
GenericPropContainerImpl* ChainBase::GpImpl()
{
  return impl_.get();
}


const GenericPropContainerImpl* ChainBase::GpImpl() const
{
  return impl_.get();
}

String ChainBase::GetName() const {
  this->CheckValidity();
  return impl_->GetName();
}

ChainType ChainBase::GetType() const {
  return impl_->GetType();
}

String ChainBase::GetDescription() const {
  return impl_->GetDescription();
}

void ChainBase::CheckValidity() const {
  if (!impl_)
    throw InvalidHandle();
}

std::ostream& operator<<(std::ostream& os, const ChainBase& chain) 
{
  if (chain.Impl()) {
    os << chain.GetName();
  } else {
    os << "invalid chain";
  }
  return os;
}

bool ChainBase::IsPolymer() const
{
  this->CheckValidity();
  return impl_->IsPolymer();
  
}

bool ChainBase::IsPolysaccharide() const
{
  this->CheckValidity();
  return impl_->IsPolysaccharide();
  
}

bool ChainBase::IsOligosaccharide() const
{
  this->CheckValidity();
  return impl_->IsOligosaccharide();

}

bool ChainBase::IsPolypeptide() const
{
  this->CheckValidity();
  return impl_->IsPolypeptide();
  
}

bool ChainBase::IsPolynucleotide() const
{
  this->CheckValidity();
  return impl_->IsPolynucleotide();
}

}} // ns

