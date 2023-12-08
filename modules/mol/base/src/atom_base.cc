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
#include "atom_base.hh"
#include <ost/mol/impl/residue_impl.hh>
#include <ost/mol/impl/atom_impl.hh>
#include <ost/invalid_handle.hh>

namespace ost { namespace mol {

AtomBase::AtomBase()
{ }

AtomBase::AtomBase(const impl::AtomImplPtr& impl):
  impl_(impl) 
{}
GenericPropContainerImpl* AtomBase::GpImpl()
{
  return Impl().get();
}

const GenericPropContainerImpl* AtomBase::GpImpl() const
{
  return Impl().get();
}

const String& AtomBase::GetName() const 
{
  this->CheckValidity();  
  return impl_->Name();
}

void AtomBase::SetName(const String& atom_name) 
{
  this->CheckValidity();
  impl_->Name()=atom_name;
}

const geom::Vec3& AtomBase::GetPos() const 
{
  this->CheckValidity();  
  return impl_->TransformedPos();
}

const geom::Vec3& AtomBase::GetOriginalPos() const
{
  this->CheckValidity();
  return impl_->OriginalPos();
}

geom::Vec3 AtomBase::GetAltPos(const String& alt_group) const
{
  this->CheckValidity();
  return impl_->GetResidue()->GetAltAtomPos(Impl(), alt_group);
}

Real AtomBase::GetAltBFactor(const String& alt_group) const
{
  this->CheckValidity();
  return impl_->GetResidue()->GetAltAtomBFactor(Impl(), alt_group);
}
Real AtomBase::GetAltOcc(const String& alt_group) const
{
  this->CheckValidity();
  return impl_->GetResidue()->GetAltAtomOcc(Impl(), alt_group);
}

std::vector<String> AtomBase::GetAltGroupNames() const
{
  this->CheckValidity();
  return impl_->GetResidue()->GetAltAtomGroupNames(Impl());
}

impl::AtomImplPtr& AtomBase::Impl() 
{
  return impl_;
}

const impl::AtomImplPtr& AtomBase::Impl() const 
{
  return impl_;
}

String AtomBase::GetQualifiedName() const 
{
  this->CheckValidity();
  return Impl()->GetQualifiedName();
}

bool AtomBase::IsValid() const
{
  return Impl().get()!=0 && impl_->GetEntity();
}

void AtomBase::CheckValidity() const 
{
  if (! IsValid())
    throw InvalidHandle();
}

std::ostream& operator<<(std::ostream& os, const AtomBase& atom) 
{
  if (atom.Impl()) {
    os << atom.GetQualifiedName();
  } else {
    os << "invalid atom";
  }
  return os;
}


Real AtomBase::GetRadius() const
{
  this->CheckValidity();
  return Impl()->GetRadius();
}


const String& AtomBase::GetElement() const
{
  this->CheckValidity();
  return Impl()->GetElement();
}


bool AtomBase::IsHetAtom() const
{
  this->CheckValidity();
  return Impl()->IsHetAtom();
}

Real AtomBase::GetBFactor() const
{
  this->CheckValidity();
  return Impl()->GetBFactor();
}
void AtomBase::SetElement(const String& element)
{
  this->CheckValidity();
  Impl()->SetElement(element);
}

void AtomBase::SetBFactor(Real factor)
{
  this->CheckValidity();
  Impl()->SetBFactor(factor);
}

void AtomBase::SetOccupancy(Real occ)
{
  this->CheckValidity();
  Impl()->SetOccupancy(occ);
}

Real AtomBase::GetMass() const
{
  this->CheckValidity();
  return Impl()->GetMass();
}

Real AtomBase::GetCharge() const
{
  this->CheckValidity();
  return Impl()->GetCharge();
}

Real AtomBase::GetOccupancy() const
{
  this->CheckValidity();
  return Impl()->GetOccupancy();
}

void AtomBase::SetCharge(Real charge)
{
  this->CheckValidity();
  Impl()->SetCharge(charge);
}

void AtomBase::SetRadius(Real radius)
{
  this->CheckValidity();
  Impl()->SetRadius(radius);
}

const geom::Mat3& AtomBase::GetAnisou() const
{
  this->CheckValidity();
  return Impl()->GetAnisou();
}

void AtomBase::SetAnisou(const geom::Mat3& anisou)
{
  this->CheckValidity();
  Impl()->SetAnisou(anisou);
}

void AtomBase::SetMass(Real mass)
{
  this->CheckValidity();
  Impl()->SetMass(mass);
}
void AtomBase::SetHetAtom(bool het)
{
  this->CheckValidity();
  Impl()->SetHetAtom(het);
}

String AtomBase::GetStringProperty(Prop::ID prop_id) const
{
  this->CheckValidity();
  return Impl()->GetStringProperty(prop_id);
}
                      

Real AtomBase::GetFloatProperty(Prop::ID prop_id) const
{
  this->CheckValidity();
  return Impl()->GetFloatProperty(prop_id);
}

int AtomBase::GetIntProperty(Prop::ID prop_id) const
{
  this->CheckValidity();
  return Impl()->GetIntProperty(prop_id);
}

unsigned long AtomBase::GetIndex() const
{
  this->CheckValidity();
  return Impl()->GetIndex();
}

void AtomBase::SetIndex(const unsigned long index) 
{
  this->CheckValidity();
  return impl_->SetIndex(index);
}

  
}} // ns

