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
#include "residue_base.hh"
#include <ost/mol/impl/residue_impl.hh>
#include <ost/mol/torsion_handle.hh>
#include <ost/invalid_handle.hh>

namespace ost { namespace mol {

ResidueBase::ResidueBase()
{ }

ResidueBase::ResidueBase(const ResidueBase& rhs):
  impl_(rhs.impl_) 
{}

ResidueBase::ResidueBase(const impl::ResidueImplPtr& impl):
  impl_(impl) 
{}

GenericPropContainerImpl* ResidueBase::GpImpl()
{
  return Impl().get();
}

const GenericPropContainerImpl* ResidueBase::GpImpl() const
{
  return Impl().get();
}

AtomHandle ResidueBase::GetCentralAtom() const 
{
  this->CheckValidity();
  return AtomHandle(Impl()->GetCentralAtom());
}

void ResidueBase::SetCentralAtom(const AtomHandle& a)
{
  this->CheckValidity();
  impl_->SetCentralAtom(a.Impl());
}

geom::Vec3 ResidueBase::GetCentralNormal() const 
{
  this->CheckValidity();
  return impl_->GetCentralNormal();
}

const ResNum& ResidueBase::GetNumber() const {
  this->CheckValidity();
  return impl_->GetNumber();
}

const ResidueKey& ResidueBase::GetKey() const {
  this->CheckValidity();
  return impl_->GetKey();
}

const String& ResidueBase::GetName() const {
  this->CheckValidity();
  return impl_->GetName();
}

String ResidueBase::GetQualifiedName() const {
  this->CheckValidity();
  return impl_->GetQualifiedName();
}

bool ResidueBase::IsPeptideLinking() const 
{
  return this->GetChemClass().IsPeptideLinking();
}

bool ResidueBase::IsNucleotideLinking() const 
{
  return this->GetChemClass().IsNucleotideLinking();
}

ChemClass ResidueBase::GetChemClass() const 
{
  this->CheckValidity();
  return Impl()->GetChemClass();
};

void ResidueBase::SetChemClass(ChemClass cc) 
{
  this->CheckValidity();
  Impl()->SetChemClass(cc);
}

void ResidueBase::SetChemType(ChemType ct)
{
  this->CheckValidity();
  Impl()->SetChemType(ct);
}

ChemType ResidueBase::GetChemType() const
{
  this->CheckValidity();
  return Impl()->GetChemType();
};

SecStructure ResidueBase::GetSecStructure() const 
{
  this->CheckValidity();
  return impl_->GetSecStructure();
}

void ResidueBase::SetSecStructure(SecStructure ss) 
{
  this->CheckValidity();
  impl_->SetSecStructure(ss);
}

TorsionHandle ResidueBase::GetPhiTorsion() const 
{
  this->CheckValidity();
  return impl_->GetPhiTorsion();
}

TorsionHandle ResidueBase::GetPsiTorsion() const 
{
  this->CheckValidity();
  return impl_->GetPsiTorsion();
}

TorsionHandle ResidueBase::GetOmegaTorsion() const
{
  this->CheckValidity();
  return impl_->GetOmegaTorsion();
}

impl::ResidueImplPtr& ResidueBase::Impl() 
{
  return impl_;
}

const impl::ResidueImplPtr& ResidueBase::Impl() const 
{
  return impl_;
}

bool ResidueBase::IsValid() const
{
  return Impl().get()!=0;
}

void ResidueBase::CheckValidity() const 
{
  if (! (impl_ && impl_->GetEntity()))
    throw InvalidHandle();
}

std::ostream& operator<<(std::ostream& os, const ResidueBase& residue) {
  if (residue.Impl()) {
    os << residue.GetQualifiedName();
  } else {
    os << "invalid residue"; 
  }
   return os;                                                                              
}

char ResidueBase::GetOneLetterCode() const 
{
  this->CheckValidity();
  return Impl()->GetOneLetterCode();
}

void ResidueBase::SetOneLetterCode(char olc)  
{
  this->CheckValidity();
  Impl()->SetOneLetterCode(olc);
}

String ResidueBase::GetStringProperty(Prop::ID prop_id) const
{
  this->CheckValidity();
  return Impl()->GetStringProperty(prop_id);
}
                      

Real ResidueBase::GetFloatProperty(Prop::ID prop_id) const
{
  this->CheckValidity();
  return Impl()->GetFloatProperty(prop_id);
}

int ResidueBase::GetIntProperty(Prop::ID prop_id) const
{
  this->CheckValidity();
  return Impl()->GetIntProperty(prop_id);
}

bool ResidueBase::IsProtein() const
{
  this->CheckValidity();
  return Impl()->IsProtein();
}

void ResidueBase::SetIsProtein(bool protein)
{
  this->CheckValidity();
  Impl()->SetProtein(protein);
}

void ResidueBase::SetIsLigand(bool ligand)
{
  this->CheckValidity();
  return Impl()->SetIsLigand(ligand);
}

bool ResidueBase::IsLigand() const
{
  this->CheckValidity();
  return Impl()->IsLigand();
}

}}

