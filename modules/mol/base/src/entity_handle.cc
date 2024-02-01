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

#include <ost/log.hh>

#include "impl/entity_impl.hh"

#include "atom_view.hh"
#include "bond_handle.hh"
#include "torsion_handle.hh"
#include "entity_visitor.hh"
#include "atom_handle.hh"
#include "residue_handle.hh"
#include "chain_handle.hh"
#include "entity_handle.hh"
#include "impl/chain_impl.hh"
#include "impl/atom_impl.hh"
#include "impl/residue_impl.hh"
#include "xcs_editor.hh"
#include "ics_editor.hh"

namespace ost { namespace mol {


EntityHandle::EntityHandle() 
{
}

EntityHandle CreateEntity() {
  return EntityHandle(impl::EntityImplPtr(new impl::EntityImpl()));
}

geom::Vec3 EntityHandle::GetCenterOfMass() const {
  this->CheckValidity();
  return Impl()->GetCenterOfMass();
}

Real EntityHandle::GetMass() const {
  this->CheckValidity();
  return Impl()->GetMass();
}

geom::Vec3 EntityHandle::GetCenterOfAtoms() const {
  this->CheckValidity();
  return Impl()->GetCenterOfAtoms();
}

geom::AlignedCuboid EntityHandle::GetBounds() const 
{
  this->CheckValidity();
  return Impl()->GetBounds();
}


EntityHandle::EntityHandle(const impl::EntityImplPtr& e)
  : EntityBase(e)
{}

EntityHandle EntityHandle::Copy() const {
  this->CheckValidity();
  return EntityHandle(Impl()->Copy());
}

void EntityHandle::Swap(EntityHandle& eh) {
  this->CheckValidity();
  impl::EntityImplPtr& other=eh.Impl();
  Impl().swap(other);
}

void EntityHandle::DeepSwap(EntityHandle& eh) {
  this->CheckValidity();
  Impl()->Swap(*eh.Impl());
  Impl()->NotifyObserver();
  eh.Impl()->NotifyObserver();
}

void EntityHandle::AttachObserver(const EntityObserverPtr& o) {
  this->CheckValidity();
  Impl()->AttachObserver(o);
}

void EntityHandle::DetachObserver(const EntityObserverPtr& o) {
  this->CheckValidity();
  Impl()->DetachObserver(o);
}


void EntityHandle::Apply(EntityVisitor& v)
{
  this->CheckValidity();
  Impl()->Apply(v);
}

AtomHandleList EntityHandle::FindWithin(const geom::Vec3& pos,
                                        Real radius) const {
  this->CheckValidity();
  impl::AtomImplList atoms=Impl()->FindWithin(pos,radius);
  AtomHandleList handles;
  handles.reserve(atoms.size());
  for (size_t i=0; i<atoms.size(); ++i) {
    handles.push_back(atoms[i]);
  }
  return handles;
}

void EntityHandle::SetDefaultQueryFlags(QueryFlags flags)
{
  this->CheckValidity();
  Impl()->SetDefaultQueryFlags(flags);
}

QueryFlags EntityHandle::GetDefaultQueryFlags() const
{
  this->CheckValidity();
  return Impl()->GetDefaultQueryFlags();
}

EntityView EntityHandle::Select(const Query& q) const
{
  this->CheckValidity();
  return Impl()->Select(*this, q);
}

EntityView EntityHandle::Select(const String& q) const
{
  this->CheckValidity();
  return Impl()->Select(*this, Query(q));
}

EntityView EntityHandle::Select(const Query& q, QueryFlags flags) const
{
  this->CheckValidity();
  return Impl()->Select(*this, q, flags);
}

EntityView EntityHandle::Select(const String& q, QueryFlags flags) const 
{
  this->CheckValidity();
  return Impl()->Select(*this, Query(q), flags);
}

EntityView EntityHandle::CreateFullView() const
{
  this->CheckValidity();
  return Impl()->CreateFullView(*this);
}

bool EntityHandle::operator==(const EntityHandle& ref) const
{
  return Impl()==ref.Impl();
}

bool EntityHandle::operator!=(const EntityHandle& ref) const
{
  return Impl()!=ref.Impl();
}

ChainHandleList EntityHandle::GetChainList() const
{
  this->CheckValidity();
  ChainHandleList chainlist;
  const impl::ChainImplList& cl = Impl()->GetChainList();
  impl::ChainImplList::const_iterator i;
  chainlist.reserve(cl.size());
  for (i=cl.begin(); i!=cl.end(); ++i) {
    chainlist.push_back(ChainHandle(*i));
  }
  return chainlist;
}

TorsionHandle EntityHandle::FindTorsion(const AtomHandle& a1,
                                        const AtomHandle& a2,
                                        const AtomHandle& a3,
                                        const AtomHandle& a4) const {
  this->CheckValidity();
  return Impl()->FindTorsion(a1.Impl(), a2.Impl(), a3.Impl(), a4.Impl());
}

ChainHandle EntityHandle::FindChain(const String& name) const {
  this->CheckValidity();
  return Impl()->FindChain(name);
}

EntityView EntityHandle::CreateEmptyView() const {
  this->CheckValidity();
  return EntityView(*this);
}

int EntityHandle::GetAtomCount() const {
  this->CheckValidity();
  return Impl()->GetAtomCount();
}

mol::BondHandleList EntityHandle::GetBondList() const {
  this->CheckValidity();
  return Impl()->GetBondList();
}

int EntityHandle::GetBondCount() const {
  this->CheckValidity();
  return Impl()->GetBondCount();
}

int EntityHandle::GetResidueCount() const {
  this->CheckValidity();
  return Impl()->GetResidueCount();
}

int EntityHandle::GetChainCount() const {
  this->CheckValidity();
  return Impl()->GetChainCount();
}

Real EntityHandle::GetAngle(const AtomHandle& a1, const AtomHandle& a2,
                              const AtomHandle& a3) const {
  this->CheckValidity();
  return Impl()->GetAngle(a1.Impl(), a2.Impl(), a3.Impl());
}

Real EntityHandle::GetAngle(const AtomView& a1, const AtomView& a2,
                              const AtomView& a3) const
{
  return GetAngle(a1.GetHandle(), a2.GetHandle(), a3.GetHandle());
}

geom::Mat4 EntityHandle::GetTransformationMatrix() const
{
  static bool warn=true;
  if(warn) {
    LOG_WARNING("Entity::GetTransformationMatrix is deprecated, use GetTransform instead");
    warn=false;
  }
  this->CheckValidity();
  return Impl()->GetTransform().GetMatrix();
}


geom::Mat4 EntityHandle::GetInvTransformationMatrix() const
{
  static bool warn=true;
  if(warn) {
    LOG_WARNING("Entity::GetInvTransformationMatrix is deprecated, use GetTransform instead");
    warn=false;
  }
  this->CheckValidity();
  return Impl()->GetTransform().GetInvertedMatrix();
}

bool EntityHandle::IsTransformationIdentity() const
{
  static bool warn=true;
  if(warn) {
    LOG_WARNING("Entity::IsTransformationIdentity is deprecated, use HasTransform instead");
    warn=false;
  }
  this->CheckValidity();
  return !Impl()->HasTransform();
}

geom::Transform EntityHandle::GetTransform() const
{
  this->CheckValidity();
  return Impl()->GetTransform();  
}

void EntityHandle::SetTransform(const geom::Transform& tf)
{
  this->CheckValidity();
  Impl()->SetTransform(tf);  
}

bool EntityHandle::HasTransform() const
{
  this->CheckValidity();
  return Impl()->HasTransform();  
}

void EntityHandle::ClearTransform()
{
  this->CheckValidity();
  Impl()->ClearTransform();  
}

void EntityHandle::FixTransform()
{
  this->CheckValidity();
  Impl()->FixTransform();  
}

ResidueHandle EntityHandle::FindResidue(const String& chain_name,
                                        const ResNum& residue) const {
  this->CheckValidity();
  return Impl()->FindResidue(chain_name, residue);
}

AtomHandle EntityHandle::FindAtom(const String& chain_name,
                                  const ResNum& residue,
                                  const String& atom_name) const {
  this->CheckValidity();
  return Impl()->FindAtom(chain_name, residue, atom_name);
}


XCSEditor EntityHandle::EditXCS(EditMode mode) const
{
  this->CheckValidity();
  return XCSEditor(*this, mode);
}

ICSEditor EntityHandle::EditICS(EditMode mode) const
{
  this->CheckValidity();
  Impl()->EnableICS();
  return ICSEditor(*this, mode);
}


ResidueHandleList EntityHandle::GetResidueList() const
{
  this->CheckValidity();
  ResidueHandleList residues;
  residues.reserve(this->GetResidueCount());
  for (impl::ChainImplList::const_iterator i = Impl()->GetChainList().begin(),
       e = Impl()->GetChainList().end(); i != e; ++i) {
    for (impl::ResidueImplList::const_iterator j = (*i)->GetResidueList().begin(),
         e2 = (*i)->GetResidueList().end(); j != e2; ++j) {
      residues.push_back(ResidueHandle(*j));
    }
  }
  return residues;
}


AtomHandleList EntityHandle::GetAtomList() const
{
  this->CheckValidity();
  AtomHandleList atoms;
  atoms.reserve(this->GetAtomCount());
  for (impl::ChainImplList::const_iterator i = Impl()->GetChainList().begin(),
       e = Impl()->GetChainList().end(); i != e; ++i) {
    for (impl::ResidueImplList::const_iterator j = (*i)->GetResidueList().begin(),
         e2 = (*i)->GetResidueList().end(); j != e2; ++j) {
      for (impl::AtomImplList::const_iterator k = (*j)->GetAtomList().begin(),
           e3 = (*j)->GetAtomList().end(); k != e3; ++k) {
        atoms.push_back(AtomHandle(*k));
      }
    }
  }
  return atoms;
}

bool less_index(const mol::AtomHandle& a1, const mol::AtomHandle& a2)
{
  return a1.GetIndex()<a2.GetIndex();
}

geom::Vec3List EntityHandle::GetAtomPosList(bool ordered_by_index) const {
  this->CheckValidity();
  geom::Vec3List atom_pos_list;
  atom_pos_list.reserve(this->GetAtomCount());
  AtomHandleList atom_list=this->GetAtomList();
  if (ordered_by_index){
    std::sort(atom_list.begin(),atom_list.end(),less_index);
  }
  for (AtomHandleList::const_iterator a=atom_list.begin(), e=atom_list.end(); a!=e; ++a) {
    atom_pos_list.push_back(a->GetPos());
  }
  return atom_pos_list;
}
  
EntityHandle EntityHandle::GetHandle() const
{
  return *this;
}

unsigned long EntityHandle::GetHashCode() const
{
  this->CheckValidity();
  return reinterpret_cast<unsigned long>(Impl().get());
}

}} // ns
