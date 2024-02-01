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
#include <ost/mol/residue_handle.hh>
#include <ost/mol/entity_visitor.hh>
#include <ost/mol/entity_handle.hh>
#include <ost/mol/chain_handle.hh>
#include <ost/mol/torsion_handle.hh>
#include <ost/mol/impl/residue_impl.hh>
#include <ost/mol/impl/chain_impl.hh>
#include <ost/mol/impl/entity_impl.hh>
#include <ost/mol/bond_handle.hh>

namespace ost { namespace mol {

ResidueHandle::ResidueHandle()
  : ResidueBase() {
}

EntityHandle ResidueHandle::GetEntity() const 
{
  this->CheckValidity();
  return EntityHandle(Impl()->GetEntity());
}

ChainHandle ResidueHandle::GetChain() const 
{
  this->CheckValidity();
  return ChainHandle(Impl()->GetChain());
}

void ResidueHandle::Apply(EntityVisitor& v)
{
  this->CheckValidity();
  Impl()->Apply(v);
}

ResidueHandle ResidueHandle::GetPrev() const 
{
  this->CheckValidity();
  return Impl()->GetChain()->GetPrev(Impl());
}
TorsionHandle ResidueHandle::FindTorsion(const String& torsion_name) const 
{
  this->CheckValidity();
  return TorsionHandle(Impl()->FindTorsion(torsion_name));
}

ResidueHandle ResidueHandle::GetNext() const 
{
  this->CheckValidity();
  return Impl()->GetChain()->GetNext(Impl());
}

AtomHandle ResidueHandle::FindAtom(const String& aname) const
{
  this->CheckValidity();
  return AtomHandle(Impl()->FindAtom(aname));
}

int ResidueHandle::GetAtomCount() const 
{
  this->CheckValidity();
  return Impl()->GetAtomCount();
}

int ResidueHandle::GetBondCount() const 
{
  this->CheckValidity();
  return Impl()->GetBondCount();
}

AtomHandleList ResidueHandle::GetAtomList() const
{
  this->CheckValidity();
  AtomHandleList atomlist;
  const impl::AtomImplList& al = Impl()->GetAtomList();
  for (impl::AtomImplList::const_iterator it=al.begin();
       it!=al.end();++it) {
    atomlist.push_back(AtomHandle(*it));
  }
  return atomlist;

}




const String& ResidueHandle::GetCurrentAltGroupName() const 
{
  this->CheckValidity();
  return Impl()->GetCurrentAltGroupName();
}

TorsionHandleList ResidueHandle::GetTorsionList() const 
{
  this->CheckValidity();
  const impl::TorsionImplList& impl_list=Impl()->GetTorsionList();
  TorsionHandleList torsions;
  torsions.reserve(impl_list.size());
  impl::TorsionImplList::const_iterator src=impl_list.begin();
  for (; src!=impl_list.end(); ++src) {
    torsions.push_back(TorsionHandle(*src));
  }
  return torsions;
}

bool ResidueHandle::operator==(const ResidueHandle& ref) const
{
  return Impl()==ref.Impl();
}

bool ResidueHandle::operator!=(const ResidueHandle& ref) const
{
  return Impl()!=ref.Impl();
}

int ResidueHandle::GetIndex() const 
{
  this->CheckValidity();
  return Impl()->GetIndex();
}

bool ResidueHandle::HasAltAtomGroup(const String& group) const 
{
  this->CheckValidity();
  return Impl()->HasAltAtomGroup(group);
}


std::vector<String> ResidueHandle::GetAltAtomGroupNames() const 
{
  this->CheckValidity();
  return Impl()->GetAltAtomGroupNames();
}


bool ResidueHandle::HasAltAtoms() const 
{
  this->CheckValidity();
  return Impl()->HasAltAtoms();
}

bool ResidueHandle::SwitchAtomPos(const String& group) 
{
  this->CheckValidity();
  return Impl()->SwitchAtomPos(group);
}


bool InSequence(const ResidueHandle& r1, const ResidueHandle& r2) 
{
  if (!r1 || !r2) return false; // protect against invalid handles
  if (r1.GetNext()!=r2) {
    return false;
  }
  if (r1.IsPeptideLinking() && r2.IsPeptideLinking()) {
    return (BondExists(r1.FindAtom("C"), r2.FindAtom("N")));    
  } else if (r1.IsNucleotideLinking() && r2.IsNucleotideLinking()) {
    return BondExists(r1.FindAtom("O3'"), r2.FindAtom("P"));
  }
  return false;
}

geom::Vec3 ResidueHandle::GetAltAtomPos(const AtomHandle& atom, 
                                        const String& group) const
{
  this->CheckValidity();
  return Impl()->GetAltAtomPos(atom.Impl(), group);
}

geom::AlignedCuboid ResidueHandle::GetBounds() const
{
  this->CheckValidity();
  return Impl()->GetBounds();
}

geom::Vec3 ResidueHandle::GetCenterOfMass() const
{
  this->CheckValidity();
  return Impl()->GetCenterOfMass();
}

geom::Vec3 ResidueHandle::GetCenterOfAtoms() const
{
  this->CheckValidity();
  return Impl()->GetCenterOfAtoms();
}

double ResidueHandle::GetMass() const
{
  this->CheckValidity();
  return Impl()->GetMass();
}

unsigned long ResidueHandle::GetHashCode() const
{
  this->CheckValidity();
  return reinterpret_cast<unsigned long>(Impl().get());
}

ResidueHandle ResidueHandle::GetHandle() const
{
  return *this;
}

EntityView ResidueHandle::Select(const Query& q, QueryFlags flags) const
{
  this->CheckValidity();
  std::stringstream s;
  if (q.GetQueryString() != "") {
    s << "cname='" << Impl()->GetChain()->GetName() << "' and rindex=" << Impl()->GetIndex()
      << " and (" << q.GetQueryString() << ")";
  } else {
    s << "cname='" << Impl()->GetChain()->GetName() << "' and rindex=" 
      << Impl()->GetIndex();
  }
  return this->GetEntity().Select(Query(s.str()), flags);
}

EntityView ResidueHandle::Select(const String& q, QueryFlags flags) const {
  this->CheckValidity();
  std::stringstream s;
  if (q != "") {
     s << "cname='" << Impl()->GetChain()->GetName() << "' and rindex=" << Impl()->GetIndex()
       << " and (" << q << ")";
  }
  else s << "cname='" << Impl()->GetChain()->GetName() << "' and rindex=" << Impl()->GetIndex();
  return this->GetEntity().Select(Query(s.str()), flags);
}

}} //ns

