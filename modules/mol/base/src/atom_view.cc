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

/*
   Author: Marco Biasini
 */

#include "atom_view.hh"
#include "residue_view.hh"
#include <ost/mol/bond_handle.hh>
#include <ost/mol/atom_handle.hh>
#include <ost/mol/impl/atom_impl.hh>
#include <ost/mol/entity_visitor.hh>
#include <ost/mol/entity_view.hh>
#include <ost/mol/chain_view.hh>

namespace ost { namespace  mol {

class DLLEXPORT_OST_MOL AtomViewData {
public:
  AtomViewData() {}
  AtomViewData(const ResidueViewDataPtr& res)
    : residue(res) {
  }
  ResidueViewDataWeakPtr residue;
  BondHandleList         bonds;
};


AtomView::AtomView() 
{

}

bool AtomView::operator==(const AtomView& rhs) const 
{
  return data_==rhs.data_;
}

bool AtomView::operator!=(const AtomView& rhs) const 
{
  return data_!=rhs.data_;
}

AtomView::AtomView(const ResidueView& residue_view,
                   const AtomHandle& atom)
 : AtomBase(atom.Impl()), data_(new AtomViewData(residue_view.ViewData())) 
{
  this->GetEntity().AddAtomInternal(*this);
}

ResidueView AtomView::GetResidue() const 
{
  this->CheckValidity();
  if (!data_->residue.expired()) {
    return ResidueView(data_->residue.lock(), Impl()->GetResidue());    
  }
  throw InvalidHandle();
}

void AtomView::Apply(EntityVisitor& visitor) 
{
  this->CheckValidity();
  visitor.VisitAtom(this->GetHandle());
}

void AtomView::Apply(EntityViewVisitor& visitor) 
{
  this->CheckValidity();
  visitor.VisitAtom(*this);
}

AtomHandle AtomView::GetHandle() const 
{
  return AtomHandle(Impl());
}

bool AtomView::AddBond(const BondHandle& bond) 
{
  this->CheckValidity();
  if (std::count(data_->bonds.begin(), data_->bonds.end(), bond))
    return false;
  data_->bonds.push_back(bond);
  return true;
}

BondHandleList AtomView::GetBondList() const 
{
  this->CheckValidity();
  return data_->bonds;
}

int AtomView::GetBondCount() const 
{
  this->CheckValidity();
  return data_->bonds.size();
}


mol::AtomViewList AtomView::GetBondPartners() const
{
  this->CheckValidity();
  mol::AtomViewList avl;
  mol::BondHandleList::const_iterator i;
  for (i=data_->bonds.begin();i!=data_->bonds.end();++i) {
    if (i->GetFirst().GetHandle()!=this->GetHandle()) {
      avl.push_back(this->GetEntity().ViewForHandle(i->GetFirst()));
    } else {
      avl.push_back(this->GetEntity().ViewForHandle(i->GetSecond()));
    }
  }
  return avl;
}

void AtomView::RemoveBonds() 
{
  this->CheckValidity();
  EntityView ent=this->GetEntity();
  // copy bonds into a new list 
  BondHandleList bonds=data_->bonds;
  for (BondHandleList::iterator i=bonds.begin(); i!=bonds.end(); ++i) {
    BondHandle b=*i;
    AtomView av;
    if (b.GetFirst()==this->GetHandle()) {
      av=ent.ViewForHandle(b.GetSecond());
    } else {
      av=ent.ViewForHandle(b.GetFirst());
    }
    ent.RemoveBond(b);
  }
  data_->bonds.clear();
}

EntityView AtomView::GetEntity() const 
{
  return this->GetResidue().GetChain().GetEntity();
}

void AtomView::RemoveBondInternal(const BondHandle& bond)
{
  BondHandleList::iterator i=std::find(data_->bonds.begin(), 
                                       data_->bonds.end(), bond);
  if (i!=data_->bonds.end()) {
    data_->bonds.erase(i);    
  }
}

unsigned long AtomView::GetHashCode() const
{
  this->CheckValidity();  
  return reinterpret_cast<unsigned long>(data_.get());
} 

bool AtomView::IsValid() const
{
  return Impl().get()!=0 && Impl()->GetEntity();
}

}} // ns


