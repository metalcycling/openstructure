//------------------------------------------------------------------------------
//
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
#include <algorithm>

#include <boost/bind/bind.hpp>
#include <boost/logic/tribool.hpp>

#include <ost/log.hh>

#include "entity_view.hh"
#include "entity_handle.hh"
#include "chain_view.hh"
#include "entity_visitor.hh"
#include "residue_view.hh"
#include "atom_view.hh"
#include "bond_table.hh"
#include "query_state.hh"
#include "entity_property_mapper.hh"

#include "impl/entity_impl.hh"
#include "impl/query_impl.hh"
#include "impl/residue_impl.hh"
#include "impl/chain_impl.hh"

using boost::logic::tribool;
using boost::logic::indeterminate;



namespace ost { namespace mol {

typedef std::map<unsigned long, AtomView> HandleToViewMap;

class DLLEXPORT_OST_MOL EntityViewData {
public:
  ChainViewList           chains;
  BondHandleList          bonds;
  HandleToViewMap         handle_to_view;
  
  AtomView ViewForHandle(const AtomHandle& atom)
  {
    HandleToViewMap::iterator i=handle_to_view.find(atom.GetHashCode());
    if (i!=handle_to_view.end()) {
      return i->second;
    }
    return AtomView();
  }
};

EntityView::EntityView() {

}

EntityView::EntityView(const EntityHandle& entity,
                       ViewAddFlags flags)
  : EntityBase(entity.Impl()), data_(new EntityViewData()) {
  if (flags & ViewAddFlag::INCLUDE_CHAINS) {
    impl::ChainImplList::const_iterator i=entity.Impl()->GetChainList().begin();
    for (; i!=entity.Impl()->GetChainList().end(); ++i) {
      this->AddChain(ChainHandle(*i), flags);
    }
  }
}

EntityView::EntityView(const EntityViewDataPtr& data,
                       const impl::EntityImplPtr impl):
 EntityBase(impl), data_(data) 
{}

const EntityViewDataPtr& EntityView::ViewData() const 
{
  return data_;
}

EntityViewDataPtr& EntityView::ViewData() 
{
  return data_;
}
void EntityView::Apply(EntityVisitor& visitor) 
{
  this->CheckValidity();
  ChainViewList::iterator i;
  visitor.OnEntry();
  for (i=data_->chains.begin(); i!=data_->chains.end(); ++i) {
    (*i).Apply(visitor);
  }
  BondHandleList::iterator j;
  for(j=data_->bonds.begin();j!=data_->bonds.end();++j) {
    (*j).Apply(visitor);
  }
  visitor.OnExit();
}

void EntityView::Apply(EntityViewVisitor& visitor) 
{
  this->CheckValidity();
  ChainViewList::iterator i;
  visitor.OnEntry();
  for (i=data_->chains.begin(); i!=data_->chains.end(); ++i) {
    (*i).Apply(visitor);
  }
  BondHandleList::iterator j;
  for(j=data_->bonds.begin();j!=data_->bonds.end();++j) 
  {
    (*j).Apply(visitor);
  }
  visitor.OnExit();
}

EntityHandle EntityView::GetHandle() const 
{
  return EntityHandle(Impl());
}

int EntityView::GetAtomCount() const 
{
  this->CheckValidity();
  ChainViewList::const_iterator it=data_->chains.begin();
  int count=0;
  while(it!=data_->chains.end()) {
    count+=(*it).GetAtomCount();
    ++it;
  }
  return count;
}

int EntityView::GetResidueCount() const 
{
  this->CheckValidity();
  ChainViewList::const_iterator it=data_->chains.begin();
  int count=0;
  while (it!=data_->chains.end()) {
    count+=(*it).GetResidueCount();
    ++it;
  }
  return count;
}

geom::Vec3 EntityView::GetCenterOfAtoms() const 
{
  geom::Vec3 center;
  
  if (!this->HasAtoms()) {
    return center;
  }
  unsigned int counter = 0;
  for (ChainViewList::const_iterator ci = this->GetChainList().begin(),
       ce = this->GetChainList().end(); ci != ce; ++ci) {
    for (ResidueViewList::const_iterator ri = ci->GetResidueList().begin(),
         re = ci->GetResidueList().end(); ri != re; ++ri) {
      for (AtomViewList::const_iterator ai = ri->GetAtomList().begin(),
           ae = ri->GetAtomList().end(); ai != ae; ++ai) {
        center += ai->GetPos();
        counter += 1;
      }
    }
  }
  return center/static_cast<Real>(counter);
}

geom::Vec3 EntityView::GetCenterOfMass() const 
{
  geom::Vec3 center;
  Real mass = this->GetMass();
  if (mass==0) {
    return center;
  }

  for (ChainViewList::const_iterator ci = this->GetChainList().begin(),
       ce = this->GetChainList().end(); ci != ce; ++ci) {
    for (ResidueViewList::const_iterator ri = ci->GetResidueList().begin(),
         re = ci->GetResidueList().end(); ri != re; ++ri) {
      for (AtomViewList::const_iterator ai = ri->GetAtomList().begin(),
           ae = ri->GetAtomList().end(); ai != ae; ++ai) {
        center += ai->GetPos()*ai->GetMass();
      }
    }
  }
  return center / mass;
}

Real EntityView::GetMass() const 
{
  this->CheckValidity();
  double mass = 0;
  for (ChainViewList::iterator i=data_->chains.begin();i!=data_->chains.end(); ++i) {
    mass+=i->GetMass();
  }
  return mass;
}

int EntityView::GetChainCount() const 
{
  this->CheckValidity();
  return static_cast<int>(data_->chains.size());
}

ChainView EntityView::AddChain(const ChainHandle& chain_handle,
                               ViewAddFlags flags) 
{
  this->CheckValidity();
  ChainView cv;
  if (flags & ViewAddFlag::CHECK_DUPLICATES &&
      (cv=this->ViewForHandle(chain_handle)))
    return cv;
  cv=ChainView(*this, chain_handle);
  data_->chains.push_back(cv);
  if (flags & ViewAddFlag::INCLUDE_RESIDUES) {
    const impl::ResidueImplList& l=chain_handle.Impl()->GetResidueList();
    for (impl::ResidueImplList::const_iterator i=l.begin(); i!=l.end(); ++i) {
      cv.AddResidue(ResidueHandle(*i), true);
    }
  }
  return cv;
}


ResidueView EntityView::AddResidue(const ResidueHandle& residue_handle,
                                   ViewAddFlags flags) {
  this->CheckValidity();
  ChainView cv=this->AddChain(residue_handle.GetChain(),
                              ViewAddFlag::CHECK_DUPLICATES);
  return cv.AddResidue(residue_handle, flags);
}

ResidueView EntityView::AddResidue(const ResidueView& residue_view,
                                   ViewAddFlags flags)
{
  this->CheckValidity();
  ChainView cv=this->AddChain(residue_view.GetChain(),
                              ViewAddFlag::CHECK_DUPLICATES);
  return cv.AddResidue(residue_view, flags);
}

AtomView EntityView::AddAtom(const AtomHandle& atom_handle,
                             ViewAddFlags flags) {
  ChainView cv=this->AddChain(atom_handle.GetResidue().GetChain(),
                              ViewAddFlag::CHECK_DUPLICATES);
  ResidueView rv=cv.AddResidue(atom_handle.GetResidue(),
                               ViewAddFlag::CHECK_DUPLICATES);
  return rv.AddAtom(atom_handle, flags);
}

ChainView EntityView::FindChain(const String& chain_name) const{
  this->CheckValidity();
  ChainViewList::const_iterator i;
  i=std::find_if(data_->chains.begin(), data_->chains.end(),
                 boost::bind(&ChainView::GetName, boost::placeholders::_1)==chain_name);
  return (i!=data_->chains.end()) ? *i : ChainView();
}


ResidueView EntityView::ViewForHandle(const ResidueHandle& residue) const {
  ChainHandle chain=residue.GetChain();
  ChainView v=this->ViewForHandle(chain);
  return v.IsValid() ? v.ViewForHandle(residue) : ResidueView();
}


AtomView EntityView::ViewForHandle(const AtomHandle& atom) const
{
  this->CheckValidity();
  return data_->ViewForHandle(atom);
}

EntityView EntityView::Select(const String& query_string, 
                              QueryFlags flags) const 
{
  return this->Select(Query(query_string), flags);
}

namespace {

  typedef BondTable<AtomView> BondTableType;
  void update_bond_table(BondTableType& bond_table,
                       const AtomView& ai,
                       const AtomView& view)
  {
    const BondHandleList& bonds = ai.GetBondList();
    BondHandleList::const_iterator it;
    for (it=bonds.begin(); it!=bonds.end(); ++it) {
      bond_table.Update(*it, view);
    }
  }
}

EntityView EntityView::Select(const Query& query, QueryFlags flags) const 
{
  this->CheckValidity();
  // Do the actual selection.
  bool c_added, r_added;
  EntityView view(this->GetHandle());
  BondTableType bond_table;
  
  QueryState qs(query.CreateQueryState(*this));
  //BondHashTable bond_hash_table;
  ChainViewList::const_iterator chain_it=data_->chains.begin();
  for (; chain_it!=data_->chains.end(); ++chain_it) {
      ChainView src_chain=*chain_it;
      c_added=false;
      tribool c=qs.EvalChain(src_chain.Impl());
      if (c==true) {
        // add everything , i.e. all residues and atoms of this chain.
        // We can't use AddChain(chain, ViewAddFlag::INCLUDE_ALL) since we need
        // to keep track of all atoms added to get the bonds right.
        ChainView dst_chain=view.AddChain(src_chain.GetHandle());
        const ResidueViewList& residues=src_chain.GetResidueList();
        ResidueViewList::const_iterator res_it;
        for (res_it=residues.begin() ;res_it!=residues.end(); ++res_it) {
          ResidueView src_res=*res_it;
          ResidueView dst_res=dst_chain.AddResidue(src_res.GetHandle());
          AtomViewList::const_iterator at_it = src_res.GetAtomList().begin();
          for (;at_it != src_res.GetAtomList().end(); ++at_it) {
            if (flags & QueryFlag::NO_BONDS) {
              dst_res.AddAtom(*at_it);
            } else {
              update_bond_table(bond_table, 
                                *at_it, dst_res.AddAtom((*at_it).GetHandle()));              
            }
          }
        }
      } else if (indeterminate(c)) {
        // Test residues
        ChainView dst_chain;
        ResidueViewList::const_iterator res_it=src_chain.GetResidueList().begin();
        for (; res_it!=src_chain.GetResidueList().end(); ++res_it) {
          ResidueView src_res=*res_it;
          tribool r=qs.EvalResidue(src_res.Impl());
          if (r==true) {
            // Include all atoms
            if (!c_added) {
              c_added=true;
              dst_chain=view.AddChain(src_chain.GetHandle());
            }
            ResidueView dst_res=dst_chain.AddResidue(src_res.GetHandle());
            AtomViewList::const_iterator at_it = src_res.GetAtomList().begin();
            for (;at_it != src_res.GetAtomList().end(); ++at_it) {
              if (flags & QueryFlag::NO_BONDS) {
                dst_res.AddAtom(*at_it);
              } else {
                update_bond_table(bond_table, *at_it, 
                                  dst_res.AddAtom((*at_it).GetHandle()));              
              }
            }
          } else if (indeterminate(r)) {
            // Test atoms
            r_added = false;
            ResidueView dst_res;
            AtomViewList::const_iterator at_it = src_res.GetAtomList().begin();
            for (;at_it != src_res.GetAtomList().end(); ++at_it) {
              tribool a = qs.EvalAtom((*at_it).Impl());
              if (indeterminate(a) || a==true) {
                if (!c_added) {
                  c_added=true;
                  dst_chain=view.AddChain(src_chain.GetHandle());
                }
                if (!r_added) {
                  r_added = true;
                  dst_res=dst_chain.AddResidue(src_res.GetHandle());
                }
                if (flags & QueryFlag::MATCH_RESIDUES) {
                  AtomViewList::const_iterator at_it2;
                  
                  for (at_it2=src_res.GetAtomList().begin(); 
                       at_it2!=src_res.GetAtomList().end(); ++at_it2) {
                    if (flags & QueryFlag::NO_BONDS) {
                      dst_res.AddAtom(*at_it2);
                    } else {
                      update_bond_table(bond_table, *at_it2, 
                                        dst_res.AddAtom((*at_it2).GetHandle()));              
                    }                    
                  }
                  break;
                } else {
                  if (flags & QueryFlag::NO_BONDS) {
                    dst_res.AddAtom(*at_it);
                  } else {
                    update_bond_table(bond_table, *at_it, 
                                      dst_res.AddAtom((*at_it).GetHandle()));              
                  }                  
                }

              }

              qs.Reset(Prop::ATOM);
            }
          }
          qs.Reset(Prop::RESIDUE);
        }
      }
  }
  if (!(flags & QueryFlag::NO_BONDS)) {
    typedef BondTableType::MapType::const_iterator ConstIt;

    for (ConstIt it=bond_table.bonds.begin();it!=bond_table.bonds.end(); ++it) {
      if (flags & QueryFlag::EXCLUSIVE_BONDS || it->second.IsComplete()) {
        BondHandle bh=it->second.bond;
        view.AddBond(bh);
      }
    }    
  }  

  return view;
}

ChainView EntityView::ViewForHandle(const ChainHandle& chain) const {
  this->CheckValidity();
  ChainViewList::const_iterator i;
  i=std::find_if(data_->chains.begin(), data_->chains.end(),
                 bind(&ChainView::GetHandle, boost::placeholders::_1)==chain);
  return (i!=data_->chains.end()) ? *i : ChainView();
}


bool EntityView::IsChainIncluded(const ChainHandle& chain) const 
{
  this->CheckValidity();
  return this->ViewForHandle(chain).IsValid();
}


bool EntityView::AddBond(const BondHandle& bh) 
{
  this->CheckValidity();
  AtomView av1=data_->ViewForHandle(bh.GetFirst());
  AtomView av2=data_->ViewForHandle(bh.GetSecond());
  bool add=true;
  if (av1) {
    add=av1.AddBond(bh);
  }
  if (av2) {
    add&=av2.AddBond(bh);
  }
  if (add) {
    data_->bonds.push_back(bh);
    return true;
  }
  return false;
}

const BondHandleList& EntityView::GetBondList() const 
{
  this->CheckValidity();
  return data_->bonds;
}

const ChainViewList& EntityView::GetChainList() const 
{
  this->CheckValidity();
  return data_->chains;
}

void EntityView::RemoveAtom(AtomView view) {
  this->CheckValidity();
  if (!view.IsValid())
    return;
  view.GetResidue().RemoveAtom(view);
}


void EntityView::RemoveResidue(ResidueView view) {
  this->CheckValidity();
  if (!view.IsValid())
    return;
  view.GetChain().RemoveResidue(view);
}


void EntityView::RemoveChain(ChainView view) {
  this->CheckValidity();
  if (!view.IsValid())
    return;
  view.RemoveResidues();
  for (ChainViewList::iterator i=data_->chains.begin(), 
       e=data_->chains.end(); i!=e; ++i) {
    if (*i==view) {
      data_->chains.erase(i);
      break;
    }
  }
}


void EntityView::RemoveBond(const BondHandle& bond) 
{
  this->CheckValidity();
  BondHandleList::iterator i=std::find(data_->bonds.begin(), 
                                        data_->bonds.end(), bond);
                                        
  if (i!=data_->bonds.end()) {
    data_->bonds.erase(i);
    AtomView av1=data_->ViewForHandle(bond.GetFirst());
    AtomView av2=data_->ViewForHandle(bond.GetSecond());
    if (av1) {
      av1.RemoveBondInternal(bond);
    }
    if (av2) {
      av2.RemoveBondInternal(bond);
    }    
  }
}

int EntityView::GetBondCount() const 
{
  this->CheckValidity();
  return static_cast<int>(data_->bonds.size());
}

EntityView EntityView::Copy() const 
{
  this->CheckValidity();
  EntityView ev(this->GetHandle());
  ChainViewList::const_iterator it=data_->chains.begin();
  for (; it!=data_->chains.end(); ++it) {
    ev.AddChain(*it, ViewAddFlag::INCLUDE_ALL);
  }
  for (BondHandleList::const_iterator i=data_->bonds.begin(), 
    e=data_->bonds.end(); i!=e; ++i) {
    ev.AddBond(*i);    
  }  
  return ev;
}

ChainView EntityView::AddChain(const ChainView& chain_view,
                               ViewAddFlags flags) 
{
  this->CheckValidity();
  ChainView cv;
  if (flags & ViewAddFlag::CHECK_DUPLICATES &&
     (cv=this->ViewForHandle(chain_view.GetHandle()))) {
    if (!(flags & ViewAddFlag::INCLUDE_RESIDUES)) {
      return cv;
    }
  } else {
    cv=ChainView(*this, chain_view.GetHandle());
  }
  data_->chains.push_back(cv);
  if (flags & ViewAddFlag::INCLUDE_RESIDUES) {
    ResidueViewList::const_iterator i=chain_view.GetResidueList().begin();
    for (; i!=chain_view.GetResidueList().end(); ++i) {
      cv.AddResidue(*i, flags);
    }
  }
  return cv;
}

ResidueViewList EntityView::GetResidueList() const
{
  this->CheckValidity();
  ResidueViewList residues;
  residues.reserve(this->GetResidueCount());
  ChainViewList::const_iterator i;
  for (i=data_->chains.begin(); i!=data_->chains.end(); ++i) {
    std::copy(i->GetResidueList().begin(), i->GetResidueList().end(), 
              std::back_inserter(residues));
  }
  return residues;
}


AtomViewList EntityView::GetAtomList() const
{
  this->CheckValidity();
  AtomViewList atoms;
  atoms.reserve(this->GetAtomCount());
  ChainViewList::const_iterator i;
  for (i=data_->chains.begin(); i!=data_->chains.end(); ++i) {
    for (ResidueViewList::const_iterator j = (*i).GetResidueList().begin(), 
         e = (*i).GetResidueList().end(); j!=e; ++j) {
      std::copy(j->GetAtomList().begin(), j->GetAtomList().end(), 
                std::back_inserter(atoms));
    }
  }
  return atoms;
}

AtomViewList EntityView::FindWithin(const geom::Vec3& center, Real radius) const
{
  // could use some heuristic here to decide whether it is useful
  // to use the FindWithin method of the handle or do a full search of
  // all atoms included in the view.
  this->CheckValidity();
  AtomHandleList ahl=this->GetHandle().FindWithin(center, radius);
  AtomViewList avl;
  for (AtomHandleList::iterator i=ahl.begin(), e=ahl.end(); i!=e; ++i) {
    if (AtomView v=this->ViewForHandle(*i)) {
      avl.push_back(v);
    }
  }
  return avl;
}

namespace {
  
class InclusiveBondAdder: public EntityVisitor
{
public:
  InclusiveBondAdder(EntityView& view, EntityViewData& data): 
     view_(view), data_(data)
  { }
  virtual bool VisitBond(const BondHandle& bh) {
    if (data_.ViewForHandle(bh.GetFirst()) && 
        data_.ViewForHandle(bh.GetSecond())) {
      view_.AddBond(bh);
    }
    return true;
  }
private:
  EntityView&     view_;
  EntityViewData& data_;
};

}

void EntityView::AddAllInclusiveBonds()
{
  this->CheckValidity();
  for (BondHandleList::iterator i=data_->bonds.begin(),
       e=data_->bonds.end(); i!=e; ++i) {
    BondHandle bond=*i;
    AtomView av1=data_->ViewForHandle(bond.GetFirst());
    AtomView av2=data_->ViewForHandle(bond.GetSecond());
    if (av1) {
      av1.RemoveBondInternal(bond);
    }
    if (av2) {
      av2.RemoveBondInternal(bond);
    }    
  }
  data_->bonds.clear();
  InclusiveBondAdder bond_adder(*this, *data_.get());
  Impl()->Apply(bond_adder);
}

AtomView EntityView::AddAtom(const AtomView& atom_view,
                             ViewAddFlags flags)
{
  return this->AddAtom(atom_view.GetHandle(), flags);
}

std::pair<Real,Real> EntityView::GetMinMax(const String& prop, 
                                           Prop::Level hint) const
{
  EntityPropertyMapper epm(prop, hint);

  Real min_v=std::numeric_limits<Real>::max();
  Real max_v=-std::numeric_limits<Real>::max();  

  for (ChainViewList::const_iterator ci = this->GetChainList().begin(),
       ce = this->GetChainList().end(); ci != ce; ++ci) {
    for (ResidueViewList::const_iterator ri = ci->GetResidueList().begin(),
         re = ci->GetResidueList().end(); ri != re; ++ri) {
      for (AtomViewList::const_iterator ai = ri->GetAtomList().begin(),
           ae = ri->GetAtomList().end(); ai != ae; ++ai) {
        try {
          Real v=epm.Get(*ai);
          max_v=std::max(v, max_v);
          min_v=std::min(v, min_v);      
        } catch(...) {
          // do nothing in case of missing property
          continue;
        }
      }
    }
  }
  return std::make_pair(min_v,max_v);
}

Real EntityView::GetAngle(const AtomHandle& a1, const AtomHandle& a2,
                const AtomHandle& a3) const
{
  return this->GetHandle().GetAngle(a1,a2,a3);
}

Real EntityView::GetAngle(const AtomView& a1, const AtomView& a2,
                const AtomView& a3) const
{
  return this->GetHandle().GetAngle(a1.GetHandle(),a2.GetHandle(),a3.GetHandle());
}


#ifdef _MSC_VER
AtomView EntityView::AddXAtom(const AtomHandle& ah, ViewAddFlags flags)
{
  return this->AddAtom(ah,flags);
}
AtomView EntityView::FindXAtom(const AtomHandle& ah)
{
  return this->ViewForHandle(ah);
}
#endif


AtomView EntityView::FindAtom(const String& chain_name, 
                              const ResNum& num,
                              const String& atom_name) const
{
  this->CheckValidity();
  ChainView ch=this->FindChain(chain_name);
  if (ch.IsValid()) {
    return ch.FindAtom(num, atom_name);
  }
  return AtomView();
}

ResidueView EntityView::FindResidue(const String& chain,
                                    const ResNum& num) const
{
  this->CheckValidity();
  ChainView ch=this->FindChain(chain);
  if (ch.IsValid()) {
    return ch.FindResidue(num);
 }
  return ResidueView();
}

geom::AlignedCuboid EntityView::GetBounds() const
{
  this->CheckValidity();
  geom::Vec3 mmin( std::numeric_limits<Real>::max());
  geom::Vec3 mmax(-std::numeric_limits<Real>::max());

  if (!this->HasAtoms()) {
    return geom::AlignedCuboid(mmin, mmax);
  }

  for (ChainViewList::const_iterator ci = this->GetChainList().begin(),
       ce = this->GetChainList().end(); ci != ce; ++ci) {
    for (ResidueViewList::const_iterator ri = ci->GetResidueList().begin(),
         re = ci->GetResidueList().end(); ri != re; ++ri) {
      for (AtomViewList::const_iterator ai = ri->GetAtomList().begin(),
           ae = ri->GetAtomList().end(); ai != ae; ++ai) {
        mmax=geom::Max(mmax, ai->GetPos());
        mmin=geom::Min(mmin, ai->GetPos());
      }
    }
  }
  return geom::AlignedCuboid(mmin, mmax);
}


EntityView EntityView::CreateFullView() const
{
  return this->Copy();
}

EntityView EntityView::CreateEmptyView() const
{
  return this->GetHandle().CreateEmptyView();
}

void EntityView::AddAtomInternal(const AtomView& av)
{
  data_->handle_to_view.insert(std::make_pair(av.GetHandle().GetHashCode(), av));
}

void EntityView::RemoveAtomInternal(const AtomView& av)
{
  data_->handle_to_view.erase(av.GetHandle().GetHashCode());
}

EntityView EntityView::ExtendViewToResidues() const
{
  this->CheckValidity();
  EntityView view=this->CreateEmptyView();
  ResidueViewList residues=this->GetResidueList();
  ResidueViewList::const_iterator res_it;
  for (res_it=residues.begin(); res_it!=residues.end(); ++res_it) {
    view.AddResidue((*res_it).GetHandle(),
                    mol::ViewAddFlag::INCLUDE_ALL|mol::ViewAddFlag::CHECK_DUPLICATES);
  }
  view.AddAllInclusiveBonds();
  return view;
}

EntityView EntityView::ExtendViewToSurrounding(Real gap) const
{
  this->CheckValidity();
  EntityView view=this->CreateEmptyView();
  AtomViewList atoms=this->GetAtomList();
  AtomViewList::const_iterator atm_it;
  Real max_dist=5+gap;
  for (atm_it=atoms.begin(); atm_it!=atoms.end(); ++atm_it) {
    view.AddAtom((*atm_it),mol::ViewAddFlag::INCLUDE_ALL|mol::ViewAddFlag::CHECK_DUPLICATES);
    AtomHandleList prot_atoms=this->GetHandle().FindWithin((*atm_it).GetPos(),max_dist);
    AtomHandleList::const_iterator protatm_it;
    for (protatm_it=prot_atoms.begin(); protatm_it!=prot_atoms.end(); ++protatm_it){
      Real dist=geom::Distance((*atm_it).GetPos(), (*protatm_it).GetPos());
      if (dist <= (*atm_it).GetRadius() + (*protatm_it).GetRadius() + gap) {
        view.AddAtom((*protatm_it),mol::ViewAddFlag::INCLUDE_ALL|mol::ViewAddFlag::CHECK_DUPLICATES);
      }
    }
  }
  view.AddAllInclusiveBonds();
  return view;
}

namespace {

class EntityViewDumper : public EntityVisitor {
public:
  EntityViewDumper(std::stringstream& stream, EntityView view):
    stream_(stream), first_atom_(false),
    first_res_(true)
  {
    view.Apply(*this);
  }

  virtual bool VisitAtom(const AtomHandle& atom)
  {
    if (!first_atom_) {
      stream_ << " ";
    } else {
      first_atom_=false;
    }
    stream_ << atom.GetName();
    return true;
  }

  virtual bool VisitResidue(const ResidueHandle& res)
  {
    if (first_res_) {
      first_res_=false;
    } else {
      stream_ << "]" << std::endl;
    }
    stream_ << res.GetName() << res.GetNumber() << " [";
    first_atom_=true;
    return true;
  }

  virtual bool VisitChain(const ChainHandle& chain)
  {
    if (first_res_==false) {
      stream_ << "]" << std::endl;
    }
    stream_ << "chain " << chain.GetName() << std::endl;
    first_res_=true;
    return true;
  }
  virtual void OnExit()
  {
    if (!first_res_) {
      stream_ << "]" << std::endl;
    }
  }
  virtual bool VisitBond(const BondHandle& bond)
  {
    return true;
  }
private:
  std::stringstream& stream_;
  bool first_atom_;
  bool first_res_;
};

}

String EntityView::Dump() const
{
  std::stringstream stream;
  EntityViewDumper dumper(stream, *this);
  return stream.str();
}

ChainView EntityView::FindChain(const ChainHandle& chain) const
{
  LOG_WARNING("EntityView::FindChain is deprecated. "
              "Use EntityView::ViewForHandle instead");
  return this->ViewForHandle(chain);
}

ResidueView EntityView::FindResidue(const ResidueHandle& residue) const
{
  LOG_WARNING("EntityView::FindResidue is deprecated. "
              "Use EntityView::ViewForHandle instead");
  return this->ViewForHandle(residue);
}

AtomView EntityView::FindAtom(const AtomHandle& atom) const
{
  LOG_WARNING("EntityView::FindAtom(handle) is deprecated. "
              "Use EntityView::ViewForHandle instead");
  return this->ViewForHandle(atom);
}

bool EntityView::HasAtoms() const 
{
  this->CheckValidity();
  for (ChainViewList::const_iterator it=data_->chains.begin(), 
       e=data_->chains.end();
       it!=e; ++it) {
    if ((*it).HasAtoms()) {
      return true;
    }
  }
  return false;
}

}} // ns

