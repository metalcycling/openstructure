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
#include <cmath>

#include <ost/log.hh>
#include <ost/mol/atom_handle.hh>
#include "atom_impl.hh"
#include "residue_impl.hh"
#include "connector_impl.hh"
#include "fragment_impl.hh"
#include "entity_impl.hh"
#include "residue_impl.hh"
#include "chain_impl.hh"
#include "torsion_impl.hh"
#include <ost/mol/entity_visitor.hh>
#include <boost/bind/bind.hpp>
#include <algorithm>


namespace ost {  namespace mol { namespace impl {

AtomImpl::AtomImpl(const EntityImplPtr& e,
                   const ResidueImplPtr& r,
                   const String& n,
                   const geom::Vec3& p,
                   const String& ele,
                   unsigned long index):
  res_(r),
  name_(n),
  pos_(p),
  occupancy_(0.0),
  b_factor_(1.0),
  prop_(NULL),  
  is_hetatm_(false),  
  element_(ele),
  prim_connector_(),
  connector_list_(),
  fragment_(),
  state_(0),
  index_(index)
{
  EntityImplPtr eip = GetEntity();
  tf_pos_ = eip->HasTransform() ? eip->GetTransform().Apply(p) : p;
  prop_=AtomProp::GetDefaultProps(element_);
}

void AtomImpl::AddSecondaryConnector(const ConnectorImplP& bp)
{
  connector_list_.push_back(bp);
}

void AtomImpl::Apply(EntityVisitor& v)
{
  LOG_TRACE("visitor @" << &v << " visiting atom impl @" << this);
  v.VisitAtom(AtomHandle(shared_from_this()));
}

ResidueImplPtr AtomImpl::GetResidue() const
{
  if (!res_.lock()) {
    throw InvalidHandle();
  }
  return res_.lock();
}

void AtomImpl::ClearDirectionality()
{
  if (prim_connector_) {
    connector_list_.push_back(prim_connector_);
  }
  prim_connector_=ConnectorImplP();
}

void AtomImpl::TraceDirectionality(FragmentImplP frag, ConnectorImplP conn, 
                                   int n, unsigned int& count)
{
  if (conn) {
#if !defined(NDEBUG)    
    if (conn->GetFirst()==shared_from_this()) {
      LOG_TRACE("dir:" << String(n,' ') << " atom " << res_.lock()->GetNumber()
                << "." << Name() << "  [" << conn->GetSecond()->GetQualifiedName()
                << " ]");      
    } else {
      LOG_TRACE("dir:" << String(n,' ') << " atom " << res_.lock()->GetNumber()
                << "." << Name() << "  [" << conn->GetFirst()->GetQualifiedName()
                << " ]");
    }
    

#endif              
  } else {
    LOG_TRACE("dir:" << String(n,' ') << " atom " << res_.lock()->GetNumber()
              << "." << Name() << "  [ ]");
  }
  
  // presence of a primary connector indicates ring closure
  if (this->IsTraced()) return;

  // recursive count
  count+=1;
  SetTraced(true);
  fragment_=frag;
  
  AtomImplPtr pthis=shared_from_this();
  if (conn) {
    ConnectorImplList::iterator prim_it=connector_list_.end();
    for (ConnectorImplList::iterator it=connector_list_.begin();
         it!=connector_list_.end();
         ++it) {
      ConnectorImplP c=*it;
      if (c->GetFirst() == pthis) {
        c->GetSecond()->TraceDirectionality(frag, c, n+1, count);
      } else {
        if (c==conn) {
          prim_it=it; // found primary connector
        } else {
          c->Switch(); // reverse order of connector
          c->GetSecond()->TraceDirectionality(frag, c, n+1, count);
        }
      }
    }
    if (prim_it!=connector_list_.end()) {
      prim_connector_=*prim_it;
      connector_list_.erase(prim_it);
    }    
  } else {
    for (ConnectorImplList::iterator it=connector_list_.begin();
         it!=connector_list_.end();
         ++it) {
      ConnectorImplP c=*it;
      if (c->GetFirst() == pthis) {
        c->GetSecond()->TraceDirectionality(frag, c, n+1, count);
      } else {
        c->Switch(); // reverse order of connector
        c->GetSecond()->TraceDirectionality(frag, c, n+1, count);
      }
    }
  }
}


void AtomImpl::UpdateFromICS()
{
  if(this->IsVisited()) return;

  this->SetVisited(true);
  for (ConnectorImplList::iterator i=connector_list_.begin();
       i!=connector_list_.end();++i) {
    if ((*i)->GetSecond().get()==this) continue;
    // since this method is called recursively, with recursion depth of up to
    // 10k, we have to be careful with stack memory.
    // start new scope to tell the compiler he can remove the local vars from the
    // stack before calling UpdateFromICS() on the next atom.
    {
      AtomImplPtr atom=(*i)->GetSecond();
      if (atom.get()==this || atom->IsVisited()==true) {
        continue;
      }
      atom->pos_=pos_+total_rot_*((*i)->GetDir()*(*i)->GetLength());
      atom->total_rot_=total_rot_*(*i)->GetLocalRot();
    }
    (*i)->GetSecond()->UpdateFromICS();
  }
}

void AtomImpl::UpdateFromXCS()
{
  if (this->IsVisited()) return; // ring closure
  this->SetVisited(true);

  // note: local_rot_ and total_rot_ were set in caller.
  for (ConnectorImplList::iterator i=connector_list_.begin();
       i!=connector_list_.end(); ++i) {
    if ((*i)->GetSecond().get()==this)
      continue;
    // since this method is called recursively, with recursion depth of up to
    // 10k, we have to be careful with stack memory.
    // start new scope to tell the compiler he can remove the local vars from the
    // stack before calling UpdateFromICS() on the next atom.
    {
      // Derive direction and length of connector from atom positions.
      geom::Vec3 global_d=((*i)->GetSecond()->OriginalPos()-this->OriginalPos());
      // Set direction and length of connector. Direction is relative to
      // local coordinate system of this atom.
      // Note the order of arguments for the matrix multiplication. This is the
      // same as Tranpose(global_rot_)*global_d <-> Invert(global_rot_)*global_d
      // but avoids the explicit transpose operation.
      (*i)->SetDirAndLength(global_d*total_rot_);
      (*i)->GetSecond()->total_rot_=total_rot_*(*i)->GetLocalRot();
    } // end of new scope
    (*i)->GetSecond()->UpdateFromXCS();
  }
}

AtomImpl::~AtomImpl() {
  if (prop_ && !prop_->is_default) {
    delete prop_;
    prop_=NULL;
  }
}

EntityImplPtr AtomImpl::GetEntity() const
{
  if (!res_.lock()) {
    throw InvalidHandle();
  }
  return res_.lock()->GetEntity();
}

std::ostream& operator<<(std::ostream& o, const AtomImplPtr ap)
{
  o << ap->GetResidue()->GetChain()->GetName() << ".";
  o << ap->GetResidue()->GetKey() << ap->GetResidue()->GetNumber() << ".";
  o << ap->Name();
  return o;
}

bool ConnectorExists(const AtomImplPtr& a, const AtomImplPtr& b) {
  return static_cast<bool>(GetConnector(a, b));
}

ConnectorImplP GetConnector(const AtomImplPtr& a, const AtomImplPtr& b) {
  if (!(a && b))
    return ConnectorImplP();
  ConnectorImplList::const_iterator i;
  AtomImpl* xxx=a.get();

  if (a->GetPrimaryConnector() && a->GetPrimaryConnector()->IsConnectorOf(a, b))
    return a->GetPrimaryConnector();

  if (b->GetPrimaryConnector() && b->GetPrimaryConnector()->IsConnectorOf(a, b))
    return b->GetPrimaryConnector();
  if (a->GetSecondaryConnectors().size()>b->GetSecondaryConnectors().size()) {
    xxx=b.get();
  }
  for (i=xxx->GetSecondaryConnectors().begin();
       i!=xxx->GetSecondaryConnectors().end(); ++i) {
    if ((*i)->GetFirst()==a && (*i)->GetSecond()==b)
      return *i;
    if ((*i)->GetFirst()==b && (*i)->GetSecond()==a)
      return *i;
  }
  return ConnectorImplP();
}

String AtomImpl::GetQualifiedName() const {
  return this->GetResidue()->GetQualifiedName()+"."+this->Name();
}

void AtomImpl::DeleteAllConnectors() {
  EntityImplPtr ent=this->GetEntity();
  if (prim_connector_) {
    ent->DeleteFromConnMap(prim_connector_);
    if (prim_connector_->GetFirst().get()==this) {
      prim_connector_->GetSecond()->DeleteConnector(prim_connector_, false);
    } else {
      prim_connector_->GetFirst()->DeleteConnector(prim_connector_, false);
    }
  }
  prim_connector_=ConnectorImplP();
  ConnectorImplList::iterator i=connector_list_.begin();
  for (; i!=connector_list_.end(); ++i) {
    ent->DeleteFromConnMap(*i);
    if ((*i)->GetFirst().get()==this) {
      (*i)->GetSecond()->DeleteConnector((*i), false);
    } else {
      (*i)->GetFirst()->DeleteConnector((*i), false);
    }
  }
  connector_list_.clear();
}

void AtomImpl::DeleteConnector(const ConnectorImplP& conn,
                               bool delete_other) {
  EntityImplPtr ent=this->GetEntity();                                 
  if (conn==prim_connector_) {
    if (delete_other) {
      ent->DeleteFromConnMap(prim_connector_);
      if (conn->GetFirst().get()==this) {
        conn->GetSecond()->DeleteConnector(conn, false);
      } else {
        conn->GetFirst()->DeleteConnector(conn, false);
      }
    }
    prim_connector_=ConnectorImplP();
    return;
  }
  ConnectorImplList::iterator i;
  i=std::find(connector_list_.begin(), connector_list_.end(), conn);
  if (i!=connector_list_.end()) {
    if (delete_other) {
      ent->DeleteFromConnMap(conn);
      if (conn->GetFirst().get()==this) {
        conn->GetSecond()->DeleteConnector(conn, false);
      } else {
        conn->GetFirst()->DeleteConnector(conn, false);
      }
    }
    connector_list_.erase(i);
  }
}

int AtomImpl::GetIntProperty(Prop::ID prop_id) const
{
  if (prop_id==Prop::AINDEX) {
    return index_;
  }
  throw PropertyError(prop_id);
}

Real AtomImpl::GetFloatProperty(Prop::ID prop_id) const
{
  switch (prop_id) {
    case Prop::AX:
      return pos_[0];
    case Prop::AY:
      return pos_[1];
    case Prop::AZ:
      return pos_[2];            
    case Prop::OCC:
      return occupancy_;
    case Prop::ABFAC:
      return b_factor_;
    case Prop::ACHARGE:
      return prop_->charge;
    default:
      throw PropertyError(prop_id);
  }
}

String AtomImpl::GetStringProperty(Prop::ID prop_id) const
{
  switch (prop_id) {
    case Prop::ANAME:
      return name_;
    case Prop::ELE:
      return element_;
    default:
      throw PropertyError(prop_id);
  }
}

void AtomImpl::DeleteAllTorsions() {
  EntityImplPtr e=this->GetEntity();
  TorsionImplMap::iterator i;
  std::vector<TorsionImplMap::iterator> t_rm_vec;
  std::set<ResidueImpl*> involved_residues;
  for (i=e->GetTorsionMap().begin(); i!=e->GetTorsionMap().end(); ++i) {
     if (i->second->IsAtomInvolved(shared_from_this())) {
        t_rm_vec.push_back(i);
        involved_residues.insert(i->second->GetFirst()->GetResidue().get());
        involved_residues.insert(i->second->GetSecond()->GetResidue().get());
        involved_residues.insert(i->second->GetThird()->GetResidue().get());
        involved_residues.insert(i->second->GetFourth()->GetResidue().get());
     }
  }
  std::vector<TorsionImplMap::iterator>::iterator it_rm;
  for (it_rm=t_rm_vec.begin(); it_rm!=t_rm_vec.end(); ++it_rm) {
     e->GetTorsionMap().erase(*it_rm);
  }
  
  for(std::set<ResidueImpl*>::iterator i = involved_residues.begin(); 
      i != involved_residues.end(); ++i){
    TorsionImplList& l = (*i)->GetTorsionList();
    TorsionImplList::iterator j;
    j=std::remove_if(l.begin(), l.end(),
                     bind(&TorsionImpl::IsAtomInvolved, boost::placeholders::_1, shared_from_this()));
    l.erase(j, l.end());
  }
}

}}} // ns
