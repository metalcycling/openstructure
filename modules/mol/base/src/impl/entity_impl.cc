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
#ifndef _MSC_VER
#include <sys/time.h>
#endif
#include <map>

#include <boost/shared_ptr.hpp>
#include <boost/version.hpp>

#include <boost/version.hpp>
#if BOOST_VERSION>=103900
#include <boost/smart_ptr/make_shared.hpp>
#define MAKE_SHARED_AVAILABLE 1
#else
#define MAKE_SHARED_AVAILABLE 0
#endif
#include <boost/logic/tribool.hpp>
#include <ost/geom/geom.hh>

#include <ost/log.hh>
#include <ost/profile.hh>
#include "atom_impl.hh"
#include "entity_impl.hh"
#include <ost/mol/entity_visitor.hh>
#include <ost/mol/entity_observer.hh>
#include "atom_impl.hh"
#include "residue_impl.hh"
#include "chain_impl.hh"
#include "connector_impl.hh"
#include "torsion_impl.hh"
#include "fragment_impl.hh"
#include "query_impl.hh"
#include <ost/mol/mol.hh>
#include <ost/mol/impl/residue_impl.hh>
#include <ost/mol/bond_table.hh>
#include <ost/mol/atom_handle.hh>
#include <ost/mol/not_connected_error.hh>
#include <ost/mol/query_state.hh>
#include <ost/profile.hh>
#include <ost/integrity_error.hh>
#include "torsion_impl.hh"

using boost::logic::tribool;
using boost::logic::indeterminate;


namespace ost { namespace mol { namespace impl {

EntityImpl::EntityImpl():
  atom_map_(),
  chain_list_(),
  connector_map_(),
  torsion_map_(),
  transform_(),
  has_transform_(false),
  atom_organizer_(5.0),
  fragment_list_(),
  observer_map_(),
  xcs_editor_count_(0),
  ics_editor_count_(0),  
  dirty_flags_(DisableICS),
  name_(""),
  next_index_(0L),
  default_query_flags_(0)
{    
}

int EntityImpl::GetAtomCount() const 
{
  return static_cast<int>(atom_map_.size());
}

int EntityImpl::GetBondCount() const
{
  return static_cast<int>(connector_map_.size());
}

mol::BondHandleList EntityImpl::GetBondList() const
{
  mol::BondHandleList bhl;
  mol::impl::ConnectorImplMap::const_iterator i;
  for (i=connector_map_.begin();i!=connector_map_.end();++i) {
    bhl.push_back(mol::BondHandle((i->second)));
  }
  return bhl;
}

int EntityImpl::GetResidueCount() const 
{
  int count=0;
  for (ChainImplList::const_iterator i=chain_list_.begin(), 
       e=chain_list_.end(); i!=e; ++i) {
    count+=(*i)->GetResidueCount();
  }
  return count;
}


EntityImplPtr EntityImpl::Copy()
{
  Profile prof("EntityImpl::Copy");
#if MAKE_SHARED_AVAILABLE
  EntityImplPtr ent_p=boost::make_shared<EntityImpl>();  
#else
  EntityImplPtr ent_p(new EntityImpl);
#endif

  this->DoCopy(ent_p);
  return ent_p;
}

ChainImplPtr EntityImpl::InsertChain(const String& chain_name, const ChainImplPtr& chain, bool deep)
{
  ChainImplPtr dst_chain=this->InsertChain(chain_name);
  dst_chain->Assign(*chain.get());
  if(deep)
  {
      ResidueImplList::iterator it=chain->GetResidueList().begin(),
                                it_end=chain->GetResidueList().end();
      for(;it!=it_end;++it){
        ResidueImplPtr res=*it;
        dst_chain->AppendResidue(res, true);
      }
  }
  return dst_chain;
}

void EntityImpl::ReplicateHierarchy(EntityImplPtr dest)
{
  for (ChainImplList::const_iterator i=chain_list_.begin(), 
       e1=chain_list_.end(); i!=e1; ++i) {
    ChainImplPtr src_chain=*i;
    dest->InsertChain(src_chain->GetName(),src_chain,true);
  }
}

namespace {


AtomImplPtr lookup_atom(const AtomImplPtr& atom, EntityImplPtr dest, 
                        int atom_index=-1)
{
  if (atom_index==-1) {
    ResidueImplPtr res=atom->GetResidue();
    AtomImplList::const_iterator i=std::find(res->GetAtomList().begin(),
                                             res->GetAtomList().end(),
                                             atom);
    if (i==res->GetAtomList().end())
      return AtomImplPtr();
    atom_index=i-res->GetAtomList().begin();
  }
  int res_index=atom->GetResidue()->GetIndex();
  ChainImplPtr chain=dest->FindChain(atom->GetResidue()->GetChain()->GetName());
  return ((chain->GetResidueList())[res_index]->GetAtomList())[atom_index];
} 
 
}

void EntityImpl::DoCopyBondsAndTorsions(EntityImplPtr dest)
{
  for (ChainImplList::const_iterator i=chain_list_.begin(), 
       i2=dest->chain_list_.begin(), e1=chain_list_.end(); i!=e1; ++i, ++i2) {
    ChainImplPtr src_chain=*i;
    ChainImplPtr dst_chain=*i2;
    for (ResidueImplList::iterator j=src_chain->GetResidueList().begin(),
         j2=dst_chain->GetResidueList().begin(),
         e2=src_chain->GetResidueList().end(); j!=e2; ++j, ++j2) {
      ResidueImplPtr src_res=*j, dst_res=*j2;
      int res_index=j-src_chain->GetResidueList().begin();
      for (AtomImplList::iterator k=src_res->GetAtomList().begin(),
           k2=dst_res->GetAtomList().begin(),
           e3=src_res->GetAtomList().end(); k!=e3; ++k, ++k2) {
        AtomImplPtr src_atom=*k;
        AtomImplPtr dst_atom=*k2;
        const ConnectorImplList& sec_conn=src_atom->GetSecondaryConnectors();
        for (ConnectorImplList::const_iterator l=sec_conn.begin(), 
             e4=sec_conn.end(); l!=e4; ++l) {
          ConnectorImplP conn=*l;
          if (conn->GetFirst()==src_atom) {
            AtomImplPtr dst_atom2=lookup_atom(conn->GetSecond(), dest);
            assert(dst_atom2);
            dest->Connect(dst_atom, dst_atom2, 0, 0, 0, conn->GetBondOrder());
          }
        }
      }
      
      const TorsionImplList& torsions=src_res->GetTorsionList();
      for (TorsionImplList::const_iterator k=torsions.begin(), 
           e3=torsions.end(); k!=e3; ++k) {
        TorsionImplP t=*k;
        // only add torsion if all involved atoms have a residue index 
        // smaller/equal to our residue index
        if (t->GetFourth()->GetResidue()->GetIndex()>res_index ||
            t->GetThird()->GetResidue()->GetIndex()>res_index ||
            t->GetSecond()->GetResidue()->GetIndex()>res_index  ||
            t->GetFirst()->GetResidue()->GetIndex()>res_index) {
          continue;      
        }
        dest->AddTorsion(t->GetName(), lookup_atom(t->GetFirst(), dest),
                         lookup_atom(t->GetSecond(), dest),
                         lookup_atom(t->GetThird(), dest),
                         lookup_atom(t->GetFourth(), dest));
      }
    }
  }
  dest->MarkTraceDirty();
}

void EntityImpl::DoCopy(EntityImplPtr dest)
{
  // first copy the chain - residue - atom hierarchy before replicating the 
  // bond network and the torsions
  dest->SetName(this->GetName());
  dest->next_index_=0;
  this->ReplicateHierarchy(dest);
  this->DoCopyBondsAndTorsions(dest);
}


int EntityImpl::GetChainCount() const 
{
  return static_cast<int>(chain_list_.size());
}

EntityImpl::~EntityImpl()
{
  // notify all observers of pending destruct
  EntityObserverMap::iterator it;
  for (it=observer_map_.begin();it!=observer_map_.end();++it) {
    it->second->OnDestroy();
  }
}
geom::AlignedCuboid EntityImpl::GetBounds() const
{

  if (this->GetAtomCount()>0) {
    geom::Vec3 mmin, mmax;    
    AtomImplMap::const_iterator it=atom_map_.begin();
    mmin=mmax=it->second->TransformedPos();
    for (++it; it!=atom_map_.end();++it) {
      mmin=geom::Min(mmin,it->second->TransformedPos());
      mmax=geom::Max(mmax,it->second->TransformedPos());
    }
    return geom::AlignedCuboid(mmin, mmax);    
  } else {
    return geom::AlignedCuboid(geom::Vec3(), geom::Vec3());
  }

}

geom::Vec3 EntityImpl::GetCenterOfAtoms() const {
  geom::Vec3 center;
  if (this->GetAtomCount()>0) {
    for (AtomImplMap::const_iterator it = atom_map_.begin();
         it!=atom_map_.end();++it) {
      center+=it->second->TransformedPos();
    }
    center/=static_cast<Real>(atom_map_.size());
  }
  return center;
}

geom::Vec3 EntityImpl::GetCenterOfMass() const {
  geom::Vec3 center;
  Real mass = this->GetMass();
  if (this->GetAtomCount()>0 && mass>0) {
    for(AtomImplMap::const_iterator it = atom_map_.begin();it!=atom_map_.end();++it) {
      center+=it->second->TransformedPos()*it->second->GetMass();
    }
    center/=mass;
  }
  return center;
}

  
Real EntityImpl::GetMass() const {
  double mass=0.0;
  for (ChainImplList::const_iterator i=chain_list_.begin(), 
       e=chain_list_.end(); i!=e; ++i) {
    ChainImplPtr chain=*i;
    mass+=chain->GetMass();
  }
  return mass;
}

AtomImplPtr EntityImpl::CreateAtom(const ResidueImplPtr& rp,
                                   const String& name,
                                   const geom::Vec3& pos2,
                                   const String& ele)
{
#if MAKE_SHARED_AVAILABLE
  AtomImplPtr ap=boost::make_shared<AtomImpl>(shared_from_this(), rp, name, 
                                              pos2, ele,next_index_++);
#else
  AtomImplPtr ap(new AtomImpl(shared_from_this(), rp, name, pos2, ele, next_index_++));
#endif
  geom::Vec3 pos = has_transform_ ? transform_.Apply(pos2) : pos2;
  ap->TransformedPos()=pos;
  atom_organizer_.Add(ap,pos);
  atom_map_.insert(AtomImplMap::value_type(ap.get(),ap));
  return ap;
}

void EntityImpl::DeleteAtom(const AtomImplPtr& atom) {
  atom_map_.erase(atom.get());
  atom_organizer_.Remove(atom);
}

ResidueImplPtr EntityImpl::CreateResidue(const ChainImplPtr& cp,
                                       const ResNum& n,
                                       const ResidueKey& k) {
#if MAKE_SHARED_AVAILABLE
  return boost::make_shared<ResidueImpl>(shared_from_this(), cp, n, k);
#else
  return ResidueImplPtr(new ResidueImpl(shared_from_this(), cp, n, k));
#endif  
}

ChainImplPtr EntityImpl::InsertChain(const String& cname)
{
  for (ChainImplList::iterator 
       i=chain_list_.begin(), e=chain_list_.end(); i!=e; ++i) {
      if ((*i)->GetName()==cname) {
        throw IntegrityError("Can't insert chain. A chain with name '"+cname+
                             "' already exists");
      }
  }
#if MAKE_SHARED_AVAILABLE
  ChainImplPtr cp=boost::make_shared<ChainImpl>(shared_from_this(), cname);
#else
  ChainImplPtr cp(new ChainImpl(shared_from_this(), cname));
#endif
  chain_list_.push_back(cp);
  return cp;
}

ConnectorImplP EntityImpl::Connect(const AtomImplPtr& first,
                                   const AtomImplPtr& second, Real len,
                                   Real theta, Real phi,
                                   unsigned char bond_order)
{
  if (!first || !second)
    return ConnectorImplP();
  // check for already existing connection
  if(first->GetPrimaryConnector() &&
     first->GetPrimaryConnector()->GetFirst()==second) {
    LOG_DEBUG("connect: returning existing bond (1)");
    return first->GetPrimaryConnector();
  } else if(second->GetPrimaryConnector() &&
            second->GetPrimaryConnector()->GetFirst()==first) {
    LOG_DEBUG("connect: returning existing bond (2)");
    return second->GetPrimaryConnector();
  } else {
    ConnectorImplList clist = first->GetSecondaryConnectors();
    for(ConnectorImplList::const_iterator it=clist.begin();it!=clist.end();++it) {
      if((*it)->GetSecond()==second) {
        LOG_DEBUG("connect: returning existing bond (3)");
        return *it;
      }
    }
    clist = second->GetSecondaryConnectors();
    for(ConnectorImplList::const_iterator it=clist.begin();it!=clist.end();++it) {
      if((*it)->GetSecond()==first) {
        LOG_DEBUG("connect: returning existing bond (4)");
        return *it;
      }
    }
  }
  // no primary connectors are generated, this is left for the directionality 
  // trace
#if MAKE_SHARED_AVAILABLE
  ConnectorImplP bp=boost::make_shared<ConnectorImpl>(shared_from_this(), first, 
                                                      second, len, theta, phi, 
                                                      bond_order);
  assert(bp->GetFirst());
  assert(bp->GetSecond());
#else
  ConnectorImplP bp(new ConnectorImpl(shared_from_this(), first, second, len, 
                                      theta, phi, bond_order));
#endif
  first->AddSecondaryConnector(bp);
  second->AddSecondaryConnector(bp);

  LOG_DEBUG("adding new connector " << bp);
  connector_map_.insert(ConnectorImplMap::value_type(bp.get(),bp));
  return bp;
}

TorsionImplP EntityImpl::AddTorsion(const String& name, const AtomImplPtr& a1,
                                    const AtomImplPtr& a2, const AtomImplPtr& a3,
                                    const AtomImplPtr& a4)
{
  if (ConnectorExists(a2, a3) && ConnectorExists(a1, a2) &&
      ConnectorExists(a3, a4)) {
#if MAKE_SHARED_AVAILABLE
    TorsionImplP impl=boost::make_shared<TorsionImpl>(name, a1, a2, a3, a4);
#else
    TorsionImplP impl(new TorsionImpl(name, a1, a2, a3, a4));
#endif
    torsion_map_.insert(TorsionImplMap::value_type(impl.get(), impl));
    std::set<ResidueImpl*> involved_residues;
    involved_residues.insert(a1->GetResidue().get());
    involved_residues.insert(a2->GetResidue().get());
    involved_residues.insert(a3->GetResidue().get());
    involved_residues.insert(a4->GetResidue().get());
    std::set<ResidueImpl*>::const_iterator i=involved_residues.begin();
    for (; i!=involved_residues.end(); ++i) {
      (*i)->AddTorsion(impl);
    }
    return impl;
  }
  return TorsionImplP();
}


bool EntityImpl::SetAngle(const AtomImplPtr& a1, const AtomImplPtr& a2,
                          const AtomImplPtr& a3, Real angle) 
{
  ConnectorImplP c12=GetConnector(a1, a2);
  ConnectorImplP c23=GetConnector(a2, a3);
  if (c12==c23)
    return false;
  geom::Vec3 v12, v23;
  if (c12 && c12->GetFirst() && c12->GetSecond() && c23 && c23->GetFirst() &&
      c23->GetSecond()) {
    if (c12==a2->GetPrimaryConnector()) {
      v12=geom::Vec3(0, 0, -1);
      v23=c23->GetDir();
    } else if (c23==a2->GetPrimaryConnector()) {
      v12=c12->GetDir();
      v23=geom::Vec3( 0, 0, -1);
    } else {
      v12=c12->GetDir();
      v23=c23->GetDir();
    }
    geom::Vec3 normal;
    if (std::abs(Dot(v12, v23)) > 0.999) {
      //The two vectors are colinear. This means we have the freedom to choose
      //the normal we want to rotate around, as long as it is orthogonal to the
      // two vectors.
      normal=OrthogonalVector(v12);
    } else {
      normal=geom::Cross(v12, v23);
    }
    geom::Mat3 rot=AxisRotation(normal, angle);

    if (c23==a2->GetPrimaryConnector()) {
      assert(0 && "not implemented yet.");
    }  else {
      c23->SetDir(rot*v12);
    }
    return true;
  }
  return false;
}

Real EntityImpl::GetAngleXCS(const AtomImplPtr& a1, const AtomImplPtr& a2,
                               const AtomImplPtr& a3) const 
{
  ConnectorImplP c12=GetConnector(a1, a2);
  ConnectorImplP c23=GetConnector(a2, a3);
  if (c12 && c12->GetFirst() && c12->GetSecond()) {
    if (c23 && c23->GetFirst() && c23->GetSecond()) {
      return Angle(a2->TransformedPos()-a1->TransformedPos(),
                   a2->TransformedPos()-a3->TransformedPos());
    } else {
      AtomHandle ah2(a2), ah3(a3);
      throw NotConnectedError(ah2, ah3);
    }
  } else {
    AtomHandle ah1(a1), ah3(a3);
    throw NotConnectedError(ah1, ah3);
  }
}

Real EntityImpl::GetAngleICS(const AtomImplPtr& a1, const AtomImplPtr& a2,
                               const AtomImplPtr& a3) const 
{
   ConnectorImplP c12=GetConnector(a1, a2);
   ConnectorImplP c23=GetConnector(a2, a3);
   if (c12 && c12->GetFirst() && c12->GetSecond()) {
     if (c23 && c23->GetFirst() && c23->GetSecond()) {
        // all that remains is to calculate the angle of the direction vectors.
        // in case one of the connector is a2's primary connector we set it's
        // direction vector to (0, 0, -1).
        if (c12==a2->GetPrimaryConnector()) {
          return Angle(geom::Vec3(0,0,-1), c23->GetDir());
         } else if (c23==a2->GetPrimaryConnector()) {
           return Angle(c12->GetDir(), geom::Vec3(0,0,-1));
         } else {
           // both are secondary connectors.
           return Angle(c12->GetDir(),c23->GetDir());
         }
     } else {
       throw NotConnectedError(AtomHandle(a2), AtomHandle(a3));
     }
   } else {
     throw NotConnectedError(AtomHandle(a1), AtomHandle(a2));
   }
}

Real EntityImpl::GetAngle(const AtomImplPtr& a1, const AtomImplPtr& a2,
                            const AtomImplPtr& a3) const 
{
  if (dirty_flags_ & DirtyXCS) {
    return this->GetAngleICS(a1, a2, a3);    
  } else {
    return this->GetAngleXCS(a1, a2, a3);
  }
}

namespace {

unsigned long ai_ptr_diff(const AtomImplPtr& aip1, const AtomImplPtr& aip2)
{
  return std::abs(reinterpret_cast<long>(aip1.get())-reinterpret_cast<long>(aip2.get()));
}

}

/*
  Trace and verify directionality, as encoded by the first and second atom
  in each connector. The hints from the original order as given by Connect
  are taken into account.
*/
void EntityImpl::TraceDirectionality()
{
  fragment_list_.clear();
  Profile profile_trace("trace directionality");
  // some resetting on the atom level
  for(AtomImplMap::iterator it=atom_map_.begin();it!=atom_map_.end();++it) {
    it->second->ClearDirectionality();
    it->second->SetVisited(false);
    it->second->SetTraced(false);
  }

  /*
    the beginning of a fragment is given by two conditions,
    namely that it does not point to a previous atom, and
    that none of the connectors have it as second.
  */

  // mark all second atoms
  for(ConnectorImplMap::iterator it=connector_map_.begin();
      it!=connector_map_.end();
      ++it) {
    it->second->GetSecond()->SetVisited(true);
  }

  // search atoms that fulfill above criteria, and start recursive trace with them
  unsigned int traced_atom_count=0;
  unsigned int traced_atom_count_prev=atom_map_.size()+1;
  while (true) {
    for(AtomImplMap::iterator it=atom_map_.begin();it!=atom_map_.end();++it) {
      if (!it->second->IsVisited() && !it->second->HasPrevious()) {
        it->second->SetVisited(true);
        LOG_DEBUG("dir: new fragment");
#if MAKE_SHARED_AVAILABLE
        FragmentImplP frag=boost::make_shared<FragmentImpl>(it->second);
#else
        FragmentImplP frag( new FragmentImpl(it->second));
#endif
        fragment_list_.push_back(frag);
        it->second->TraceDirectionality(frag,ConnectorImplP(), 0, 
                                        traced_atom_count);
      }
    }
    if(traced_atom_count<atom_map_.size()) {
      LOG_DEBUG("entering closed loops search, since only " << traced_atom_count 
                 << " from " << atom_map_.size() << " atoms were traced");
      /*
  identify closed loops that prohibit
  the trace from completing
      */
      AtomImplPtr aip1;
      AtomImplPtr aip2;
      // step 1: find the connector with the largest atom order gap
      // (this is a heuristic guess for the best atom to isolate)
      for (ConnectorImplMap::iterator it=connector_map_.begin();
           it!=connector_map_.end();
           ++it) {
        if(!it->second->GetFirst()->HasPrevious() &&
           !it->second->GetSecond()->HasPrevious() ) {
          // candidates for loop removal
          if(aip1) {
            if(ai_ptr_diff(it->second->GetFirst(),it->second->GetSecond())>ai_ptr_diff(aip1,aip2)) {
              aip1=it->second->GetFirst();
              aip2=it->second->GetSecond();
            }
          } else {
            aip1=it->second->GetFirst();
            aip2=it->second->GetSecond();
          }
        }
      }
      LOG_DEBUG("found pair [" << aip1 << "] [" << aip2 << "]");
      // step 2: swap all connectors pointing to aip2
      // thus removing its 'visited' state
      for(ConnectorImplMap::iterator it=connector_map_.begin();
         it!=connector_map_.end();
         ++it) {
        if(it->second->GetSecond()==aip2) {
          it->second->Switch();
          it->second->GetSecond()->SetVisited(true);
        }
      }
      aip2->SetVisited(false);
      // and now back to the recursive trace
    } else {
      // finished trace
      break;
    }

    // infinite loop safeguard
    if(traced_atom_count == traced_atom_count_prev) {
      LOG_VERBOSE("Encountered unbreakable locked state during directionality trace");
      break;
    }
    traced_atom_count_prev=traced_atom_count;
  }
  LOG_VERBOSE("Directionality trace completed with " << fragment_list_.size()
               << " fragment(s)");
}

bool EntityImpl::IsXCSDirty() const
{
  return (dirty_flags_ & DirtyXCS)!=0;
}

void EntityImpl::UpdateFromICS()
{
  LOG_DEBUG("updating external from internal coordinate system");
  Profile prof_conv("update external from internal coordinate system");    
  // traverse through all atoms
  for(AtomImplMap::iterator it=atom_map_.begin();it!=atom_map_.end();++it) {
    it->second->SetVisited(false);
  }

  for (FragmentImplList::iterator it=fragment_list_.begin();
       it!=fragment_list_.end();++it) {
    LOG_DEBUG("entering root of tree for ICS update");
    assert(*it);
    AtomImplPtr ap=(*it)->GetAtom();
    // recursive update, using position of first atom
    /*geom::Mat4 pmat=(*it)->GetTMat();
    pmat=tmat*pmat;
    ap->SetICSMat(pmat);*/
    ap->UpdateFromICS();
    LOG_DEBUG("done with tree traversal for ICS update");
  }

  LOG_DEBUG("refreshing all connector positions");
  // refresh all connector positions
  /*
  for (ConnectorImplMap::iterator it=connector_map_.begin();
       it!=connector_map_.end();++it) {
    it->second->UpdatePosition();
  }*/
  this->UpdateTransformedPos();
}


void EntityImpl::UpdateFromXCS()
{
  LOG_DEBUG("rebuilding internal from external coordinate system");
  Profile prof_conv("updating internal from external coordinate system");
  if(fragment_list_.empty()) {
    // this should have been handled by conopology after the
    // initial connectivity assignments - if not, then do it here
    Profile profile_trace("finished directionality trace");    
    LOG_DEBUG("warning: no fragment list, re-running directionality trace");
    TraceDirectionality();
  }

  if(Logger::Instance().GetVerbosityLevel()>Logger::DEBUG) {
    LOG_TRACE("dumping directionality");
    for(AtomImplMap::iterator it=atom_map_.begin();it!=atom_map_.end();++it) {
      LOG_TRACE(" " << it->second << ":");
      ConnectorImplList clist = it->second->GetSecondaryConnectors();
      for(ConnectorImplList::const_iterator cit=clist.begin();cit!=clist.end();++cit) {
        LOG_TRACE("  " << (*cit)->GetFirst() << " -> " << (*cit)->GetSecond());
      }
    }
  }


  for(AtomImplMap::iterator it=atom_map_.begin();it!=atom_map_.end();++it) {
    it->second->SetVisited(false);
  }

  for (FragmentImplList::iterator it=fragment_list_.begin();
       it!=fragment_list_.end();++it) {
    LOG_DEBUG("entering fragment tree traversal for XCS update");
    FragmentImplP frag= *it;
    AtomImplPtr atom = frag->GetAtom();
    {
      atom->UpdateFromXCS();
    }
  }
}

void EntityImpl::Apply(EntityVisitor& v)
{
  LOG_TRACE("visitor @" << &v << " visiting entity impl @" << this);
  v.OnEntry();
  for(ChainImplList::iterator 
      it = chain_list_.begin();it!=chain_list_.end();++it) {
    (*it)->Apply(v);
  }

  for(ConnectorImplMap::iterator it = connector_map_.begin();it!=connector_map_.end();++it) {
    it->second->Apply(v);
  }

  for(TorsionImplMap::iterator it = torsion_map_.begin();it!=torsion_map_.end();++it) {
    it->second->Apply(v);
  }
  v.OnExit();
}

void EntityImpl::ApplyTransform(const geom::Transform& tf)
{
  SetTransform(transform_.Apply(tf));
}

void EntityImpl::SetTransform(const geom::Transform& tf)
{
  transform_=tf;
  has_transform_=true;
  this->UpdateTransformedPos();
  this->MarkOrganizerDirty();
}

void EntityImpl::ClearTransform()
{
  SetTransform(geom::Transform());
  has_transform_=false;
}

void EntityImpl::FixTransform()
{
  if(!has_transform_) return;
  for(AtomImplMap::iterator it = atom_map_.begin();it!=atom_map_.end();++it) {
    it->second->OriginalPos()=it->second->TransformedPos();
  }
  transform_=geom::Transform();
  has_transform_=false;
  this->UpdateTransformedPos();
  this->MarkOrganizerDirty();
}

void EntityImpl::AttachObserver(const EntityObserverPtr& o)
{
  observer_map_.insert(EntityObserverMap::value_type(o.get(),o));
}

void EntityImpl::DetachObserver(const EntityObserverPtr& o)
{
  EntityObserverMap::iterator it=observer_map_.find(o.get());

  if(it != observer_map_.end()) {
    observer_map_.erase(it);
  }
}

void EntityImpl::NotifyObserver()
{
  for(EntityObserverMap::iterator it=observer_map_.begin();it!=observer_map_.end();++it) {
    it->second->OnUpdate();
  }
}

void EntityImpl::Swap(EntityImpl& impl)
{
  atom_map_.swap(impl.atom_map_);
  chain_list_.swap(impl.chain_list_);
  connector_map_.swap(impl.connector_map_);
  torsion_map_.swap(impl.torsion_map_);
  std::swap(transform_,impl.transform_);
  std::swap(has_transform_,impl.has_transform_);
  atom_organizer_.Swap(impl.atom_organizer_);
  fragment_list_.swap(impl.fragment_list_);
  observer_map_.swap(impl.observer_map_);
}


TorsionImplP EntityImpl::FindTorsion(const AtomImplPtr& a1,
                                     const AtomImplPtr& a2,
                                     const AtomImplPtr& a3,
                                     const AtomImplPtr& a4) const {
  TorsionImplMap::const_iterator i;
  
  for (i=torsion_map_.begin(); i!=torsion_map_.end(); ++i) {
    TorsionImplP t = i->second;
    if (t->Matches(a1, a2, a3, a4)) {
      return t;
    }
  }
  return TorsionImplP();
}

AtomImplList EntityImpl::FindWithin(const geom::Vec3& pos, Real radius)
{
  this->UpdateOrganizerIfNeeded();
  SpatialAtomOrganizer::ItemList alist = atom_organizer_.FindWithin(pos,radius);
  return alist;
}

namespace {


  typedef BondTable<AtomView> BondTableType;
  void update_bond_table(BondTableType& bond_table,
                       const AtomImplPtr& ai,
                       const AtomView& view)
  {
    bond_table.Update(BondHandle(ai->GetPrimaryConnector()), view);
    const ConnectorImplList& clist = ai->GetSecondaryConnectors();
    ConnectorImplList::const_iterator it;
    for (it=clist.begin(); it!=clist.end(); ++it) {
      bond_table.Update(BondHandle(*it), view);
    }
  }
}

template<bool always_true>
EntityView EntityImpl::do_selection(const EntityHandle& eh,
                                    const Query& query,
                                    QueryFlags flags) const
{

  BondTableType bond_table;
  LOG_DEBUG("entering do selection");

  size_t chain_count = 0, residue_count = 0, atom_count = 0;
  bool c_added, r_added;
  LOG_DEBUG("creating view");
  EntityView view(eh);
  LOG_DEBUG("creating query state");
  EntityHandle myself(const_cast<impl::EntityImpl*>(this)->shared_from_this());
  QueryState qs(query.CreateQueryState(myself));
  LOG_DEBUG("entering chain loop");
  for (ChainImplList::const_iterator 
       ch_it=chain_list_.begin(); ch_it!=chain_list_.end();++ch_it) {
    LOG_DEBUG("checking chain " << (*ch_it)->GetName());
    c_added = false;
    tribool c = always_true ? tribool(true) : qs.EvalChain(*ch_it);
    if (c == true) {
      LOG_DEBUG("chain is selected");
      // Include all residues
      const ChainImplPtr& ci=*ch_it;
      ++chain_count;
      ChainView chain=view.AddChain(ci);
      ResidueImplList::const_iterator re_it = ci->GetResidueList().begin();
      for ( ;re_it!=ci->GetResidueList().end();++re_it) {
        LOG_DEBUG(" adding residue " << (*re_it)->GetNumber());
        ResidueImplPtr& ri = const_cast<ResidueImplPtr&>(*re_it);
        ++residue_count;
        ResidueView res =chain.AddResidue(ri);
        AtomImplList::const_iterator at_it = ri->GetAtomList().begin();
        for (;at_it != ri->GetAtomList().end(); ++at_it) {
          LOG_DEBUG("  adding atom " << (*at_it)->GetName());
          ++atom_count;
          if (flags & QueryFlag::NO_BONDS) {
            res.AddAtom(*at_it);
          } else {
            update_bond_table(bond_table, *at_it, res.AddAtom(*at_it));              
          }
        }
      }
    } else if (indeterminate(c)) {
      // Test residues
      r_added = false;
      const ChainImplPtr& ci = *ch_it;

      ChainView chain;
      ResidueImplList::const_iterator re_it = ci->GetResidueList().begin();
      for ( ;re_it!=ci->GetResidueList().end();++re_it) {
        LOG_DEBUG(" checking residue " << (*re_it)->GetNumber());
        tribool r = qs.EvalResidue(*re_it);
        ResidueImplPtr& ri = const_cast<ResidueImplPtr&>(*re_it);
        if (r == true) {
          LOG_DEBUG(" residue is selected");
          // Include all atoms
          ResidueImplPtr& ri = const_cast<ResidueImplPtr&>(*re_it);
          AtomImplList::const_iterator at_it = ri->GetAtomList().begin();
          if (!c_added) {
            c_added=true;
            ++chain_count;
            chain=view.AddChain(ci);
          }
          ++residue_count;
          ResidueView res=chain.AddResidue(ri);
          for (;at_it != ri->GetAtomList().end(); ++at_it) {
            LOG_DEBUG("  adding atom " << (*at_it)->GetName());
            ++atom_count;
            if (flags & QueryFlag::NO_BONDS) {
              res.AddAtom(*at_it);
            } else {
              update_bond_table(bond_table, *at_it, res.AddAtom(*at_it));              
            }
          }
        } else if (indeterminate(r)) {
          // Test atoms
          r_added = false;

          ResidueView res;

          AtomImplList::const_iterator at_it = ri->GetAtomList().begin();
          for (;at_it != ri->GetAtomList().end(); ++at_it) {
            tribool a = qs.EvalAtom(*at_it);
            if (indeterminate(a) || a==true) {
              if (!c_added) {
                c_added = true;
                ++chain_count;
                chain = view.AddChain(ci);
              }
              if (!r_added) {
                r_added = true;
                ++residue_count;
                res = chain.AddResidue(ri);
              }
              if (flags & QueryFlag::MATCH_RESIDUES) {
                AtomImplList::const_iterator at_it2 = ri->GetAtomList().begin();
                for (;at_it2 != ri->GetAtomList().end(); ++at_it2) {
                  LOG_DEBUG("  adding atom " << (*at_it2)->GetQualifiedName());
                  ++atom_count;                  
                  if (flags & QueryFlag::NO_BONDS) {
                    res.AddAtom(*at_it2);
                  } else {
                    update_bond_table(bond_table, *at_it2, 
                                      res.AddAtom(*at_it2));
                  }                  
                }          
                break;
              } else {
                LOG_DEBUG("  adding atom " << (*at_it)->GetName());
                ++atom_count;                
                if (flags & QueryFlag::NO_BONDS) {
                  res.AddAtom(*at_it);
                } else {
                  update_bond_table(bond_table, *at_it, res.AddAtom(*at_it));
                }
              }
            }
            qs.Reset(Prop::ATOM);
          }
        }
        qs.Reset(Prop::RESIDUE);
      }
    }
    if(!always_true) {
      qs.Reset(Prop::CHAIN);
    }
  }
  int bond_count=0;
  if (!(flags & QueryFlag::NO_BONDS)) {
    LOG_DEBUG("adding bonds");
    typedef BondTableType::MapType::const_iterator ConstIt;
    for (ConstIt it=bond_table.bonds.begin();it!=bond_table.bonds.end(); ++it) {
      if ((QueryFlag::EXCLUSIVE_BONDS & flags) || it->second.IsComplete()) {
        BondHandle bh=it->second.bond;
        view.AddBond(bh);
        bond_count++;
      }
    }    
  }
  LOG_VERBOSE("selected " << chain_count << " chain(s), "
               << residue_count << " residue(s), " << atom_count
               << " atom(s)" << " " << bond_count << " bond(s)");
  return view;
}

EntityView EntityImpl::Select(const EntityHandle& h, const Query& q) const
{
  return do_selection<false>(h, q, default_query_flags_);
}

EntityView EntityImpl::Select(const EntityHandle& h, const Query& q, 
                              QueryFlags flags) const
{
  return do_selection<false>(h, q, flags);
}

EntityView EntityImpl::CreateFullView(const EntityHandle& h) const
{
  return do_selection<true>(h, Query(), 0);
}

ChainImplPtr EntityImpl::FindChain(const String& name) const {
  ChainImplList::const_iterator i;
  for(i=chain_list_.begin(); i!=chain_list_.end();++i) {
    if ((*i)->GetName()==name)
      return *i;
  }
  return ChainImplPtr();
}

ResidueImplPtr EntityImpl::FindResidue(const String& chain_name,
                                     const ResNum& residue) const {
  ChainImplPtr chain=this->FindChain(chain_name);
  if (!chain)
    return ResidueImplPtr();
  return chain->FindResidue(residue);
}

AtomImplPtr EntityImpl::FindAtom(const String& chain_name,
                                 const ResNum& residue,
                                 const String& atom_name) const {
  ResidueImplPtr res=this->FindResidue(chain_name, residue);
  if (!res)
    return AtomImplPtr();
  return res->FindAtom(atom_name);
}

void EntityImpl::DeleteFromConnMap(const ConnectorImplP& conn) {
  connector_map_.erase(conn.get());
}

TorsionImplMap& EntityImpl::GetTorsionMap() {
  return torsion_map_;
}

void EntityImpl::DeleteChain(const ChainImplPtr& chain) {
  if (chain && chain->GetEntity().get()==this) {
    chain->DeleteAllResidues();
    chain_list_.erase(this->GetChain(chain->GetName()));
  }
}


void EntityImpl::UpdateOrganizer()
{
  atom_organizer_.Clear();
  for (AtomImplMap::const_iterator i=atom_map_.begin(), 
       e=atom_map_.end(); i!=e; ++i) {
    atom_organizer_.Add(i->second, i->second->TransformedPos());
  }
  dirty_flags_&=~DirtyOrganizer;
}

void EntityImpl::IncXCSEditorCount() 
{
  xcs_editor_count_++;
}

void EntityImpl::UpdateXCSIfNeeded()
{
  if (dirty_flags_ & DirtyTrace) {
    this->TraceDirectionality();    
    dirty_flags_&=~DirtyTrace;
  }
  if (dirty_flags_ & DirtyXCS) {
    LOG_VERBOSE("XCS marked dirty. Updating");
    this->UpdateFromICS();
    dirty_flags_&=~DirtyXCS;
  }
}

void EntityImpl::MarkXCSDirty() 
{
  dirty_flags_|=DirtyXCS|DirtyOrganizer;
}

void EntityImpl::MarkICSDirty()
{
  dirty_flags_|=DirtyICS;  
}

void EntityImpl::MarkTraceDirty()
{
    dirty_flags_|=DirtyTrace;
}

void EntityImpl::MarkOrganizerDirty()
{
    dirty_flags_|=DirtyOrganizer;
}

void EntityImpl::EnableICS()
{
  dirty_flags_=(dirty_flags_ & ~DisableICS) | DirtyTrace;
  this->UpdateICSIfNeeded();
}

void EntityImpl::UpdateICSIfNeeded()
{
  if (dirty_flags_ & DisableICS) {
    return;    
  }
  if (dirty_flags_ & DirtyTrace) {
    LOG_VERBOSE("rerunning directionality trace");
    this->TraceDirectionality();
    this->UpdateFromXCS();
    dirty_flags_&=~DirtyTrace;
  } else if (dirty_flags_ & DirtyICS) {
    LOG_VERBOSE("updating internal from external coordinates");    
    this->UpdateFromXCS();    
    dirty_flags_&=~DirtyICS;    
  }
}

void EntityImpl::DecXCSEditorCount()
{
  assert(xcs_editor_count_>0);  
  xcs_editor_count_--;
  if (xcs_editor_count_==0) {    
    this->UpdateICSIfNeeded();
    this->UpdateOrganizerIfNeeded();    
  }
}

void EntityImpl::UpdateOrganizerIfNeeded()
{
  if (dirty_flags_ & DirtyOrganizer) {
    LOG_DEBUG("atom organizer marked as dirty. updating");
    this->UpdateOrganizer();
    dirty_flags_&=~DirtyOrganizer;
  }
}

void EntityImpl::IncICSEditorCount()
{
  ics_editor_count_++;
}

bool EntityImpl::HasICS() const
{
  return !(dirty_flags_ & DisableICS);
}
void EntityImpl::DecICSEditorCount()
{
  assert(ics_editor_count_>0);
  ics_editor_count_--;
  if (ics_editor_count_==0) {
    this->UpdateXCSIfNeeded();
    this->UpdateOrganizerIfNeeded();
  }
}

impl::ChainImplList::iterator EntityImpl::GetChain(const String& name)
{
  impl::ChainImplList& cc=this->GetChainList();
  for (impl::ChainImplList::iterator i=cc.begin(), e=cc.end(); i!=e; ++i) {
    if ((*i)->GetName()==name) {
      return i;
    }
  }
  return  cc.end();
}

pointer_it<ChainImplPtr> EntityImpl::GetChainIter(const String& name)
{
  impl::ChainImplList& cc=this->GetChainList();
  for (size_t i=0; i<cc.size(); ++i) {
    if (cc[i]->GetName()==name) {
      return &cc.front()+i;
    }
  }
  return  pointer_it<ChainImplPtr>(NULL);
}

void EntityImpl::RenameChain(ChainImplPtr chain, const String& new_name)
{
  //ChainImplList::iterator i;
  ChainImplPtr  ch=this->FindChain(new_name);
  if ((ch) && (ch != chain)) {
    throw IntegrityError("unable to rename chain '"+chain->GetName()+
                         "' to '"+new_name+"', since there is already a chain "
                         "with that name");
  }
  chain->SetName(new_name);
}

void EntityImpl::UpdateTransformedPos(){
  for(AtomImplMap::iterator it = atom_map_.begin();it!=atom_map_.end();++it) {
    it->second->TransformedPos()=has_transform_ ? transform_.Apply(it->second->OriginalPos()) : it->second->OriginalPos();
  }
  for(ChainImplList::iterator cit=chain_list_.begin();cit!=chain_list_.end();++cit) {
    (*cit)->UpdateTransformedPos();
  }
}

const String& EntityImpl::GetName() const {
  return name_;
}

void EntityImpl::SetName(const String& ent_name){
  name_=ent_name;
}

void EntityImpl::ReorderAllResidues()
{
  for(ChainImplList::iterator cit=chain_list_.begin();cit!=chain_list_.end();++cit) {
    (*cit)->ReorderResidues();
  }
}

void EntityImpl::RenumberAllResidues(int start, bool keep_spacing)
{
  for(ChainImplList::iterator cit=chain_list_.begin();cit!=chain_list_.end();++cit) {
    (*cit)->RenumberAllResidues(start, keep_spacing);
  }

}

void EntityImpl::RenumberChain(const String& name, int start, bool keep_spacing)
{
  ChainImplPtr ch=this->FindChain(name);

  if(!ch) {
    throw Error("Could not find chain with name "+name);
  }

  ch->RenumberAllResidues(start, keep_spacing);

}

}}} // ns
