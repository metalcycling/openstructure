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
#include <ost/mol/mol.hh>

#include "ring_finder.hh"

namespace ost { namespace conop {
  RingFinder::RingFinder(mol::EntityHandle &eh) :
    ent_(eh), rings_perceived_(false)
  {}
  
  void RingFinder::PerceiveRings()
  {
    mol::AtomHandleList ahl=ent_.GetAtomList();
    for (mol::AtomHandleList::iterator i=ahl.begin();i!=ahl.end();++i) {
      if (i->IsValid()) {
        if (! i->HasProp("RFvisited")) {
          // start ring search from each unvisited atom (e.g. if there are
          // multiple unconnected fragments)
          std::map<long int,bool> ring_closings;
          VisitNext((*i), (*i), ring_closings);
        }
      }
    }
    rings_perceived_=true;
  }

  int RingFinder::GetRingAtomCount()
  {
    if(!ring_view_.IsValid()) {
      GetRings();
    }
    return ring_view_.GetAtomCount();
  }

  int RingFinder::GetRingBondCount()
  {
    if(!ring_view_.IsValid()) {
      GetRings();
    }
    return ring_view_.GetBondCount();
  }

  bool RingFinder::RingsPerceived() {
    return rings_perceived_;
  }

  mol::EntityView RingFinder::GetRings()
  {
    if (!rings_perceived_) {
      PerceiveRings();
    }
    if (!ring_view_.IsValid()) {
      mol::EntityView view=ent_.CreateEmptyView();
      mol::AtomHandleList ahl=ent_.GetAtomList();
      for (mol::AtomHandleList::iterator i=ahl.begin();i!=ahl.end();++i) {
        if (i->HasProp("RFinring")) {
          view.AddAtom(*i);
        }
      }
      mol::BondHandleList bhl=ent_.GetBondList();
      for (mol::BondHandeList::const_iterator j=bhl.begin();j!=bhl.end();++j) {
        if (j->HasProp("RFinring")) {
          view.AddBond(*j);
        }
      }
      ring_view_=view;
    }

    return ring_view_;
  }
  
  int RingFinder::VisitNext(mol::AtomHandle &curr, mol::AtomHandle &prev,
                            std::map<long int,bool> &ring_closings)
  {
    int num_rings_start=ring_closings.size();
    curr.SetBoolProp("RFvisited",true);
    mol::AtomHandleList bp=curr.GetBondPartners();
    for (mol::AtomHandleList::iterator i=bp.begin(); i!=bp.end(); ++i) {
      if ((*i)==prev) {
        // do not visit last visited atom
        continue;
      }
      mol::BondHandle b=curr.FindBondToAtom(*i);
      if (b.HasProp("RFvisited")) {
        // do not go along previously visited bonds
        continue;
      }
      b.SetBoolProp("RFvisited", true);
      if (i->HasProp("RFvisited")) {
        // we have found a ring
        ring_closings[i->GetHashCode()]=true;
        b.SetBoolProp("RFinring", true);
      } else {
        // no ring yet - continue to next atom
        int num_found=VisitNext((*i),curr,ring_closings);
        // back from the recursion
        if (num_found>0) {
          b.SetBoolProp("RFinring", true);
        }
      }
    }
    if (ring_closings.size()-num_rings_start>0) {
      // if current atom is in the ring closing list, remove it
      curr.SetBoolProp("RFinring", true);
      ring_closings.erase(curr.GetHashCode());
    }
    return ring_closings.size()-num_rings_start;
  }

}}
