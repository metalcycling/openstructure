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
#include <ost/mol/entity_view.hh>
#include <ost/mol/atom_view.hh>
#include <ost/mol/residue_view.hh>
#include <ost/mol/chain_view.hh>
#include "clash_score.hh"

namespace ost { namespace mol { namespace alg {

Real StericEnergy(const geom::Vec3& pos1, Real r1,
                  const geom::Vec3& pos2, Real r2)
{
  geom::Vec3 d=pos1-pos2;
  // delay calculation of square root as long as possible. It is only needed
  // in the rare case when 0.8254*rr <= d <= rr
  Real distance_sqr=geom::Length2(d); 
  Real rr=r1+r2;
  if (distance_sqr>rr*rr) {
    return 0.0;
  }
  if (distance_sqr<0.8254*0.8254*rr*rr) {
    return 10.0;
  }
  return 57.273*(1.0-sqrt(distance_sqr)/rr);
}

Real ClashScore(const EntityView& ent_a, const EntityView& ent_b)
{
  Real energy=0.0;
  for (ChainViewList::const_iterator ci = ent_a.GetChainList().begin(),
       ce = ent_a.GetChainList().end(); ci != ce; ++ci) {
    for (ResidueViewList::const_iterator ri = ci->GetResidueList().begin(),
         re = ci->GetResidueList().end(); ri != re; ++ri) {
      for (AtomViewList::const_iterator ai = ri->GetAtomList().begin(),
           ae = ri->GetAtomList().end(); ai!= ae; ++ai) {
        AtomViewList clashees=ent_b.FindWithin(ai->GetPos(), 
                                               ai->GetRadius()+1.7);

        for (AtomViewList::iterator j=clashees.begin(), 
            e2=clashees.end(); j!=e2; ++j) {
            energy+=StericEnergy((*j).GetPos(), (*j).GetRadius()-0.25,
                                  ai->GetPos(), ai->GetRadius()-0.25);           
        }
      }
    }
  }
  return energy;  
}

Real ClashScore(const EntityHandle& ent_a, const EntityView& ent_b)
{
  Real energy=0.0;
  for (ChainViewList::const_iterator ci = ent_b.GetChainList().begin(),
       ce = ent_b.GetChainList().end(); ci != ce; ++ci) {
    for (ResidueViewList::const_iterator ri = ci->GetResidueList().begin(),
         re = ci->GetResidueList().end(); ri != re; ++ri) {
      for (AtomViewList::const_iterator ai = ri->GetAtomList().begin(),
           ae = ri->GetAtomList().end(); ai!= ae; ++ai) {
        AtomHandleList clashees=ent_a.FindWithin(ai->GetPos(), 
                                                 ai->GetRadius()+1.7);

        for (AtomHandleList::iterator j=clashees.begin(), 
            e2=clashees.end(); j!=e2; ++j) {
            energy+=StericEnergy((*j).GetPos(), (*j).GetRadius()-0.25,
                                  ai->GetPos(), ai->GetRadius()-0.25);           
        }
      }
    }
  }
  return energy;  
}

Real ClashScore(const AtomHandle& atom, const EntityView& ent_b)
{
  Real energy=0.0;
  AtomViewList clashees=ent_b.FindWithin(atom.GetPos(), 
                                         atom.GetRadius()+2.0);
  for (AtomViewList::iterator j=clashees.begin(), 
       e2=clashees.end(); j!=e2; ++j) {
       energy+=StericEnergy((*j).GetPos(), (*j).GetRadius(),
                            atom.GetPos(), atom.GetRadius());           
  }
  return energy;                                 
}

}}}
