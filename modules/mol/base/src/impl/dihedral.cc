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
#include <ost/mol/impl/atom_impl.hh>
#include <ost/mol/impl/connector_impl.hh>
#include <ost/mol/impl/entity_impl.hh>
#include <ost/mol/not_connected_error.hh>
#include <ost/log.hh>

#include "dihedral.hh"

using namespace geom;

namespace ost { namespace mol { namespace impl {


Dihedral::Dihedral(const AtomImplPtr& atom1,
                   const AtomImplPtr& atom2,
                   const AtomImplPtr& atom3,
                   const AtomImplPtr& atom4) {
  atoms_.push_back(atom1);
  atoms_.push_back(atom2);
  atoms_.push_back(atom3);
  atoms_.push_back(atom4);
}

Dihedral::Dihedral(const AtomImplList& atoms)
  : atoms_(atoms) {

}


Real Dihedral::GetAngleXCS() const {
  return geom::DihedralAngle(atoms_[0]->TransformedPos(), atoms_[1]->TransformedPos(), 
                             atoms_[2]->TransformedPos(), atoms_[3]->TransformedPos());
}

Real Dihedral::GetAngleICS() const {
  const AtomImplPtr& a1=atoms_[0];
  const AtomImplPtr& a2=atoms_[1];
  const AtomImplPtr& a3=atoms_[2];
  const AtomImplPtr& a4=atoms_[3];
  ConnectorImplP conn2=GetConnector(a2, a3);
  geom::Vec3 v1, v2, v3;
  if (conn2->GetFirst()==a2) {
    // to calculate the angle we make use of the same method as in GetAngleXCS
    // by first bringing the direction vectors of the connectors into the local
    // coordinate system of a2.
    if (a2->GetPrimaryConnector() &&
        a2->GetPrimaryConnector()->GetFirst()==a1) {
      v1=geom::Vec3(0, 0, 1);
      v2=conn2->GetDir();
      if (a3->GetPrimaryConnector() &&
          a3->GetPrimaryConnector()->GetFirst()==a4) {
        v3=conn2->GetLocalRot()*geom::Vec3(0, 0, 1);
        throw Error("not implemented");
      } else {
        assert(a3->GetPrimaryConnector()==conn2);
        // both are secondary connectors of a3.
        ConnectorImplP c=GetConnector(a4, a3);
        v3=conn2->GetLocalRot()*c->GetDir();
      }
    } else if (a2->GetPrimaryConnector()==conn2) {
      v1=GetConnector(a1, a2)->GetDir();
      v2=geom::Vec3(0, 0, 1);
      throw Error("not implemented3");
    } else {
      throw Error("not implemented2");
    }
  } else {
    throw Error("not implemented");
  }


  Vec3 v12cross = Cross(v1, v2);
  Vec3 v23cross = Cross(v2, v3);
  return atan2(Dot(v1*Length(v2), v23cross),
               Dot(v12cross, v23cross));
}

namespace {

geom::Vec3 vec_for_angle(const geom::Vec3& v, Real angle) {
  Real circle_radius=sqrt(v[0]*v[0]+v[1]*v[1]);
  geom::Vec3 n(circle_radius*cos(angle), circle_radius*sin(angle), v[2]);
  return n;
}

}
void Dihedral::SetAngleICS(Real angle, bool update_other) {

  const AtomImplPtr& a1=atoms_[0];
  const AtomImplPtr& a2=atoms_[1];
  const AtomImplPtr& a3=atoms_[2];
  const AtomImplPtr& a4=atoms_[3];
  ConnectorImplP conn2=GetConnector(a2, a3);
  if (!conn2) {
    throw NotConnectedError(AtomHandle(a2), AtomHandle(a3));
  }
  // we operate under the assumption that conn2_->First and atom_1 are
  // connected and not conn2_->Second and atom1_.
  LOG_DEBUG("SetAngleICS: " << a1->GetQualifiedName() << " " 
             << a2->GetQualifiedName() << " " << a3->GetQualifiedName() << " " 
             << a4->GetQualifiedName());
  assert(ConnectorExists(a1, conn2->GetFirst()));
  geom::Vec3 v1, v2, v3;
  
  // Make sure conn2 is the primary connector of a3. If you hit another case
  // go ahead and implement it :)
  // not all connectors are living in the same coordinate system. We first
  // have to bring all of them into the local coordinate system of a3.
  if (a2->GetPrimaryConnector() && a2->GetPrimaryConnector()->GetFirst()==a1) {
      v1=geom::Vec3(0, 0,-1)*conn2->GetLocalRot();
      v2=geom::Vec3(0, 0, 1);
      ConnectorImplP c=GetConnector(a4, a3);
      v3=c->GetDir();
      Real phi1=atan2(v1[1], v1[0]);
      geom::Vec3 n3=vec_for_angle(v3, angle+phi1);
      c->SetDir(n3);
      // check if we have to update the other connectors
      const ConnectorImplList& sc=a3->GetSecondaryConnectors();
      if (update_other && sc.size()>1) {
        Real phi2=(angle-phi1)-atan2(v3[1], v3[0]);
        for (ConnectorImplList::const_iterator i=sc.begin(); i!=sc.end(); ++i) {
          if (*i==c)
            continue;
          geom::Vec3 v=(*i)->GetDir();
          Real phi4=phi2+atan2(v[1], v[0]);
          geom::Vec3 n=vec_for_angle(v, phi4);
          (*i)->SetDir(n);
        }
      }
  }
}

AtomImplPtr Dihedral::GetFirst() const
{
  return atoms_[0];
}

AtomImplPtr Dihedral::GetSecond() const
{
  return atoms_[1];
}

AtomImplPtr Dihedral::GetThird() const
{
  return atoms_[2];
}

AtomImplPtr Dihedral::GetFourth() const
{
  return atoms_[3];
}

geom::Vec3 Dihedral::GetPos() const
{
  return GetConnector(atoms_[1], atoms_[2])->GetPos();
}

geom::Vec3 Dihedral::GetOriginalPos() const
{
  return GetConnector(atoms_[1], atoms_[2])->GetOriginalPos();
}

bool Dihedral::IsAtomInvolved(const AtomImplPtr& atom) {
  for (size_t i=0; i <4; ++i) {
    if (atom==atoms_[i])
      return true;
  }
  return false;
}

bool Dihedral::Matches(const AtomImplPtr& a1, const AtomImplPtr& a2,
                       const AtomImplPtr& a3, const AtomImplPtr& a4) const {
  return atoms_[0]==a1 && atoms_[1]==a2 && atoms_[2]==a3 && atoms_[3]==a4;
}


}}} // ns
