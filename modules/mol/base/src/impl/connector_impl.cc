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
#include "connector_impl.hh"
#include "atom_impl.hh"
#include "entity_impl.hh"
#include <ost/mol/bond_handle.hh>
#include <ost/mol/entity_visitor.hh>

namespace ost { namespace mol { namespace impl {

ConnectorImpl::ConnectorImpl(const EntityImplPtr& e, const AtomImplPtr& first,
                             const AtomImplPtr& second,
                             Real l, Real t, Real p, unsigned char bond_order):
  first_(first),
  second_(second),
  len_(l),
  bond_order_(bond_order)
{}

void ConnectorImpl::Apply(EntityVisitor& v)
{
  LOG_TRACE("visitor @" << &v << " visiting bond impl @" << this);
  v.VisitBond(BondHandle(shared_from_this()));
}

ConnectorImpl::~ConnectorImpl() {
}


namespace {

// this function finds the rotation matrix that rotates the positive z axis
// onto the given vector. it is assumed that the vector has unit length.
geom::Mat3 find_rotation(const geom::Vec3& d) {
  // assume the vectors are already normalized
  if (std::abs(Real(1.0-Length(d)))>0.00001) {
    LOG_DEBUG("connector find_rotation() has faulty length: " << Length(d));
    assert(0 && "error is big");
  }
  Real dot=d[2];
  if (dot<-0.999999) {
    return geom::Mat3(-1, 0, 0,
                       0, 1, 0,
                       0, 0,-1);
  } else if (dot>0.999999) {
    return geom::Mat3();
  } else {
    geom::Vec3 c=Normalize(geom::Cross(geom::Vec3(0,0,1), d));
    geom::Vec3 v=Normalize(geom::Cross(c, d));
    geom::Mat3 mat(v[0], c[0], d[0],
                   v[1], c[1], d[1],
                   v[2], c[2], d[2]);
    return mat;
  }
}

}

geom::Vec3 ConnectorImpl::GetPos() const 
{
  return (GetFirst()->TransformedPos()+GetSecond()->TransformedPos())*0.5;
}

Real ConnectorImpl::GetLength() const
{
  Real length;
  if (this->GetFirst()->GetEntity()->HasICS()==false) {
    length=geom::Distance(this->GetFirst()->OriginalPos(),this->GetSecond()->OriginalPos());
  } else {
    length=len_;
  }
  return length;
}

void ConnectorImpl::SetDir(const geom::Vec3& dir) {
  if(Length(dir)==0) {
    local_rot_=geom::Mat3(1,0,0,0,1,0,0,0,1);
  } else {
    geom::Vec3 v=Normalize(dir);
    local_rot_=find_rotation(v);
  }
}

void ConnectorImpl::SetDirAndLength(const geom::Vec3& dir) {
  len_=Length(dir);
  if(len_==0) {
    local_rot_=geom::Mat3(1,0,0,0,1,0,0,0,1);
  } else {
    geom::Vec3 v=dir/len_;
    local_rot_=find_rotation(v);
  }
}

bool ConnectorImpl::IsConnectorOf(const AtomImplPtr& a,
                                  const AtomImplPtr& b) const {
  return ((first_.lock()==a && second_.lock()==b) ||
          (first_.lock()==b && second_.lock()==a));
}

geom::Vec3 ConnectorImpl::GetOriginalPos() const 
{
  return (this->GetFirst()->OriginalPos()+
          this->GetSecond()->OriginalPos())*0.5;
}

void ConnectorImpl::Switch()
{
  first_.swap(second_);
}

}}} // ns
