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
#ifndef OST_CONNECTOR_IMPL_HH
#define OST_CONNECTOR_IMPL_HH

#include <boost/enable_shared_from_this.hpp>

#include <ost/geom/geom.hh>

#include <ost/mol/impl/atom_impl_fw.hh>
#include <ost/mol/impl/connector_impl_fw.hh>
#include <ost/mol/impl/entity_impl_fw.hh>
#include <ost/mol/entity_visitor_fw.hh>

#include <ost/generic_property.hh>

namespace ost { namespace mol { namespace impl {

class ConnectorImpl: public GenericPropContainerImpl,
                     public boost::enable_shared_from_this<ConnectorImpl> {
public:
  ConnectorImpl(const EntityImplPtr& e, const AtomImplPtr& first, 
                const AtomImplPtr& second, Real l, Real t, Real p, unsigned char bond_order);
  ~ConnectorImpl();
  
  AtomImplPtr GetFirst() const {return first_.lock();}
  AtomImplPtr GetSecond() const {return second_.lock();}

  geom::Vec3 GetPos() const;
  geom::Vec3 GetOriginalPos() const;

  Real GetLength() const;
  void SetLength(Real l) {len_=l;}
  
  // Tests whether this connector is connecting atoms a and b.
  bool IsConnectorOf(const AtomImplPtr& a, const AtomImplPtr& b) const;
  


  // reverse order of atoms, does not affect angles!
  void Switch();

  void Apply(EntityVisitor& v);

  geom::Vec3 GetDir() const {
    return local_rot_*geom::Vec3(0, 0, 1);
  }
  
  void SetDirAndLength(const geom::Vec3& dir);  
  
  void SetDir(const geom::Vec3& dir);
  
  const geom::Mat3& GetLocalRot() const { return local_rot_; }  
  
  void SetBondOrder(unsigned char bo) { bond_order_ = bo; }
  
  unsigned char GetBondOrder() const { return bond_order_; }

private:
  AtomImplW first_,second_;
  Real len_;
  unsigned char bond_order_;
  // rotation matrix that brings a vector from the local coordinates of the 
  // first to the local coordinate system of the second atom.
  geom::Mat3 local_rot_;
};

}}} // ns

#endif
