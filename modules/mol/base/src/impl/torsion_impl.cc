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
#include <ost/log.hh>

#include <ost/mol/impl/torsion_impl.hh>
#include <ost/mol/impl/connector_impl.hh>
#include <ost/mol/impl/atom_impl.hh>
#include <ost/mol/entity_visitor.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/mol/torsion_handle.hh>
#include <ost/mol/impl/entity_impl.hh>
namespace ost { namespace mol { namespace impl {

using namespace geom;

TorsionImpl::TorsionImpl(const String& name, const AtomImplPtr& a1,
                         const AtomImplPtr& a2,
                         const AtomImplPtr& a3, 
                         const AtomImplPtr& a4)
  : Dihedral(a1, a2, a3, a4), name_(name) {
}

const String& TorsionImpl::GetName() const {
  return name_;
}

Real TorsionImpl::GetAngle() const {
  if (this->GetFirst()->GetEntity()->IsXCSDirty()) {
    return this->GetAngleICS();
  } else {
    return this->GetAngleXCS();    
  }
}

void TorsionImpl::SetAngle(Real angle, bool update_others) {
  this->SetAngleICS(angle, update_others);
}

void TorsionImpl::Apply(EntityVisitor& v) {
  LOG_TRACE("visitor @" << &v << " visiting torsion impl @" << this);
  v.VisitTorsion(TorsionHandle(shared_from_this()));
}

void TorsionImpl::RotateAngle(Real angle, bool update_others) {
  this->SetAngle(this->GetAngleICS()+angle, update_others);
}

}}} // ns
