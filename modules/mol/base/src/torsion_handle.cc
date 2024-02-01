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
#include <ost/mol/torsion_handle.hh>
#include <ost/mol/entity_visitor.hh>
#include <ost/invalid_handle.hh>
#include <ost/mol/impl/connector_impl.hh>
#include <ost/mol/impl/torsion_impl.hh>
#include <ost/mol/impl/dihedral.hh>

namespace ost { namespace mol {

TorsionHandle::TorsionHandle():
  impl_()
{}

TorsionHandle::TorsionHandle(const impl::TorsionImplP& im):
  impl_(im)
{}

TorsionHandle::operator bool() const
{
  return impl_.get()!=0;
}

bool TorsionHandle::IsValid() const {
  return impl_.get()!=0;
}

void TorsionHandle::Apply(EntityVisitor& v)
{
  this->CheckValidity();  
  impl_->Apply(v);
}

const String& TorsionHandle::GetName() const 
{
  this->CheckValidity();  
  return impl_->GetName();
}

AtomHandle TorsionHandle::GetFirst() const
{
  this->CheckValidity();  
  return impl_->GetFirst();
}

AtomHandle TorsionHandle::GetSecond() const
{
  this->CheckValidity();  
  return impl_->GetSecond();
}

AtomHandle TorsionHandle::GetThird() const
{
  this->CheckValidity();  
  return impl_->GetThird();
}

AtomHandle TorsionHandle::GetFourth() const
{
  this->CheckValidity();  
  return impl_->GetFourth();
}

geom::Vec3 TorsionHandle::GetPos() const
{
  this->CheckValidity();  
  return impl_->GetPos();
}

geom::Vec3 TorsionHandle::GetOriginalPos() const
{
  this->CheckValidity();
  return impl_->GetOriginalPos();
}

Real TorsionHandle::GetAngle() const
{
  this->CheckValidity();  
  return impl_->GetAngle();
}

void TorsionHandle::CheckValidity() const {
  if (!impl_)
    throw InvalidHandle();
}

Real DihedralAngle(const AtomHandle& a1, 
                     const AtomHandle& a2,
                     const AtomHandle& a3, 
                     const AtomHandle& a4) {
  impl::Dihedral d(a1.Impl(), a2.Impl(), a3.Impl(), a4.Impl());
  return d.GetAngleICS();
}

std::ostream& operator<<(std::ostream& os, const TorsionHandle& torsion)
{
  os << torsion.GetName() << " (" << torsion.GetFirst() << ", "
     << torsion.GetSecond() << ", " << torsion.GetThird() << ", "
     << torsion.GetFourth() << ")";
  return os;
}


}}  // ns
